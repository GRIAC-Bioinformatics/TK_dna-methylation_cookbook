#!/usr/bin/env Rscript
# =============================================================================
# 01_create_rgset_object.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script processes raw IDAT files to generate various minfi objects
# (RGChannelSet, MethylChannelSet, GenomicRatioSet, RGChannelSetExtended).
# These objects are essential for downstream DNA methylation data analysis,
# providing structured data containers for further processing and quality control.
# Arguments:
#    --platform               → DNA methylation platform (e.g., Illumina HumanMethylation450)
#    --manifestkey            → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly               → Assembly version of the platform
#    --idat_dir               → Path to the directory containing IDAT files
#    --metadatasheet          → Path to the metadata sheet file
#    --rg_channel_set_output  → Output path for the RGChannelSet object (.RData)
#    --methyl_channel_set_output → Output path for the MethylChannelSet object (.RData)
#    --genomic_ratio_set_output → Output path for the GenomicRatioSet object (.RData)
#    --rg_channel_set_extended_output → Output path for the RGChannelSetExtended object (.RData)
# Usage:
#   Rscript 01_create_rgset_object.R --platform <PLATFORM> --assembly <ASSEMBLY> \
#       --idat_dir <DIR> --metadatasheet <FILE> \
#       --rg_channel_set_output <FILE.RData>
# Notes:
#   - This script requires a manifest key and metadata sheet to correctly process
#     and annotate the IDAT files.
#   - The output files are saved in RData format for easy loading into subsequent scripts.
# =============================================================================

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

option_list <- list(
  make_option(c("-p", "--platform"), type="character", default=NULL, 
              help="DNA methylation platform", metavar="character"),
  make_option(c("-k", "--manifestkey"), type="character", default="data/manifest.annotation.key.csv", 
              help="Manifest key file", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default=NULL,
              help="Assembly version of the DNA methylation platform", metavar="character"),
  make_option(c("-i", "--idat_dir"), type = "character", default = NULL,
              help="Path to IDAT files directory", metavar="character"),
  make_option(c("-m", "--metadatasheet"), type = "character", default = NULL,
              help="Path to the metadata sheet file", metavar="character"),
  make_option(c("-R", "--rg_channel_set_output"),type = "character",
              help = "Output file path for RGChannelSet (.RData)",metavar = "FILE"),
  make_option(c("-M", "--methyl_channel_set_output"),type = "character",
              help = "Output file path for MethylChannelSet (.RData)",metavar = "FILE"),
  make_option(c("-G", "--genomic_ratio_set_output"),type = "character",
              help = "Output file path for GenomicRatioSet (.RData)",metavar = "FILE"),
  make_option(c("-E", "--rg_channel_set_extended_output"), type = "character",
              help = "Output file path for RGChannelSetExtended (.RData)", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("platform", "assembly", "idat_dir", "metadatasheet" ,"rg_channel_set_output", 
                    "methyl_channel_set_output", "genomic_ratio_set_output", "rg_channel_set_extended_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}

# LOAD MANIFEST AND ANNOTATION
message("Platform: ", opt$platform)
message("Assembly: ", opt$assembly)
message("Reading ",opt$manifestkey," ...")
manifest.key <- tryCatch({
  read.csv(opt$manifestkey, header = TRUE)
}, error = function(e) {
  stop("Failed to read manifest file: ", e$message)
})

Nmanifest <- manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$manifest
Nannotation <- manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$annotation

message("According to the ",opt$manifestkey," file, the following manifest and annotation will be used:")
message("Manifest: ", Nmanifest)
message("Annotation: ", Nannotation)

message("Loading the manifest and annotation packages...")
tryCatch({
    library(Nmanifest, character.only = TRUE)
    library(Nannotation, character.only = TRUE)
}, error = function(e) {
    stop("Failed to load manifest or annotation package: ", e$message)
})

message("Manifest and annotation packages loaded successfully.")

# META DATA FILE PROCESSING
message("Creating list of idats files using the filename column in ", opt$metadatasheet, " and the IDAT dir " ,opt$idat_dir)
message("Reading metadata from ", opt$metadatasheet, "...")
metadata <- tryCatch({
  read.csv(opt$metadatasheet, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
  stop("Failed to read metadata file: ", e$message)
})
metadata$Basename <- file.path(opt$idat_dir,metadata$Basename)
message("Found ", dim(metadata)[1], " IDAT pair files in the directory.")


# IDAT FILE PROCESSING
message("Reading IDAT files (this may take several minutes)...")
# read.metharray.exp() fully loads all IDATs into memory and constructs dense matrices.
# This is memory-intensive due to duplicated quant matrices and large final RGChannelSet.
# Expect RAM usage to exceed 2x the IDAT file size. For 22 samples (450k), >60GB is normal.
# List of Quant matrices → Combined Matrix:
# GreenMean <- do.call(cbind, lapply(G.Quants, ...))
# RedMean   <- do.call(cbind, lapply(R.Quants, ...))
# This constructs large in-memory double-precision numeric matrices.
# These matrices are dense and not compressed. For ~485,000 probes × 22 samples × 2 channels (Green + Red), this is a lot of RAM.
# https://rdrr.io/bioc/minfi/man/read.metharray.exp.html 
####################################################################################################################################################
# Instead of using read.metharray.exp(), we will use read.metharray.exp.par() which constructs RGChannelSet in parallel using BiocParallel.
# BiocParallel uses MulticoreParam on Unix-like systems (Linux/macOS), which distributes the workload across multiple CPU cores on a single machine.
# The function splits the IDAT files into chunks and processes them in parallel. Each core reads and parses a subset of the IDAT files independently. 
# After all files are read, the results are combined into a single RGChannelSet object.
####################################################################################################################################################
# Well, that did not work. the par version fo read.metharray.exp() always runs out of memory. If the files are too large, it will crash as you need
# quite some RAM to read the IDAT files. So , you do not want to do it parallel.

message("Reading IDAT files into RGChannelSet object...")

RGset <- tryCatch({
  read.metharray.exp(targets = metadata,
                     extended = FALSE, 
                     verbose = TRUE)
}, error = function(e) {
  stop("IDAT processing failed : ", e$message)
})

rownames(RGset@colData) = RGset@colData$Sample_ID
annotation(RGset)["array"] = manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$name_manifest
annotation(RGset)["annotation"] = manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$name_annotation

message("RGChannelSet object created. Now, saving RGChannelSet object...")

tryCatch({
  save(RGset, file = opt$rg_channel_set_output)
}, error = function(e) {
  stop("Failed to save RGChannelSet object output file : ", e$message)
})
message("Saved RGChannelSet object.")

message("Converting RGChannelSet object into MethylChannelSet object...")

tryCatch({
  MSet <- minfi::preprocessRaw(RGset)
}, error = function(e) {
  stop("Failed to convert RGChannelSet object into MethylChannelSet object : ", e$message)
})

message("MethylChannelSet object created. Now, saving MethylChannelSet object...")

tryCatch({
  save(MSet, file = opt$methyl_channel_set_output)
}, error = function(e) {
  stop("Failed to save MethylChannelSet object output file: ", e$message)
})

message("Saved MethylChannelSet object.")

message("Removing RGChannelSet object to make judicial use of the space.")
rm(RGset)

message("Converting MethylChannelSet object to GenomicRatioSet object...")

tryCatch({
  ratioSet <- minfi::ratioConvert(MSet, what = "both", keepCN = TRUE)
  grSet <- minfi::mapToGenome(ratioSet)
}, error = function(e) {
  stop("Failed to convert MethylChannelSet object to GenomicRatioSet object: ", e$message)
})

message("GenomicRatioSet object created. Now, saving GenomicRatioSet object...")

tryCatch({
  save(grSet, file = opt$genomic_ratio_set_output)
}, error = function(e) {
  stop("Failed to save GenomicRatioSet object output file: ", e$message)
})

message("Saved GenomicRatioSet object.")

message("Removing GenomicRatioSet object to make judicial use of the space.")
rm(ratioSet,grSet)

message("Reading IDAT files into RGChannelSetEXT object...")

RGsetEXT <- tryCatch({
  read.metharray.exp(targets = metadata,
                     extended = TRUE, 
                     verbose = TRUE)
}, error = function(e) {
  stop("IDAT processing failed : ", e$message)
})

rownames(RGsetEXT@colData) = RGsetEXT@colData$Sample_ID
annotation(RGsetEXT)["array"] = manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$name_manifest
annotation(RGsetEXT)["annotation"] = manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$name_annotation

message("RGChannelSetEXT object created. Now, saving RGChannelSetEXT object...")

tryCatch({
  save(RGsetEXT, file = opt$rg_channel_set_extended_output)
}, error = function(e) {
  stop("Failed to save RGChannelSetEXT object output file: ", e$message)
})
message("Saved RGChannelSetEXT object.")
