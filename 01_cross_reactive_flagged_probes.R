#!/usr/bin/env Rscript
# Rscript 01_cross_reactive_flagged_probes.R --crossreactprobes 48639-non-specific-probes-Illumina450k.csv --multimapprobes HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt --flagprobe cross_reactive_flagged_probes.csv
# ------------------------------------------------------------------------------
# Remove cross-reactive and unreliable probes
# Background:
# - Approximately 6% of array probes may generate spurious signals due to co-hybridization
# - These probes target sequences with high homology to unintended genomic regions
# - Primary reference: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1

# ------------------------------------------------------------------------------
# Platform-Specific Filtering Resources
# ------------------------------------------------------------------------------

# (1) For EPIC V2 array:
# - Reference: Peters et al. 2024 (https://doi.org/10.1186/s12864-024-10027-5)
# - Contains custom manifest with identified off-target CpG sites
# - Cross-reactive probes available in Additional file 8 from the publication

# (2) For 450k array:
# - Comprehensive filtering resources available at:
#   https://github.com/sirselim/illumina450k_filtering
#
# Key 450k filtering files:
# - Chen et al. 2016 cross-reactive probes: "cross_reactive_probes_file.csv"
# - BOWTIE2 multi-mapping probes (hg19): "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt"
#   - 33,457 probes with potential hybridization issues
# - Chen et al. 2013 non-specific probes: "48639-non-specific-probes-Illumina450k.csv"
#   - 29,233 probes identified
#   - Note: Overlap exists between cross-reactive and non-specific probe sets

# (3) For 850k (EPIC) array:
# - Reference: Pidsley et al. 2016 (https://doi.org/10.1186/s13059-016-1066-1)
# - Supplementary data suggests filtering cross-reactive and variant-containing probes
# - Note: Some overlap with 450k probe lists (no adverse effects expected)

# -------------------------------------------------------------------------------------------------
# Feature               | Cross-Reactive Probes               | Multi-Mapping Probes
# ----------------------|-------------------------------------|------------------------------------
# Cause                 | Sequence similarity (homology)      | Short probe length + repeats
# Detection             | Wet-lab hybridization               | Bioinformatics alignment
# Example               | Binds Gene A and its pseudogene     | Aligns to Chr1:1000 and Chr5:2000
# -------------------------------------------------------------------------------------------------
# Some cross-reactive probes are multi-mapping, but not all vice versa
# ------------------------------------------------------------------------------

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
library(minfi)
library(ggplot2)
library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# Parse command line arguments 
option_list <- list(
  make_option(c("--crossreactprobes"), type="character",
              help="non-specific probes across the 450k design identified by Chen et al.", metavar="character"),
  make_option(c("--multimapprobes"), type="character",
              help="All probe sequences were mapped to the human genome (hg19) using BOWTIE2 to identify potential hybridisation issues. 33,457 probes were identified as aligning greater than once and these were made available in HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt. Check sirselim/illumina450k_filtering:master", metavar="character"),
  make_option(c("--flagprobe"), type="character",
              help="Cross reactive flagged probe", metavar="character")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

message("Loading cross-reactive probes...")
cross_reactive <- read.csv(file = opt$crossreactprobes, head = T, as.is = T)
cross_reactive_probes <- as.character(cross_reactive$TargetID)

# BOWTIE2 multi-mapped
print("Loading BOWTIE2 multi-mapping probes...")
multi_map <- read.csv(opt$multimapprobes, head = F, as.is = T)
multi_map_probes <- as.character(multi_map$V1)

# determine unique probes
failed.probes <- data.frame(unique(c(cross_reactive_probes, multi_map_probes)))

# Write output
write.csv(
  na.omit(failed.probes),
  file = opt$flagprobe,
  row.names = FALSE,
  col.names = FALSE
)

# Note for cross-reactive probes in EPIC V1 check: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1


