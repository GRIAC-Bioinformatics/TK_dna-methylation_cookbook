# =============================================================================
# check.my.input.R
#
# Author: Vartika Bisht
# Created: 13-Aug-2025
#
# Description:
# This script validates input parameters and files required for a DNA 
# methylation data analysis workflow. It performs a series of checks to 
# ensure that the provided inputs meet expected formats and exist in the 
# specified locations before proceeding with downstream analysis.
#
# Checks performed:
#    --platform        → Must be listed in manifest key file
#    --assembly        → Must be listed in manifest key file
#    --manifestkey     → Manifest key file ( defult : data/manifest.annotation.key.csv)
#    --idat_dir        → Directory must exist
#    --metadatasheet   → File must exist and contain a Basename, Sample_ID, Plate, Well, Gender, Array, Slide
#    --outdir          → Directory must exist
#    If both --idat_dir and --metadatasheet are provided:
#        → Verify that all IDAT files referenced in the metadata sheet exist 
#          in the specified directory (both red and green channel files).
#
# Usage:
#   Rscript check.my.input.R --platform <PLATFORM> --assembly <ASSEMBLY> \
#       --idat_dir <DIR> --metadatasheet <FILE> --outdir <DIR>
#
# Notes:
#   - This script stops execution immediately if any check fails.
#   - Intended to be run before launching a full methylation analysis to avoid 
#     runtime errors due to missing or invalid input files.
# =============================================================================

library(optparse)

# Set up command line options
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
    make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help="Output directory for the results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# --- Check for required parameters (default=NULL) ---
message("Checking for required parameters...")
params <- c()
checks <- c(
    "platform" = opt$platform,
    "assembly" = opt$assembly,
    "idat_dir" = opt$idat_dir,
    "metadatasheet" = opt$metadatasheet,
    "outdir" = opt$outdir
)

# Check --platform
if (!is.null(opt$platform)) {
    params <- c(params, "--platform")
}

# Check --assembly
if (!is.null(opt$assembly)) {
    params <- c(params, "--assembly")
}

# Check --idat_dir
if (!is.null(opt$idat_dir)) {
    params <- c(params, "--idat_dir")
}   

# Check --metadatasheet
if (!is.null(opt$metadatasheet)) {
    params <- c(params, "--metadatasheet")
}

# Check --outdir
if (!is.null(opt$outdir)) {
    params <- c(params, "--outdir")
}

if(length(params) > 0){
    message("The following parameters are provided: ", paste(params, collapse = ", "))
} else {
    message("No parameters are provided.")
    message("Please provide at least one parameter to proceed.")
    stop("Exiting the script due to missing parameters.")
}


check <- function(check_name, check_value){
    
    # https://github.com/achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38.git
    # https://github.com/jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38.git
    # https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylationEPICanno.ilm10b4.hg19.html
    # https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html

    # Platform checks
    if(check_name == "platform"){
        message("Checking platform compatibility...")
        message("Opening ",opt$manifestkey," to check platform compatibility...")
        tryCatch({
            manifest <- read.csv(opt$manifestkey, stringsAsFactors = FALSE)
            if(!(check_value %in% unique(manifest$platform))){
                stop(paste("Platform", check_value, "is not supported. Please check the manifest file. Please use one of the following platforms:", 
                           paste(unique(manifest$platform), collapse = ", ")))
            }
        }, error = function(e) {
            stop("Error reading the manifest file: ", e$message)
        })
        message("Platform is compatible.")
    }

    # Assembly checks
    if(check_name == "assembly"){
        message("Checking assembly version compatibility...")
        message("Opening ",opt$manifestkey," to check platform compatibility...")
        tryCatch({
            manifest <- read.csv(opt$manifestkey, stringsAsFactors = FALSE)
            if(!(check_value %in% unique(manifest$assembly))){
                stop(paste("Assembly", check_value, "is not supported. Please check the manifest file. Please use one of the following platforms:", 
                           paste(unique(manifest$assembly), collapse = ", ")))
            }
        }, error = function(e) {
            stop("Error reading the manifest file: ", e$message)
        })
        message("Assembly version is compatible.")
    }

    # IDAT directory checks
    if(check_name == "idat_dir"){
        message("Checking IDAT directory...")
        if(!dir.exists(check_value)){
            stop(paste("IDAT directory does not exist:", check_value))
        }
        message("IDAT directory exists.")
    }

    # Metadata sheet checks
    if(check_name == "metadatasheet"){
        message("Checking metadata sheet...")
        if(!file.exists(check_value)){
            stop(paste("Metadata sheet file does not exist:", check_value))
        }
        message("Metadata sheet file exists. Reading the file...")
        metadata <- read.csv(check_value, stringsAsFactors = FALSE)
        required_columns <- c("Basename", "Sample_ID","Plate", "Well", "Gender", "Array", "Slide")
        for(col in required_columns){
            if(!col %in% colnames(metadata)){
                stop(paste("Metadata sheet must contain a", col, "column. With exact spelling and case."))
            }
            message(paste("Metadata sheet has", col, "column."))
        }
        message("Metadata sheet is valid.")
    }

    # Output directory checks
    if(check_name == "outdir"){
        message("Checking output directory...")
        if(!dir.exists(check_value)){
            stop(paste("Output directory does not exist:", check_value))
        }
        message("Output directory exists.")
    }

}

for(p in names(checks)){
    check(p, checks[[p]])
}

if( "idat_dir" %in% names(checks) && "metadatasheet" %in% names(checks) ){

    message("As both metadata sheet and IDAT directory are provided, checking if IDAT files exist in the directory...")
    metadata <- read.csv(opt$metadatasheet, stringsAsFactors = FALSE)
    idat_files <- file.path(opt$idat_dir,metadata$Basename)
    check_result <- unlist(lapply(idat_files, function(x) as.numeric(system(paste0("ls ", x, "*.idat | wc -l"), intern = TRUE)) == 2 ))
    if(all(check_result)){
        message("All IDAT files exist in the directory specified.")
    } else {
        missing_files <- idat_files[!check_result]
        stop("The following IDAT files are missing: ", paste(missing_files, collapse = ", "))

    } 
}   else {
    message("IDAT directory and metadata sheet are required to check if IDAT files exist.")    
    }