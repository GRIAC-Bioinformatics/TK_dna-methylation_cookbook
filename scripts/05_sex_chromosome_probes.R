#!/usr/bin/env Rscript
# =============================================================================
# 05_sex_chromosome_probes.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script generates a CSV file containing all probes located on chrX and chrY
# from a given GenomicRatioSet (GRset) object.
# Arguments:
# --grset : Path to the GRset file (RDS format)
# --flagged : Output CSV file path for saving list of probes associated with ChrX and ChrY
# Usage:
#   Rscript 05_sex_chromosome_probes.R --grset <FILE.RData> \
#     --flagged <sex_mismatch_failures.csv> 
# =============================================================================


suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# Set up command line options
option_list <- list(
  make_option(c("-g", "--grset"), type="character", default=NULL,
              help="Path to GRset file", metavar="FILE"),
  make_option(c("-f", "--flagged"), type="character", default=getwd(),
              help="CSV file with probes associated to sex chromosome", metavar="DIR")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("grset", "flagged" )

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}

load(opt$grset)  

flagged = data.frame("sexprobes" = c(rownames(grSet)[which(seqnames(grSet) == "chrX")] , rownames(grSet)[which(seqnames(grSet) == "chrY")]))

write.csv(flagged, file = opt$flagged, row.names = FALSE)

message("All sex probes are saved at ", opt$flagged)