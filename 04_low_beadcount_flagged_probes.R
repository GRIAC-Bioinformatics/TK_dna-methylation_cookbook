# Author: Vartika Bisht + Tatiana Karp
# Date: 24 May 2025
# Description: This script filters probes with low bead counts from an Extended RG Channel Set.
# Probes are flagged if more than a specified percentage of samples have a bead count below 3.
# Usage: Rscript 04_low_beadcount_flagged_probes.R -r RGsetEXT_sample_filtered.RData -o low_beadcount_flagged_probes.csv -p 5
# Input:
#   -r, --rgsetext: Path to the RGset object file (ExtendedRGChannelSet).
#
# Output:
#   -o, --output: Path to save the list of flagged probes (CSV file).
#
# Parameters:
#   -p, --percentage: Percentage threshold for filtering probes (default: 5).
#                     Probes with a percentage of samples having bead count < 3
#                     exceeding this threshold will be flagged.

suppressPackageStartupMessages({
  tryCatch({
    library(wateRmelon)
    library(minfi)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# Parse command line arguments
option_list <- list(
  make_option(c("-r", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGset extended object file"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save the list of failed probes"),
  make_option(c("-p", "--percentage"), type = "numeric", default = 5,
              help = "Percentage threshold for filtering probes (default: 5)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$rgsetext) || !file.exists(opt$rgsetext)) {
  stop("Please provide the path to the RGset object file using -r or --rgsetext option.")
}
if (is.null(opt$output)) {
  stop("Please provide the path to save the list of failed probes using -o or --output option.")
}

# Load the RGset object
load(opt$rgsetext)

## Perform bead count filtering
# Generate a matrix of bead counts. NAs represent probes with bead count < 3.
df.bead.counts <- wateRmelon::beadcount(RGset)

# Define a function to calculate the percentage of NAs in each row (probe)
pct.of.nas <- function(row) {
  sum(is.na(row)) / length(row) * 100
}

# Calculate the percentage of NAs for each probe
df.bead.counts$na_pct <- apply(df.bead.counts, 1, pct.of.nas)

# Identify probes that exceed the specified percentage of low bead counts
flagged_probes <- rownames(df.bead.counts[which(df.bead.counts$na_pct > opt$percentage),])

# Save the list of flagged probes to the specified output file
write.csv(flagged_probes, file = opt$output, quote = FALSE, row.names = FALSE, col.names = FALSE)

print(paste("Number of probes flagged due to low bead counts: ", length(flagged_probes)))
print(paste("List of flagged probes saved to: ", opt$output))