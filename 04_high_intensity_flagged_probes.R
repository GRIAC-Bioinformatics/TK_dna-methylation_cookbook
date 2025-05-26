# Author: Vartika Bisht + Tatiana Karp
# Date: 24 May 2025
# Description: From this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4239800/
#  "For example, we have looked at the relation between signal intensities and DNA methylation measurements 
# [b-values, defined as the ratio of the methylated signal over the total signal (methylated + unmethylated) 
# and have observed that probes displaying a high average intensity (i.e. a high average of the methylated 
# and unmethylated signals) are more prone than probes displaying a lower average intensity to provide DNA 
# methylation measurements inconsistent with measurements obtained with other approaches, such as BPS (bisulfite pyrosequencing).
# They have a tendency to provide values close to 0.5, independently of their true methylation state. 
# Of note, type II Infinium probes seem to be less prone to this phenomenon."  
# Usage:
# Rscript 04_high_intensity_flagged_probes.R -r RGset_sample_filtered.RData -o high_intensity_flagged_probes.csv -c 10000
# Input:
#   -r, --rgset: Path to the RGset object file. This file is expected to contain
#                both an 'RGset' object (for manifest information) and an 'MSet'
#                object (MethylSet, for methylated and unmethylated signal intensities).
#
# Output:
#   -o, --output: Path to save the list of flagged probes (CSV file).
#
# Parameters:
#   -c, --high_intensity_cutoff: Numeric threshold for filtering probes (default: 10000).
#                                Probes with a median intensity greater than this
#                                value will be flagged.

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Parse command line arguments
option_list <- list(
  make_option(c("-r", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGset object file"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save the list of failed probes"),
  make_option(c("-c", "--high_intensity_cutoff"), type = "numeric", default = 10000,
              help = "High intensity cutoff threshold for filtering probes (default: 10000)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$rgset) || !file.exists(opt$rgset)) {
  stop("Please provide the path to the RGset object file using -r or --rgset option.")
}
if (is.null(opt$output)) {
  stop("Please provide the path to save the list of failed probes using -o or --output option.")
}

# Load the RGset object. This file is assumed to also contain the MSet object.
load(opt$rgset)

# Ensure both RGset and MSet objects are available after loading
if (!exists("RGset")) {
  stop("The 'RGset' object was not found in the provided file.")
}
if (!exists("MSet")) {
  stop("The 'MSet' (MethylSet) object was not found in the provided file. It is required for intensity extraction.")
}

# ---
## Extract probe information and intensities

# Get the manifest to identify Type I probes
manifest <- getManifest(RGset)
snpProbesI <- getProbeInfo(manifest, type = "I")$Name

# Extract methylated and unmethylated signal intensities from the MethylSet
Meth <- getMeth(MSet)
Unmeth <- getUnmeth(MSet)

# Calculate the median methylated and unmethylated intensities for each probe
Methmed <- apply(Meth, 1, median)
Unmethmed <- apply(Unmeth, 1, median)

# Calculate the combined median intensity (geometric mean) for each probe
MUmed <- sqrt(Methmed * Unmethmed)

# Identify all probes with combined median intensity above the specified cutoff
hi.MU <- names(which(MUmed > opt$high_intensity_cutoff))

# ---
## Filter and save flagged probes

# Select only Type I probes from the high-intensity list
flagged_probes <- hi.MU[hi.MU %in% snpProbesI]

# Save the list of flagged probes to the specified output file
write.csv(flagged_probes, file = opt$output, quote = FALSE, row.names = FALSE, col.names = FALSE)

print(paste("Number of probes flagged due to high intensity: ", length(flagged_probes)))
print(paste("List of flagged probes saved to: ", opt$output))