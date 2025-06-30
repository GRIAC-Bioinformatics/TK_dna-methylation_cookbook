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
# Rscript 04_high_intensity_flagged_probes.R -m MSet.RData -o high_intensity_flagged_probes.csv -c 10000
# Input:
#   -m, --mset: Path to the MSet object file. This file is expected to contain
#                both an 'MSet' object (for manifest information) and an 'MSet'
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
    library(ggplot2)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Parse command line arguments
option_list <- list(
  make_option(c("-m", "--mset"), type = "character", default = NULL,
              help = "Path to the RGset object file"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save the list of failed probes"),
  make_option(c("-c", "--high_intensity_cutoff"), type = "numeric", default = 10000,
              help = "High intensity cutoff threshold for filtering probes (default: 10000)"),
  make_option(c("-p", "--pdf"), type="character", default=NULL,
                help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$mset) || !file.exists(opt$mset)) {
  stop("Please provide the path to the MSet object file using -r or --mset option.")
}
if (is.null(opt$output)) {
  stop("Please provide the path to save the list of failed probes using -o or --output option.")
}

# Load the RGset object. This file is assumed to also contain the MSet object.
load(opt$mset)

# Check if the MSet object is loaded
if (!exists("MSet")) {
  stop("The 'MSet' (MethylSet) object was not found in the provided file. It is required for intensity extraction.")
}

# ---
## Extract probe information and intensities

# Get the manifest to identify Type I probes
manifest <- getManifest(MSet)
ProbeI <- getProbeInfo(manifest, type = "I")$Name

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
flagged_probes <- hi.MU[hi.MU %in% ProbeI]

# Save the list of flagged probes to the specified output file
write.csv(data.frame("FlaggedProbes" = flagged_probes), file = opt$output, quote = FALSE, row.names = FALSE)

print(paste("Number of probes flagged due to high intensity: ", length(flagged_probes)))
print(paste("List of flagged probes saved to: ", opt$output))


# Open PDF for plotting
pdf(opt$pdf, width = 10, height = 7)

# geometric mean intensities for each probe
MUmed <- data.frame(Probe = names(MUmed), GM = MUmed)

message("Generating intesity plot for type I probes...")
ggplot(MUmed, aes(x = GM)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  labs(title = "Distribution of Geometric Mean Intensity for Type I Probes",
        subtitle = paste0("Cutoff: ", opt$high_intensity_cutoff, " (red line indicates threshold)\n",
                          "Probes with geometric mean intensity above this cutoff are flagged"),
        x = "Geometric Mean Intensity",
        y = "Number of Probes") +
  theme_minimal() + 
  geom_vline(xintercept = opt$high_intensity_cutoff, color = "red", linetype = "dashed") +
  annotate("text", x = opt$high_intensity_cutoff, y = max(table(MUmed$GM)), 
            label = paste0("Cutoff ",opt$high_intensity_cutoff), color = "red", vjust = -1)

dev.off()

message("Filtering probes based on intensity...")
