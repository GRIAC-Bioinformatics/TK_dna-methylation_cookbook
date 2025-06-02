# sex mismatch sample

## You can also check for sex missmatch 
# Sex is determined by the getSex function using copy number information from X and Y chr
# Estimation of sex is based on the median values of measurements on
# If a region of the genome has a higher copy number, the total intensity will be higher, and vice versa.
# Expected CN Differences by Sex:
# Females (XX): Will generally have higher copy number values (due to two copies) on the X chromosome compared to males.
# Males (XY): Will generally have lower copy number values on the X chromosome (one copy) and detectable copy number values on the Y chromosome (one copy), whereas females should have near-zero or very low copy number values on the Y chromosome.
# The difference between the median copy number values of the Y and X chromosomes is used to determine the sex of the sample.
# Males (XY): yMed will be relatively high compared to XX (due to the presence of a Y chromosome) and xMed will be relatively low (due to one X chromosome) compared to XX. Therefore, yMed - xMed will tend to be a positive and relatively high value.
# Females (XX): yMed will be very low (close to zero, as there's no Y chromosome) and xMed will be relatively high (due to two X chromosomes). Therefore, yMed - xMed will tend to be a negative and relatively low value.
# Hence, we take min of difference between the median copy number values of the Y and X chromosomes (yMed - xMed) to determine the center for female and the max for the center for males
# the X and Y chromosomes respectively. If ‘yMed’ - ‘xMed’ is less
# than ‘cutoff’ we predict a female, otherwise male.
# Rscript 02_sample_sex_mismatch.R -g grSet.RData -c Gender -f Flagged_sample_sex_mismatch.csv -p sample_sex_mismatch.pdf


# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(ggplot2)
  library(dplyr)
  library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Set up command line options
option_list <- list(
  make_option(c("-g", "--grset"), type="character", default=NULL,
              help="Path to GRset file", metavar="FILE"),
  make_option(c("-c", "--gendercolname"), type="character", default=NULL,
              help="Column name for gender", metavar="STRING"),
  make_option(c("-f", "--flagged"), type="character", default=getwd(),
              help="CSV file with flagged samples where sex mismatch", metavar="DIR"),
  make_option(c("-p", "--pdf"), type="character", default=getwd(),
              help="PDF file with all plots", metavar="DIR")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$grset) || !file.exists(opt$grset)) {
  stop("Please provide a valid path to the GRset file using -g or --grset option.")
}
if (is.null(opt$gendercolname)) {
  stop("Please provide the column name for gender using -c or --gendercolname option.")
}
if (is.null(opt$flagged)) {
  stop("Please provide the path to save flagged samples using -f or --flagged option.")
}
if (is.null(opt$pdf)) {
  stop("Please provide the path to save the PDF output using -p or --pdf option.")
}

message("Checking sex mismatch...")

load(opt$grset)

# Observed sex is taken from the column in the colData of the GRset object
observedSex <- grSet@colData[[opt$gendercolname]]
names(observedSex) <- rownames(grSet@colData)

# Predicted Sex is determined using the copy number values of X and Y chromosomes
# Get copy number values which are defined as the sum of the methylation and unmethylation channel.
CN <- getCN(grSet)

# Calculate median values for X and Y chromosomes
xIndex <- which(seqnames(grSet) == "chrX")
yIndex <- which(seqnames(grSet) == "chrY")
xMed <- colMedians(CN, rows = xIndex, na.rm = TRUE, useNames = TRUE)
yMed <- colMedians(CN, rows = yIndex, na.rm = TRUE, useNames = TRUE)

# Calculate the difference between Y and X medians
dd <- yMed - xMed

# Remove samples with NaNs which can occur if there are no probes on the X or Y chromosome
# This can happen if the sample are controls.
message("Removing following samples as they have NaN values...")
print(dd[is.na(dd)])
dd <- dd[!is.na(dd)]
dd <- data.frame((dd))
colnames(dd) <- "CopyNumberDifference"
dd <- cbind(dd, observedSex = observedSex[rownames(dd)])

# remove samples with no sex information
dd <- dd[dd$observedSex != "-", ]

# Find local maxima to determine centers for clustering
message("Finding local maxima to determine centers for clustering...")


local_max_indices <- c()  

for (i in seq_along(density(dd$CopyNumberDifference)$y)) {
  if (i == 1 || i == length(density(dd$CopyNumberDifference)$y)) {
    next  # Skip the first and last elements
  }
  if (density(dd$CopyNumberDifference)$y[i] > density(dd$CopyNumberDifference)$y[i - 1] && 
      density(dd$CopyNumberDifference)$y[i] > density(dd$CopyNumberDifference)$y[i + 1]) {
    local_max_indices <- c(local_max_indices, i)
  }
}

female_center <- min(density(dd$CopyNumberDifference)$x[local_max_indices])
male_center <- max(density(dd$CopyNumberDifference)$x[local_max_indices])

message("Local maxima indices found: ", paste(density(dd$CopyNumberDifference)$x[local_max_indices], collapse = ", "))
message("Female center: ", female_center)
message("Male center: ", male_center)

# Check variation and run clustering
if (length(unique(dd$CopyNumberDifference)) > 1) {
  message("Clustering samples based on copy number values...")
  message("Number of samples with valid copy number values: ", length(dd$CopyNumberDifference))
  message("Centers for clustering will be set to ", female_center, " for females and ", male_center, " for males.")
  # Clusetring using k-means
  k <- kmeans(dd$CopyNumberDifference, c(female_center,male_center), iter.max = 10, nstart = 1, algorithm = "Hartigan-Wong")
  # Determine which cluster - min for females (1) - max for males (2)
  predictedSex <- k$cluster
  predictedSex <- ifelse(k$cluster == 1 , "Female","Male")
} else {
  warning("Not enough variation to estimate sex")
  predictedSex <- rep("Unknown", length(dd$CopyNumberDifference))
}

dd <- cbind(dd, predictedSex)
dd$MatchStatus <- ifelse(dd$observedSex == dd$predictedSex, "Match", "Mismatch")


pdf(opt$pdf, width = 6, height = 6)

ggplot(dd, aes(x = CopyNumberDifference)) +
  geom_density(fill = "#1E90FF", color = "#333333", alpha = 0.8, linewidth = 0.8) + 
  labs(title = "Distribution of Sex Chromosome\nCopy Number Differences", 
       x = expression(bold("Copy Number Difference (Y - X Chromosome)")),
       y = expression(bold("Density"))) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.5), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    legend.position = "bottom"
  ) +
  geom_vline(xintercept=female_center, linetype="dashed", color="black", linewidth=1) +
  geom_vline(xintercept=male_center, linetype="dashed", color="black", linewidth=1)


ggplot(dd, aes(x = CopyNumberDifference, fill = MatchStatus)) +
  geom_histogram( alpha = 0.7, binwidth=0.05) + 
  labs(title = "Distribution of Sex Chromosome\nCopy Number Differences", 
       x = expression(bold("Copy Number Difference (Y - X Chromosome)")),
       y = expression(bold("Frequency"))) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.5), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("Match" = "#90EE90", "Mismatch" = "#FF0000")) + 
  geom_vline(xintercept=female_center, linetype="dashed", color="black", linewidth=1) +
  geom_vline(xintercept=male_center, linetype="dashed", color="black", linewidth=1)


beta <- getBeta(grSet)
  
# Handle missing/infinite values
if (any(is.na(beta))) {
  message("NA values found - imputing with row medians")
  beta_imputed <- apply(beta, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
  beta <- t(beta_imputed)
}

if (any(!is.finite(beta))) {
  message("Infinite values found - replacing with row max/min")
  beta_fixed <- apply(beta, 1, function(x) {
    x[x == Inf] <- max(x[is.finite(x)], na.rm = TRUE)
    x[x == -Inf] <- min(x[is.finite(x)], na.rm = TRUE)
    return(x)
  })
  beta <- t(beta_fixed)
}

# Transpose and perform PCA
beta_t <- t(beta)
pca_result <- prcomp(beta_t, scale. = TRUE, center = TRUE)

# Calculate variance explained
pca_var <- pca_result$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

# Prepare plot data
plot_data_pca <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)
rownames(plot_data_pca) <- rownames(pca_result$x)

plot_data_pca$MatchStatus <- "Match"
plot_data_pca[rownames(dd[dd$MatchStatus == "Mismatch",]), "MatchStatus"] <- "Mismatch"
ggplot(plot_data_pca, aes(x = PC1, y = PC2), color = "gray50") +
  geom_point(alpha = 0.7) +
  geom_point(data=plot_data_pca[plot_data_pca$MatchStatus == "Mismatch",], aes(x = PC1, y = PC2), color = "red") +
  theme(legend.position = "right") +
  labs( title = "PCA with all samples",
        subtitle = "Samples with sex mismatch are highlighted",
        x = paste0("PC1 (", pca_var_per[1], "%)"),
        y = paste0("PC2 (", pca_var_per[2], "%)")) +
theme_minimal() +
theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dev.off()


write.csv(data.frame("sexmismatch" = rownames(dd[dd$MatchStatus == "Mismatch",])), file = opt$flagged, row.names = TRUE)