# Title: Sample QC for Methylation Analysis
# Date: 18-09-2024

# Load required packages
library(minfi)
library(tidyverse)
library(reshape2)
library(Polychrome)

# Set working directory (adjust as needed)
setwd("/path/to/your/working/directory")

# Load prepared data
load("rgSet.RData")
load("grSet.RData")
load("MSet.RData")

## Qc metrics for sample quality check
# - Plot log median (met) vs log median(un-met).
# - Detection p-value for every CpG site in all samples.
# - Overall distributions of Beta values for each sample.
# - Check internal control probes:
#   - hybridization
#   - bisulfite conversion(only this one with minfi?)
#   - staining control
#   - extension control


## Plot log median (met) vs log median(un-met)
#Visual inspection of all samples. It is expected to see in the median plot good samples clustering together and 
#the bad QC samples separating from the good quality samples due to lower median values.
#To retrieve QC information from MSet object is used getQC() and plotQC() for visualization. 

# **Note** The threshold to decide which samples are good or bad is totally arbitrary and might change depending on the dataset.

qc <- getQC(MSet)
png("plot_log_median.png")
plotQC(qc)
dev.off()

## Mean detection p-values
# The function detectionP() is a method created by minfi to calculate the detection
#   p-values by taking the total signal (methylated+unmethylated) for each probe
#   to the background signal level. This background signal level is estimated from
#   the negative control probes. Small p-values are indicative of a reliable signal
#   whilst large p-values (p > 0.01) generally indicate a poor quality signal.
# **Note** p-value cut-off is an arbitrary choice

detP <- detectionP(rgSet)
df_man <- melt(detP, varnames = c("CpGID", "sample"), value.name = "pvalue")

png("barplot_detection_pvalue.png")
barplot(
  colMeans(detP),
  las = 2,
  cex.names = 0.8,
  ylab = "Mean detection p-values"
)
abline(h = 0.05, col = "red")
dev.off()


## Overall distributions of Beta values for each sample
# Here we expect beta values presenting values close to zero or one. 
# This is an indication of unmethylatilated CpG sites when close to zero and methylated when close to one.

png("densityPlot.png",
    width = 3000, height = 2000)
densityPlot(
  MSet,
  legend = FALSE
)
dev.off()

## Multidimensional scaling (MDS) plot
# The MDS plot can be use to identify sample mix and major sources in variation.
# You can check which technical/biological factors are contributing to the major source of variation. Sex could be the main one. 

group2 <- MSet@colData$Gender
names <- MSet@colData$Sample_Name

png("mds_plot_sex.png", width = 800, height = 600)
mdsPlot(
  MSet,
  sampGroups = group2,
  sampNames = names,
  legendNCol = 2,
  legendPos = "bottom"
)
dev.off()

## You can also check for sex missmatch 
# Sex is determined by the getSex function using copy number information from X and Y chr
# Estimation of sex is based on the median values of measurements on
# the X and Y chromosomes respectively. If ‘yMed’ - ‘xMed’ is less
# than ‘cutoff’ we predict a female, otherwise male.

object = getSex(grSet) 
metadata <- metadata[match(grSet@colData$Sample_Name, metadata$Sample_Name), ]

png("getSex_plot.png",
    width = 1000, height = 1000)
plot(x = object$xMed, y = object$yMed, type = "n", xlab = "X chr, median total intensity(log2)", ylab = "Y chr, median total intensity(log2)")
text(x = object$xMed, y = object$yMed, labels = as.character(rownames(pData(grSet))), 
     col = dplyr::case_when(
       metadata$GENDER == "M" ~ "deepskyblue",
       metadata$GENDER == "F" ~ "deeppink3",
       metadata$GENDER == NA ~ "grey"))
legend("bottomleft", c("M", "F"), col=c("deepskyblue", "deeppink3"), pch=16)
dev.off()

## Check control probes
# The plots of the different control probes can be exported into a pdf file in one step using the function qcReport

png("control_probes_bisulfite_conversion.png", width = 800, height = 600)
controlStripPlot(rgSet, controls = c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II"))
dev.off()

# Generate QC report
qcReport(rgSet, pdf = "qcReport.pdf")

# Remove poor quality samples: create a list of poor quality sampels:
to.rm <- c("")

# Remove from rgSet
keep <- which(!colnames(rgSet) %in% to.rm)
sf.rgSet <- rgSet[, keep]
save(sf.rgSet, file = "sf_rgSet.RData")

# Remove from MSet
keep <- which(!colnames(MSet) %in% to.rm)
sf.MSet <- MSet[, keep]
save(sf.MSet, file = "sf_MSet.RData")

# Remove from grSet
keep <- which(!colnames(grSet) %in% to.rm)
sf.grSet <- grSet[, keep]
save(sf.grSet, file = "sf_grSet.RData")

# Print summary of filtered objects
print(sf.rgSet)
print(sf.MSet)
print(sf.grSet)

# References
# - https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html
# - doi: 10.1186/s13059-015-0600-x