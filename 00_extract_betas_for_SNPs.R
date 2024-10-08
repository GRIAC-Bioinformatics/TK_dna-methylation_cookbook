## This script will extract data for probes from EPIC array that were developed for known common SNPs (their id starts with rs)
## You can use this data to check for sample swaps if there is also genotype data available for the same samples.

library(minfi)

baseDir <- "/folder/with/idats/"
head(list.files(baseDir, recursive = TRUE))

#Read sample sheet 

targets <- read.metharray.sheet(
  "/folder/with/idats/",
  pattern="metafile.csv",
  recursive = T,
  verbose = TRUE
)

RGset <- read.metharray.exp(targets = targets) #upload data to RGset

snp_betas_raw <- getSnpBeta(RGset) #extract betas for snps (EPIC v2 contains 63 SNPs)

write.csv(snp_betas_raw, "/folder/to/save/snp_betas_raw.csv", quote = F)

