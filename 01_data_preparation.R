# Title: Data Preparation for Methylation Analysis
# Date: 18-09-2024

# Load required packages
library(minfi)
library(missMethyl) # Contains annot850K

# if you were using EPIC V2 array: 
library("IlluminaHumanMethylationEPICv2manifest") # https://github.com/jokergoo/IlluminaHumanMethylationEPICv2manifest

# Set working directory (adjust as needed)
setwd("/path/to/your/working/directory")

###### Read files ######

# Read intensity data (IDAT)
dataDirectory <- "path/to/idat/files/"
list.files(dataDirectory, recursive = TRUE)[1:5]  # Display first 5 files

# Read sample sheet
targets <- read.metharray.sheet(
  dataDirectory,
  pattern = "your_sample_sheet.csv",
  recursive = TRUE,
  verbose = TRUE
)
head(targets)

###### Add annotations ######
# The next step is to add annotations. Initial annotation infomration is available in illumina manifest files.
# Annotation contains most of the important information about the probes: location in the genome, probe types, etc
# Depending on the array type and genome build you want to use you need to adjust this part: 

ann850k <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # 850 array
head(ann850k[, 1:5])

# Read methylation data - Using the function **read.metharray.exp()** the raw intensity data is stored in the 
# RGChannelSet object that contains the intensities from the green and red channels. 

rgSet_ext <- read.metharray.exp(targets = targets, extended = TRUE) # extended version for wm::beadcount()
save(rgSet_ext, file = "rgSet_ext.RData")
rm(rgSet_ext) # remove extended object

rgSet <- read.metharray.exp(targets = targets)

# In case you are using EPIC V2 array:
# add array info and annotation name
rgSet@annotation <- c(array='IlluminaHumanMethylationEPICv2', annotation='20a1.hg38')

# Get an overview of the data
rgSet
# Phenotype data can be access pData
pData(rgSet)
# Overview of probe design
getManifest(rgSet)

# Check and modify sample names if needed
head(sampleNames(rgSet))

###### Preprocess data ######
# The preprocessRaw() function from Minfi package convert RGChannelSet to MethylSet, this is done by converting green and red channels into methylated and unmethylated matrices. 
# Important: note that no normalization is performed.
# RGChannelSet to MethylSet
MSet <- minfi::preprocessRaw(rgSet)
MSet

# Access unmethylated and methylated matrices
head(getMeth(MSet)[, 1:3])
head(getUnmeth(MSet)[, 1:3])

# Convert to RatioSet
# "A RatioSet object is class designed to store Beta and/or M-values instead of the (un)methylated signals. 
#  An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. 
#  Mapping a MethylSet to a RatioSet is irreversible, i.e. one cannot technically retrieve the methylated and unmethylated signals from a RatioSet. A RatioSet can be created with the function ratioConvert"
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
ratioSet

# Map to genome
# "The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation information. 
# The output object is a GenomicRatioSet"
grSet <- mapToGenome(ratioSet)
grSet

# Access Beta and M values
beta <- getBeta(grSet)
head(beta[1:4, 1:3])

m <- getM(grSet[1:4, 1:3])
head(m)

cn <- getCN(grSet)
head(cn[1:4, 1:3])

# Save objects
save(MSet, file = "MSet.RData")
save(rgSet, file = "rgSet.RData")
save(grSet, file = "grSet.RData")

# References
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html