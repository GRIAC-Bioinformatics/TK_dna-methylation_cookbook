# Author: Vartika Bisht
# Date: 24 June 2025
# Description: This R script, `07_get_betamatrix_norm_methods.R`, was created on June 25, 2025, by Vartika Bisht. Its primary purpose is to compare the results of various normalization methods applied to DNA methylation data, specifically focusing on how these methods impact the beta matrix when control samples are utilized for analysis. The script processes raw methylation data, applies different normalization techniques, and outputs normalized beta matrices for downstream analysis.

# Input Files
# The script requires the following input files:
# Male_Positive_Controls.csv: A CSV file containing a list of controls to be used for the analysis.
# RGset_sample_filtered_probe_filtered.RData: An RData file containing the `RGChannelSet` object, which is the raw methylation data after sample and probe filtering.
# Existing Normalized RData files (optional):
#   SWAN.RData
#   BMIQ.RData
#   PBC.RData
#   dyeCorrFalse_NOOB.RData
#   dyeCorrReference_NOOB.RData
#   dyeCorrSingle_NOOB.RData
#   dyeCorrSingle_NOOB_BMIQ.RData
# These files represent intermediate or final normalized data from specific methods that can be referenced or compared.
# probe_type_raw_controls.csv: A CSV file containing information about probe types and raw controls, used to facilitate normalization or quality control.

# Output Files
# The script generates the following output: RData file containing the final normalized beta matrix.
# Other RData files corresponding to the individual normalization methods if they are generated as part of the script's execution.

# Parameters
# The script accepts several command-line parameters to specify input and output files, as well as to enable different normalization methods:
# -c <file_path> --controls <file_path>: Specifies the path to the CSV file containing the list of controls (e.g., `Male_Positive_Controls.csv`).
# -r <file_path> --rgset <file_path>: Specifies the path to the RData file containing the raw `RGChannelSet` (e.g., `RGset_sample_filtered_probe_filtered.RData`).
# -W <file_path> --swan <file_path>: Specifies the path for the output RData file after **SWAN** normalization.
#   SWAN (Subset-quantile within array normalization): This within-array normalization method corrects for technical differences between Type I and Type II array designs. It matches the Beta-value distributions of Type I and Type II probes by applying a within-array quantile normalization separately for different subsets of probes (divided by CpG content). It takes an `RGChannelSet` or `MethylSet` as input and returns a `MethylSet`.
# -B <file_path> --bmiq <file_path>: Specifies the path for the output RData file after **BMIQ** normalization.
#   BMIQ (Beta-mixture quantile dilation): An intra-sample normalization procedure that corrects the bias of Type II probe values. It involves a 3-step process: (i) fitting a 3-state beta mixture model, (ii) transforming state-membership probabilities of Type II probes into quantiles of the Type I distribution, and (iii) a conformal transformation for hemi-methylated probes. It is often called internally by functions like `champ.norm()`.
# -P <file_path> --pbc <file_path>: Specifies the path for the output RData file after **PBC** normalization.
#   PBC (Peak-based correction): This method normalizes methylation data by rescaling the methylation levels of Type II probes to match the distributions of Type I probes. It achieves this by estimating the modes (peaks) of the methylation level distributions for both probe types and adjusting Type II probe values accordingly.
# -F <file_path> --dyecorrfalse_noob <file_path>: Specifies the path for the output RData file after **NOOB** normalization without dye correction.
#   NOOB (Normal-exponential convolution on out of band probes): This background correction method uses the normexp model (assuming observed intensities are a sum of normally distributed background noise and exponentially distributed true signal) on out-of-band (OOB) probe data. OOB data from Type I probes (e.g., red channel fluorescence for a green channel probe) is used to estimate and correct background noise. `preprocessNoob` takes an `RGChannelSet` and returns a `MethylSet`.
# -R <file_path> --dyecorrreference_noob <file_path>: Specifies the path for the output RData file after **NOOB** normalization with **reference dye correction**.
#   Dye Bias Correction (Reference): This approach corrects dye bias by selecting a "reference" sample from the batch (the sample with the Red-to-Green ratio closest to 1, calculated from background-corrected internal normalization control probes) and normalizes all other samples relative to it.
# -S <file_path> --dyecorrsingle_noob <file_path>: Specifies the path for the output RData file after **NOOB** normalization with **single sample dye correction**.
#   Dye Bias Correction (Single Sample): This approach performs dye bias correction on a sample-by-sample basis. It uses internal normalization control probes (e.g., NORM_C, NORM_G, NORM_A, NORM_T on 450k/EPIC arrays) to calculate an average Red-to-Green ratio for each sample, which reflects the relative efficiency of the dyes. This ratio is then used to correct methylation values.
# -N <file_path> --dyecorrsingle_noob_bmiq <file_path>: Specifies the path for the output RData file after **NOOB** with single dye correction and **BMIQ** normalization.
# -O <file_path> --output <file_path>: Specifies the path for the main output RData file containing the final normalized beta matrix (e.g., `pOOBAH_ssNoob_BMIQ_normalised_beta_matrix.RData`).
# -p <file_path> --probetype <file_path>: Specifies the path to the CSV file containing probe type raw controls (e.g., `probe_type_raw_controls.csv`).

# Usage
# To run the script, use the following command structure in your terminal:
# Rscript ${SCRIPTDIR}/07_get_betamatrix_norm_methods.R \
#   -c Male_Positive_Controls.csv \
#   -r RGset_sample_filtered_probe_filtered.RData \
#   -W SWAN.RData \
#   -B BMIQ.RData \
#   -P PBC.RData \
#   -F dyeCorrFalse_NOOB.RData \
#   -R dyeCorrReference_NOOB.RData \
#   -S dyeCorrSingle_NOOB.RData \
#   -N dyeCorrSingle_NOOB_BMIQ.RData \
#   -O pOOBAH_ssNoob_BMIQ_normalised_beta_matrix.RData \
#   -p probe_type_raw_controls.csv
#
# Replace `${SCRIPTDIR}` with the actual directory where your script is located, and adjust file paths as necessary for your environment.


suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
    library(ggplot2)
    library(matrixStats)
    library(dplyr)
    library(tidyr)
    library(data.table)
    }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})



# Set up command line options
option_list <- list(
    make_option(c("-c", "--control_sample"), default=NULL,
                help="Path to file with list of control sample names", metavar="FILE"),
    make_option(c("-r", "--rgset"), default=NULL,
                help="Path to RGset RData file", metavar="FILE"),
    make_option(c("-p", "--probetype"), default="probetype.RData",
                help="Path to output file for probe type information", metavar="FILE"),
    make_option(c("-W", "--SWAN_normalised_MSet_Rdata"),  default=".",
                help="Path to MethylSet Object after SWAN normalisation", metavar="PATH"),
    make_option(c("-B", "--BMIQ_normalised_beta_matrix"),  default=".",
                help="Path to beta matrix after BMIQ normalisation", metavar="PATH"),
    make_option(c("-N", "--ssNoob_BMIQ_normalised_beta_matrix"),  default=".",
                help="Output beta matrix after ssNOOB + BMIQ normalisation", metavar="PATH"),
    make_option(c("-P", "--PBC_normalised_beta_matrix"),  default=".",
                help="Path to beta matrix after PBC normalisation", metavar="PATH"),
    make_option(c("-F", "--NOOB_dyeCorrFalse_MSet_Rdata"),  default=".",
                help="Path to MethylSet Object after NOOB normalisation", metavar="PATH"),
    make_option(c("-R", "--NOOB_dyeCorrRefernce_MSet_Rdata"),  default=".",
                help="Path to MethylSet Object after NOOB normalisation", metavar="PATH"),
    make_option(c("-O", "--pOOBAH_NOOB_dyeCorrSingle_MSet_Rdata"),  default=".",
                help="Path to MethylSet Object after NOOB normalisation", metavar="PATH"),
    make_option(c("-S", "--NOOB_dyeCorrSingle_MSet_Rdata"),  default=".",
                help="Path to MethylSet Object after NOOB normalisation", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$control_sample)) {
  stop("Control sample file is required. Use -c or --control_sample.")
}
if (is.null(opt$rgset)) {
  stop("RGset RData file is required. Use -r or --rgset.")
}

# Load control samples
controls = read.csv(opt$control_sample, header= TRUE)[,1]

# Get the beta values before normalisation ( control samples )
load(opt$rgset)
rawBetas <- getBeta(RGset[,controls])
save(rawBetas, file = file.path(dirname(opt$rgset), "raw_control_betas.RData"))

typeIprobes <- getProbeInfo(RGset[,controls], type = "I")$Name
typeIIprobes <- getProbeInfo(RGset[,controls], type = "II")$Name
probetype <- data.frame(
    probe = c(typeIprobes, typeIIprobes),
    type = c(rep("I", length(typeIprobes)), rep("II", length(typeIIprobes)))
)
save(probetype, file = opt$probetype)

# Remove variables not in use to make judicial use of space
rm(RGset)

# SWAN

# Errors with useNames = NA handeled with :
# https://github.com/satijalab/seurat/issues/7501#issuecomment-1854571904

# SWAN
# Perform Subset-quantile within array normalization (SWAN), a within-array normalization correction for the technical differences between the Type I and Type II array designs. The algorithm matches the Beta-value distributions of the Type I and Type II probes by applying a within-array quantile normalization separately for different subsets of probes (divided by CpG content). The input of SWAN is a MethylSet, and the function returns a MethylSet as well. If an RGChannelSet is provided instead, the function will first call preprocessRaw on the RGChannelSet, and then apply the SWAN normalization.

# Input: RGChannelSet or MethylSet
# Output: MethylSet

# MSet, an optional object of class MethylSet. If set to NULL preprocessSwan uses preprocessRaw on the rgSet argument. In case MSet is supplied, make sure it is the result of preprocessing the rgSet argument.
if (!is.null(opt$SWAN_normalised_MSet_Rdata)) {
    message("--- Starting SWAN normalization ---")
    load(opt$SWAN_normalised_MSet_Rdata)
    message("Loaded SWAN_normalised_MSet_Rdata.")
    SWAN_normalisedBetas <- getBeta(SWAN_MSet[,controls])
    message("Extracted beta values for SWAN.")
    rm(SWAN_MSet)
    message("Removed SWAN_MSet from memory.")
    save(SWAN_normalisedBetas, file = file.path(dirname(opt$SWAN_normalised_MSet_Rdata), "SWAN_normalised_control_betas.RData"))
}



# Noob
# Normal-exponential convolution on out of band probes (NOOB): The normal-exponential (normexp) convolution model was developed as part of the RMA algorithm for Affymetrix microarray data (see recommended reading section below for more details on this). The model assumes that the observed intensities are the sum of background noise and true signal components. The background is normally distributed and the signal is exponentially distributed. NOOB is simply using the normexp model on out of band probe data.
# Out-of-band (OOB) probe data is commonly used as part of the preprocessing pipleline. Recall that type I probes utilize only one color channel but have two probes for one locus. The probe that does not match the actual methylation state of the sample still captures a measure of fluorescence that can be used to help estimate background noise. In other words, if a type I probe is operating on the green channel to capture the methylation state of the locus, we will still have flourescence measures from the red channel (and that fluorescence from the red channel is “OOB” data).
# preprocessNoob is a function that implements the NOOB background correction method. It uses the OOB probe data to estimate the background noise and corrects the methylation values accordingly. The function can be applied to either an RGChannelSet or a MethylSet object, and it returns a MethylSet object with the corrected methylation values.
# Input: RGChannelSet
# Output: MethylSet

# How should dye bias correction be done: use a single sample approach (ssNoob), or a reference array?
# The main idea is to use internal control probes on the array, specifically the normalization controls (NORM_C, NORM_G for the green channel, and NORM_A, NORM_T for the red channel on 450k/EPIC arrays, or Normalization-Green and Normalization-Red for older array types). These probes are designed to have consistent, known signals, making them ideal for assessing and correcting dye bias.
# Background Correction of Controls: The intensities of the red and green internal control probes are first background corrected.
# Calculate Average Control Intensities: The average intensities of the relevant green and red control probes are calculated for each sample. Let's call these Green.avg and Red.avg.
# Determine Red-to-Green Ratio: The ratio R.G.ratio = Red.avg / Green.avg is computed for each sample. This ratio reflects the relative efficiency of the red dye compared to the green dye in that specific sample.
# Single --> sample by sample dye bias correction
# Reference --> selects a "reference" sample from the batch and normalizes all other samples relative to it.
# The reference sample is chosen as the sample with the R.G.ratio closest to 1.
# The output a MethylSet object with the corrected methylation values.

# NOOB No Dye Correction
if (!is.null(opt$NOOB_dyeCorrFalse_MSet_Rdata)) {
    message("--- Starting NOOB normalization without dye correction ---")
    load(opt$NOOB_dyeCorrFalse_MSet_Rdata)
    NOOB_dyeCorrFalse_normalisedBetas <- getBeta(NOOB_Mset[,controls])
    message("Extracted beta values for NOOB.")
    rm(NOOB_Mset)
    message("Removed NOOB_Mset from memory.")
    save(NOOB_dyeCorrFalse_normalisedBetas, file = file.path(dirname(opt$NOOB_dyeCorrFalse_MSet_Rdata), "NOOB_dyeCorrFalse_normalised_control_betas.RData"))
}


# NOOB Dye Correction Reference
# Assuming you've corrected the option name in your setup script to match (e.g., --NOOB_dyeCorrReference_MSet_Rdata)
if (!is.null(opt$NOOB_dyeCorrReference_MSet_Rdata)) {
    message("--- Starting NOOB normalization with dye correction reference ---")
    load(opt$NOOB_dyeCorrReference_MSet_Rdata)
    NOOB_dyeCorrReference_normalisedBetas <- getBeta(NOOB_Mset[,controls])
    message("Extracted beta values for NOOB.")
    rm(NOOB_Mset)
    message("Removed NOOB_Mset from memory.")
    save(NOOB_dyeCorrReference_normalisedBetas, file = file.path(dirname(opt$NOOB_dyeCorrReference_MSet_Rdata), "NOOB_dyeCorrReference_normalised_control_betas.RData"))
}


# NOOB Dye Correction Single
# Assuming you've corrected the option name in your setup script to match (e.g., --NOOB_dyeCorrSingle_MSet_Rdata)
if (!is.null(opt$NOOB_dyeCorrSingle_MSet_Rdata)) {
    message("--- Starting NOOB normalization with dye correction single ---")
    load(opt$NOOB_dyeCorrSingle_MSet_Rdata)
    NOOB_dyeCorrSingle_normalisedBetas <- getBeta(NOOB_Mset[,controls])
    message("Extracted beta values for NOOB.")
    rm(NOOB_Mset)
    message("Removed NOOB_Mset from memory.")
    save(NOOB_dyeCorrSingle_normalisedBetas, file = file.path(dirname(opt$NOOB_dyeCorrSingle_MSet_Rdata), "NOOB_dyeCorrSingle_normalised_control_betas.RData"))
}




# BMIQ is an intra-sample normalisation procedure, correcting the bias of type-2 probe values. BMIQ uses a 3-step procedure: (i) fitting of a 3-state beta mixture model, (ii) transformation of state-membership probabilities of type2 probes into quantiles of the type1 distribution, and (iii) a conformal transformation for the hemi-methylated probes. Exact details can be found in the reference below.
# Champ norm gives you opportunity to normalise using various method but on the back end it works with processXXX function. Other methods like BMIQ and PBC are called internally.
# Champ norm descrtiption -- Option to normalize data with a selection of normalization methods. There are four functions could be selected: "PBC","BMIQ","SWAN" and "FunctionalNormalize". SWAN method call for BOTH rgSet and mset input, FunctionNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. BMIQ method is the default function, which would also return normalised density plots in PDF format in results Dir. FunctionalNormalize is provided in minfi package, which ONLY support 450K data yet. Not that BMIQ function might fail if you sample's beta value distribution is not beta distribution, which occationally happen when too many CpGs are deleted while loading .idat files with champ.load() function. Also multi-cores parallel is conductable for BMIQ function, if your server or computer is good enought with more than one cores, you may assign more cores like 10 to accelerate the process. No matter what method you selected, they all will return the same result: Normalize beta matrix with effect of Type-I and Type-II probes corrected.

if (!is.null(opt$BMIQ_normalised_beta_matrix)) {
    message("--- Starting BMIQ normalization ---")
    load(opt$BMIQ_normalised_beta_matrix)
    message("Loaded BMIQ_normalised_beta_matrix.")
    BMIQ_normalisedBetas <- BMIQ_beta_matrix[,controls]
    message("Extracted beta values for BMIQ.")
    rm(BMIQ_beta_matrix)
    message("Removed BMIQ_beta_matrix from memory.")
    save(BMIQ_normalisedBetas, file = file.path(dirname(opt$BMIQ_normalised_beta_matrix), "BMIQ_normalised_control_betas.RData"))
}


# PBC
# peak-based correction
# Peak-based correction (PBC) is a method used to normalize methylation data, where the method rescales the methylation levels of Type II probes to match the distributions of Type I probes, aiming to make them more comparable. This is achieved by estimating the modes (peaks) of the methylation level distributions for both Type I and Type II probes and adjusting Type II probe values accordingly. 
if (!is.null(opt$PBC_normalised_beta_matrix)) {
    message("--- Starting PBC normalization ---")
    load(opt$PBC_normalised_beta_matrix)
    message("Loaded PBC_normalised_beta_matrix.")
    PBC_normalisedBetas <- PBC_beta_matrix[,controls]
    message("Extracted beta values for PBC.")
    rm(PBC_beta_matrix)
    message("Removed PBC_beta_matrix from memory.")
    save(PBC_normalisedBetas, file = file.path(dirname(opt$PBC_normalised_beta_matrix), "PBC_normalised_control_betas.RData"))
}


# ssNOOB + BMIQ
if (!is.null(opt$ssNoob_BMIQ_normalised_beta_matrix)) {
    message("--- Starting ssNOOB + BMIQ normalization ---")
    load(opt$ssNoob_BMIQ_normalised_beta_matrix)
    message("Loaded ssNoob_BMIQ_normalised_beta_matrix.")
    ssNoob_BMIQ_normalisedBetas <- BMIQ_beta_matrix[,controls]
    message("Extracted beta values for BMIQ.")
    rm(BMIQ_beta_matrix)
    message("Removed BMIQ_beta_matrix from memory.")
    save(ssNoob_BMIQ_normalisedBetas, file = file.path(dirname(opt$ssNoob_BMIQ_normalised_beta_matrix), "ssNoob_BMIQ_normalised_control_betas.RData"))
}


# pOOBAH NOOB Dye Correction Single
if (!is.null(opt$pOOBAH_NOOB_dyeCorrSingle_MSet_Rdata)) {
    message("--- Starting pOOBAH NOOB Dye Correction Single normalization ---")
    load(opt$pOOBAH_NOOB_dyeCorrSingle_MSet_Rdata)
    pOOBAH_NOOB_dyeCorrSingle_normalisedBetas <- BMIQ_beta_matrix[,controls]
    message("Extracted beta values for pOOBAH NOOB.")
    rm(BMIQ_beta_matrix)
    message("Removed BMIQ_beta_matrix from memory.")
    save(pOOBAH_NOOB_dyeCorrSingle_normalisedBetas, file = file.path(dirname(opt$pOOBAH_NOOB_dyeCorrSingle_MSet_Rdata), "pOOBAH_NOOB_dyeCorrSingle_normalised_control_betas.RData"))
}   

message("All beta matrices for control samples have been processed and saved.")