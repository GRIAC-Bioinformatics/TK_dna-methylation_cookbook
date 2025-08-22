#!/usr/bin/env bash

# This script assume you will run each step as a slurm job due to
# resource intensity. You can remove the slurm specification to run
# this script like a normal .sh . 

#SBATCH --export=ALL
#SBATCH --job-name=dna.methylation.QC
#SBATCH --output=%x.log
#SBATCH --error=%x.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4

### 1. Initialize Conda ###
# Source Conda initialization script (adjust path if needed)
source ${PROJECT}/anaconda3/etc/profile.d/conda.sh

### 2. Activate your Conda environment ###
conda activate dna-methylation

### 3. Define variable ###
SCRIPTDIR=${PROJECT}/projects/VB250630_piama_dna_methylation/TK_dna-methylation_cookbook/scripts
source config.sh
# Check if the OUTDIR exists
if [ ! -d "${OUTDIR}" ]; then
  echo "Creating output directory: ${OUTDIR}"
  mkdir -p "${OUTDIR}"
else
  echo "Output directory already exists: ${OUTDIR}"
fi
# Check if the LOGDIR exists
if [ ! -d "${LOGDIR}" ]; then
  echo "Creating logging directory: ${LOGDIR}"
  mkdir -p "${LOGDIR}"
else
  echo "Logging directory already exists: ${LOGDIR}"
fi

### 4. Run all scripts ###
### STEP 0 : Check Input ###
Rscript ${SCRIPTDIR}/00_check_my_input.R --platform ${PLATFORM} \
                         --assembly ${ASSEMBLY} \
                         --idat_dir ${IDATDIR} \
                         --metadatasheet ${METADATASHEET} \
                         --outdir ${OUTDIR} \
                         --manifestkey ${MANIFEST_KEY} 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/00_check_my_input.$(date +%Y-%m-%d_%H-%M-%S).log"

# Check if the rdata exists
if [ ! -d "${OUTDIR}/rdata" ]; then
  echo "Creating logging directory: ${OUTDIR}/rdata"
  mkdir -p "${OUTDIR}/rdata"
else
  echo "Logging directory already exists: ${OUTDIR}/rdata"
fi
### STEP 1 : Create required objects ###
Rscript ${SCRIPTDIR}/01_create_rgset_object.R --platform ${PLATFORM} \
                                              --assembly ${ASSEMBLY} \
                                              --idat_dir ${IDATDIR} \
                                              --metadatasheet ${METADATASHEET} \
                                              --rg_channel_set_output ${OUTDIR}/rdata/RGChannelSet.RData \
                                              --methyl_channel_set_output ${OUTDIR}/rdata/MethylChannelSet.RData \
                                              --genomic_ratio_set_output ${OUTDIR}/rdata/GenomicRatioSet.RData \
                                              --rg_channel_set_extended_output ${OUTDIR}/rdata/RGChannelSetExtended.RData \
                                              --manifestkey ${MANIFEST_KEY} 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/01_create_rgset_object.$(date +%Y-%m-%d_%H-%M-%S).log"

# Check if the pdf exists
if [ ! -d "${OUTDIR}/pdf" ]; then
  echo "Creating logging directory: ${OUTDIR}/pdf"
  mkdir -p "${OUTDIR}/pdf"
else
  echo "Logging directory already exists: ${OUTDIR}/pdf"
fi

### STEP 2 : Meta data overview ###
Rscript ${SCRIPTDIR}/02_meta_data_overview.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                             --pca "Plate,Well,Gender,Array,Slide" \
                                             --stackbarplot "Plate:Gender,Slide:Gender" \
                                             --platecol "Plate" \
                                             --slidecol "Slide" \
                                             --output ${OUTDIR}/pdf/02_meta_data_overview.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/02_meta_data_overview.$(date +%Y-%m-%d_%H-%M-%S).log"

# Check if the flag exists
if [ ! -d "${OUTDIR}/flag" ]; then
  echo "Creating logging directory: ${OUTDIR}/flag"
  mkdir -p "${OUTDIR}/flag"
else
  echo "Logging directory already exists: ${OUTDIR}/flag"
fi

# Check if the flag/sample exists
if [ ! -d "${OUTDIR}/flag/sample" ]; then
  echo "Creating logging directory: ${OUTDIR}/flag/sample"
  mkdir -p "${OUTDIR}/flag/sample"
else
  echo "Logging directory already exists: ${OUTDIR}/flag/sample"
fi

### STEP 3 : Sample flags ###

# Sample Intensity #
Rscript ${SCRIPTDIR}/03_intensity_flagged_samples.R --mset ${OUTDIR}/rdata/MethylChannelSet.RData \
                                                    --cutoff 10.5 \
                                                    --flagged ${OUTDIR}/flag/sample/sample_intensity.csv \
                                                    --pdf_output ${OUTDIR}/pdf/03_intensity_flagged_samples.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_intensity_flagged_samples.$(date +%Y-%m-%d_%H-%M-%S).log"


# Sample Sex Mismatch #
Rscript ${SCRIPTDIR}/03_sample_sex_mismatch.R --grset ${OUTDIR}/rdata/GenomicRatioSet.RData \
                                              --flagged ${OUTDIR}/flag/sample/sex_mismatch.csv \
                                              --pdf_output ${OUTDIR}/pdf/03_sample_sex_mismatch.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_sex_mismatch.$(date +%Y-%m-%d_%H-%M-%S).log"


# Detection P value #
Rscript ${SCRIPTDIR}/03_sample_detection_p_value.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                   --flagged ${OUTDIR}/flag/sample/detectionP.csv \
                                                   --detectionP ${OUTDIR}/rdata/detectionP.Rdata \
                                                   --cutoff 0.01 \
                                                   --pdf_output ${OUTDIR}/pdf/03_sample_detection_p_value.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_detection_p_value.$(date +%Y-%m-%d_%H-%M-%S).log"


# Bisulfite Conversion #
Rscript ${SCRIPTDIR}/03_sample_dependent_controls_bisulfite_conversion.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                                         --platform ${PLATFORM} \
                                                                         --assembly ${ASSEMBLY} \
                                                                         --manifestkey ${MANIFEST_KEY} \
                                                                         --bs1_output ${OUTDIR}/flag/sample/bisulfite_conversion_I.csv \
                                                                         --bs2_output ${OUTDIR}/flag/sample/bisulfite_conversion_II.csv \
                                                                         --pdf_output ${OUTDIR}/pdf/03_sample_dependent_controls_bisulfite_conversion.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_bisulfite_conversion.$(date +%Y-%m-%d_%H-%M-%S).log"


# Negative Controls #
Rscript ${SCRIPTDIR}/03_sample_dependent_controls_negative.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                             --platform ${PLATFORM} \
                                                             --assembly ${ASSEMBLY} \
                                                             --manifestkey ${MANIFEST_KEY} \
                                                             --pdf_output ${OUTDIR}/pdf/03_sample_dependent_controls_negative.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_negative.$(date +%Y-%m-%d_%H-%M-%S).log"


# Non-Polymeric Controls #
Rscript ${SCRIPTDIR}/03_sample_dependent_controls_nonpolymeric.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                                 --platform ${PLATFORM} \
                                                                 --assembly ${ASSEMBLY} \
                                                                 --manifestkey ${MANIFEST_KEY} \
                                                                 --flagged ${OUTDIR}/flag/sample/nonpolymeric.csv \
                                                                 --pdf_output ${OUTDIR}/pdf/03_sample_dependent_controls_nonpolymeric.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_nonpolymeric.$(date +%Y-%m-%d_%H-%M-%S).log"


# Specificity Controls #
Rscript ${SCRIPTDIR}/03_sample_dependent_controls_specificity.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                                --platform ${PLATFORM} \
                                                                --assembly ${ASSEMBLY} \
                                                                --manifestkey ${MANIFEST_KEY} \
                                                                --sp1_output ${OUTDIR}/flag/sample/specificity_I.csv \
                                                                --sp2_output ${OUTDIR}/flag/sample/specificity_II.csv \
                                                                --pdf_output ${OUTDIR}/pdf/03_sample_dependent_controls_specificity.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_specificity.$(date +%Y-%m-%d_%H-%M-%S).log"


# Extension Controls #
Rscript  ${SCRIPTDIR}/03_sample_independent_controls_extension.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                                 --platform ${PLATFORM} \
                                                                 --assembly ${ASSEMBLY} \
                                                                 --manifestkey ${MANIFEST_KEY} \
                                                                 --flagged ${OUTDIR}/flag/sample/extension.csv \
                                                                 --pdf_output ${OUTDIR}/pdf/03_sample_independent_controls_extension.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_extension.$(date +%Y-%m-%d_%H-%M-%S).log"


# Hybridization Controls #
Rscript  ${SCRIPTDIR}/03_sample_independent_controls_hybridization.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                                     --platform ${PLATFORM} \
                                                                     --assembly ${ASSEMBLY} \
                                                                     --manifestkey ${MANIFEST_KEY} \
                                                                     --flagged ${OUTDIR}/flag/sample/hybridization.csv \
                                                                     --pdf_output ${OUTDIR}/pdf/03_sample_independent_controls_hybridization.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_hybridization.$(date +%Y-%m-%d_%H-%M-%S).log"

# Staining Controls #
Rscript ${SCRIPTDIR}/03_sample_independent_controls_staining.R --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                               --platform ${PLATFORM} \
                                                               --assembly ${ASSEMBLY} \
                                                               --manifestkey ${MANIFEST_KEY} \
                                                               --flagged ${OUTDIR}/flag/sample/staining.csv \
                                                               --pdf_output ${OUTDIR}/pdf/03_sample_independent_controls_staining.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_staining.$(date +%Y-%m-%d_%H-%M-%S).log"

# Target removal Controls #
Rscript ${SCRIPTDIR}/03_sample_independent_controls_target_removal.R  --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                                      --platform ${PLATFORM} \
                                                                      --assembly ${ASSEMBLY} \
                                                                      --manifestkey ${MANIFEST_KEY} \
                                                                      --flagged ${OUTDIR}/flag/sample/target_removal.csv \
                                                                      --pdf_output ${OUTDIR}/pdf/03_sample_independent_controls_target_removal.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_target_removal.$(date +%Y-%m-%d_%H-%M-%S).log"


### STEP 4 : Flagged sample overview and filtering ###

# Make a CSV file with all flagged samples to avoid provided indivisual CSV files #
for file in "${OUTDIR}/flag/sample/"*.csv; do
    # Get the base filename with the extension
    filename_with_ext=$(basename "$file")
    
    # Remove the .csv extension to get the base filename
    filename_without_ext="${filename_with_ext%.csv}"
    
    # Append the filename and full path to the CSV file
    echo "${filename_without_ext},${file}" >> "${OUTDIR}/flag/all_flagged_samples_list.csv"
done

# Flagged sample Overview #
Rscript ${SCRIPTDIR}/04_sample_qc_overview.R  --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                              --platform ${PLATFORM} \
                                              --assembly ${ASSEMBLY} \
                                              --manifestkey ${MANIFEST_KEY} \
                                              --pca_vars "Plate,Well,Gender,Slide,Array" \
                                              --pca_rdata ${OUTDIR}/rdata/pca.RData \
                                              --flaggedcombined ${OUTDIR}/flag/all_flagged_samples.csv \
                                              --flaggedlist ${OUTDIR}/flag/all_flagged_samples_list.csv \
                                              --pdf_output ${OUTDIR}/pdf/04_sample_qc_overview.pdf \
                                              --min_flag_overlap 2  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/04_sample_qc_overview.$(date +%Y-%m-%d_%H-%M-%S).log"

### CREATE excluded_samples.csv ###
# The creation of the excluded_samples.csv file is not done in the script
# but is assumed to be done by the user manually or by another script.
# This is to ensure that the user has control over which samples to exclude.
# As the reason to exclude samples can vary, it is not hardcoded in the script.

# Flagged sample filter #
Rscript ${SCRIPTDIR}/04_filter_flagged_samples.R  --rgset ${OUTDIR}/rdata/RGChannelSet.RData \
                                                  --mset ${OUTDIR}/rdata/MethylChannelSet.RData \
                                                  --grset ${OUTDIR}/rdata/GenomicRatioSet.RData \
                                                  --rgsetext ${OUTDIR}/rdata/RGChannelSetExtended.RData \
                                                  --flagged ${OUTDIR}/flag/excluded_samples.csv \
                                                  --base_suffix "_sample_filtered" \
                                                  --pdf_output ${OUTDIR}/pdf/04_filter_flagged_samples.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/04_filter_flagged_samples.$(date +%Y-%m-%d_%H-%M-%S).log"


# Check if the flag/probe exists
if [ ! -d "${OUTDIR}/flag/probe" ]; then
  echo "Creating logging directory: ${OUTDIR}/flag/probe"
  mkdir -p "${OUTDIR}/flag/probe"
else
  echo "Logging directory already exists: ${OUTDIR}/flag/probe"
fi


### STEP 5 : Probe flags ###
Rscript ${SCRIPTDIR}/05_detection_pvalue_flagged_probes.R --rgset ${OUTDIR}/rdata/RGChannelSet_sample_filtered.RData \
                                                          --flagged ${OUTDIR}/flag/probe/detectionP.csv \
                                                          --cutoff 0.01 \
                                                          --threshold 0.99\
                                                          --pdf_output ${OUTDIR}/pdf/05_detection_pvalue_flagged_probes.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_detection_pvalue_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"



Rscript ${SCRIPTDIR}/05_high_intensity_flagged_probes.R --mset ${OUTDIR}/rdata/MethylChannelSet_sample_filtered.RData \
                                                        --flagged ${OUTDIR}/flag/probe/high_intensity.csv \
                                                        --cutoff 10000 \
                                                        --pdf_output ${OUTDIR}/pdf/05_high_intensity_flagged_probes.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_high_intensity_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"


Rscript ${SCRIPTDIR}/05_low_beadcount_flagged_probes.R  --rgsetext ${OUTDIR}/rdata/RGChannelSetExtended_sample_filtered.RData \
                                                        --flagged ${OUTDIR}/flag/probe/low_beadcount.csv \
                                                        --cutoff 10000 \
                                                        --pdf_output ${OUTDIR}/pdf/05_low_beadcount_flagged_probes.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_low_beadcount_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"


Rscript ${SCRIPTDIR}/05_snp_containing_flagged_probes.R --grset ${OUTDIR}/rdata/GenomicRatioSet_sample_filtered.RData \
                                                        --flagged ${OUTDIR}/flag/probe/snp_containing_SBE_CpG.csv \
                                                        --cutoff 0.1 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_snp_containing_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"



### STEP 6 : Flagged probes overview and filtering ###

# Make a CSV file with all flagged samples to avoid provided indivisual CSV files #
for file in "${OUTDIR}/flag/probe/"*.csv; do
    # Get the base filename with the extension
    filename_with_ext=$(basename "$file")
    
    # Remove the .csv extension to get the base filename
    filename_without_ext="${filename_with_ext%.csv}"
    
    # Append the filename and full path to the CSV file
    echo "${filename_without_ext},${file}" >> "${OUTDIR}/flag/all_flagged_probes_list.csv"
done

# Flagged sample Overview #
Rscript ${SCRIPTDIR}/06_probe_qc_overview.R  --flaggedcombined ${OUTDIR}/flag/all_flagged_probes.csv \
                                              --flaggedlist ${OUTDIR}/flag/all_flagged_probes_list.csv \
                                              --pdf_output ${OUTDIR}/pdf/06_probe_qc_overview.pdf \
                                              --min_flag_overlap 2  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/06_probe_qc_overview.$(date +%Y-%m-%d_%H-%M-%S).log"

### CREATE excluded_probes.csv ###
# The creation of the excluded_probes.csv file is not done in the script
# but is assumed to be done by the user manually or by another script.
# This is to ensure that the user has control over which probes to exclude.
# As the reason to exclude probes can vary, it is not hardcoded in the script.

# Flagged sample filter #
Rscript ${SCRIPTDIR}/06_filter_flagged_probes.R --rgset ${OUTDIR}/rdata/RGChannelSet_sample_filtered.RData \
                                                  --mset ${OUTDIR}/rdata/MethylChannelSet_sample_filtered.RData \
                                                  --grset ${OUTDIR}/rdata/GenomicRatioSet_sample_filtered.RData \
                                                  --rgsetext ${OUTDIR}/rdata/RGChannelSetExtended_sample_filtered.RData \
                                                  --flagged ${OUTDIR}/flag/excluded_probes.csv \
                                                  --base_suffix "_probe_filtered" \
                                                  --pdf_output ${OUTDIR}/pdf/06_filter_flagged_probes.R.pdf 2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/06_filter_flagged_probes.R.$(date +%Y-%m-%d_%H-%M-%S).log"