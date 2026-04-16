#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --job-name=dna.methylation.sample_qc
#SBATCH --output=%x.log
#SBATCH --error=%x.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4

set -euo pipefail

### 1. Initialize Conda ###
ml Anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"

### 2. Activate environment ###
conda activate /groups/umcg-lifelines/tmp02/projects/ov21_0402/conda/dna_methylation_env

### 3. Source config ###
source config.sh

### 4. Create required directories ###
mkdir -p "${OUTDIR}" "${LOGDIR}" "${OUTDIR}/rdata" "${OUTDIR}/pdf" "${OUTDIR}/flag/sample"

### STEP 0 : Check Input ###
Rscript "${SCRIPTDIR}/00_check_my_input.R" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --idat_dir "${IDATDIR}" \
  --metadatasheet "${METADATASHEET}" \
  --outdir "${OUTDIR}" \
  --manifestkey "${MANIFEST_KEY}" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/00_check_my_input.$(date +%Y-%m-%d_%H-%M-%S).log"

### STEP 1 : Create required objects ###
Rscript "${SCRIPTDIR}/01_create_rgset_object.R" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --idat_dir "${IDATDIR}" \
  --metadatasheet "${METADATASHEET}" \
  --rg_channel_set_output "${OUTDIR}/rdata/RGChannelSet.RData" \
  --methyl_channel_set_output "${OUTDIR}/rdata/MethylChannelSet.RData" \
  --genomic_ratio_set_output "${OUTDIR}/rdata/GenomicRatioSet.RData" \
  --rg_channel_set_extended_output "${OUTDIR}/rdata/RGChannelSetExtended.RData" \
  --manifestkey "${MANIFEST_KEY}" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/01_create_rgset_object.$(date +%Y-%m-%d_%H-%M-%S).log"

### STEP 2 : Meta data overview ###
Rscript "${SCRIPTDIR}/02_meta_data_overview.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --pca "Plate,Well,Gender,Array,Slide" \
  --stackbarplot "Plate:Gender,Slide:Gender" \
  --platecol "Plate" \
  --slidecol "Slide" \
  --output "${OUTDIR}/pdf/02_meta_data_overview.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/02_meta_data_overview.$(date +%Y-%m-%d_%H-%M-%S).log"

### STEP 3 : Sample flags ###

Rscript "${SCRIPTDIR}/03_intensity_flagged_samples.R" \
  --mset "${OUTDIR}/rdata/MethylChannelSet.RData" \
  --cutoff 10.5 \
  --flagged "${OUTDIR}/flag/sample/sample_intensity.csv" \
  --pdf_output "${OUTDIR}/pdf/03_intensity_flagged_samples.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_intensity_flagged_samples.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_sex_mismatch.R" \
  --grset "${OUTDIR}/rdata/GenomicRatioSet.RData" \
  --flagged "${OUTDIR}/flag/sample/sex_mismatch.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_sex_mismatch.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_sex_mismatch.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_detection_p_value.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --flagged "${OUTDIR}/flag/sample/detectionP.csv" \
  --detectionP "${OUTDIR}/rdata/detectionP.Rdata" \
  --cutoff 0.01 \
  --pdf_output "${OUTDIR}/pdf/03_sample_detection_p_value.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_detection_p_value.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_dependent_controls_bisulfite_conversion.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --bs1_output "${OUTDIR}/flag/sample/bisulfite_conversion_I.csv" \
  --bs2_output "${OUTDIR}/flag/sample/bisulfite_conversion_II.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_dependent_controls_bisulfite_conversion.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_bisulfite_conversion.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_dependent_controls_negative.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --pdf_output "${OUTDIR}/pdf/03_sample_dependent_controls_negative.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_negative.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_dependent_controls_nonpolymeric.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --flagged "${OUTDIR}/flag/sample/nonpolymeric.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_dependent_controls_nonpolymeric.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_nonpolymeric.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_dependent_controls_specificity.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --sp1_output "${OUTDIR}/flag/sample/specificity_I.csv" \
  --sp2_output "${OUTDIR}/flag/sample/specificity_II.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_dependent_controls_specificity.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_dependent_controls_specificity.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_independent_controls_extension.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --flagged "${OUTDIR}/flag/sample/extension.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_independent_controls_extension.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_extension.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_independent_controls_hybridization.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --flagged "${OUTDIR}/flag/sample/hybridization.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_independent_controls_hybridization.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_hybridization.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_independent_controls_staining.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --flagged "${OUTDIR}/flag/sample/staining.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_independent_controls_staining.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_staining.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/03_sample_independent_controls_target_removal.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --flagged "${OUTDIR}/flag/sample/target_removal.csv" \
  --pdf_output "${OUTDIR}/pdf/03_sample_independent_controls_target_removal.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/03_sample_independent_controls_target_removal.$(date +%Y-%m-%d_%H-%M-%S).log"

### STEP 4 : Sample QC overview only ###
: > "${OUTDIR}/flag/all_flagged_samples_list.csv"
for file in "${OUTDIR}/flag/sample/"*.csv; do
    filename_with_ext=$(basename "$file")
    filename_without_ext="${filename_with_ext%.csv}"
    echo "${filename_without_ext},${file}" >> "${OUTDIR}/flag/all_flagged_samples_list.csv"
done

Rscript "${SCRIPTDIR}/04_sample_qc_overview.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --platform "${PLATFORM}" \
  --assembly "${ASSEMBLY}" \
  --manifestkey "${MANIFEST_KEY}" \
  --pca_vars "Plate,Well,Gender,Slide,Array" \
  --pca_rdata "${OUTDIR}/rdata/pca.RData" \
  --flaggedcombined "${OUTDIR}/flag/all_flagged_samples.csv" \
  --flaggedlist "${OUTDIR}/flag/all_flagged_samples_list.csv" \
  --pdf_output "${OUTDIR}/pdf/04_sample_qc_overview.pdf" \
  --min_flag_overlap 2 \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/04_sample_qc_overview.$(date +%Y-%m-%d_%H-%M-%S).log"

echo "======================================================"
echo "JOB 1 COMPLETE"
echo "Now manually inspect the sample QC results and create:"
echo "  ${OUTDIR}/flag/excluded_samples.csv"
echo "Then submit Job 2."
echo "======================================================"