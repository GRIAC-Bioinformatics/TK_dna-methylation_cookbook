#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --job-name=dna.methylation.probe_qc
#SBATCH --output=%x.log
#SBATCH --error=%x.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4

set -euo pipefail

ml Anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /groups/umcg-lifelines/tmp02/projects/ov21_0402/conda/dna_methylation_env

source config.sh

mkdir -p "${OUTDIR}/flag/probe" "${OUTDIR}/pdf" "${OUTDIR}/rdata" "${LOGDIR}"

if [ ! -f "${OUTDIR}/flag/excluded_samples.csv" ]; then
  echo "ERROR: ${OUTDIR}/flag/excluded_samples.csv does not exist."
  echo "Create it manually after reviewing Job 1 outputs."
  exit 1
fi

### STEP 4 : Filter flagged samples ###
Rscript "${SCRIPTDIR}/04_filter_flagged_samples.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet.RData" \
  --mset "${OUTDIR}/rdata/MethylChannelSet.RData" \
  --grset "${OUTDIR}/rdata/GenomicRatioSet.RData" \
  --rgsetext "${OUTDIR}/rdata/RGChannelSetExtended.RData" \
  --flagged "${OUTDIR}/flag/excluded_samples.csv" \
  --base_suffix "_sample_filtered" \
  --pdf_output "${OUTDIR}/pdf/04_filter_flagged_samples.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/04_filter_flagged_samples.$(date +%Y-%m-%d_%H-%M-%S).log"

### STEP 5 : Probe flags ###
Rscript "${SCRIPTDIR}/05_detection_pvalue_flagged_probes.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet_sample_filtered.RData" \
  --flagged "${OUTDIR}/flag/probe/detectionP.csv" \
  --cutoff 0.01 \
  --threshold 0.99 \
  --pdf_output "${OUTDIR}/pdf/05_detection_pvalue_flagged_probes.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_detection_pvalue_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/05_high_intensity_flagged_probes.R" \
  --mset "${OUTDIR}/rdata/MethylChannelSet_sample_filtered.RData" \
  --flagged "${OUTDIR}/flag/probe/high_intensity.csv" \
  --cutoff 10000 \
  --pdf_output "${OUTDIR}/pdf/05_high_intensity_flagged_probes.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_high_intensity_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/05_low_beadcount_flagged_probes.R" \
  --rgsetext "${OUTDIR}/rdata/RGChannelSetExtended_sample_filtered.RData" \
  --flagged "${OUTDIR}/flag/probe/low_beadcount.csv" \
  --cutoff 10000 \
  --pdf_output "${OUTDIR}/pdf/05_low_beadcount_flagged_probes.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_low_beadcount_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"

Rscript "${SCRIPTDIR}/05_snp_containing_flagged_probes.R" \
  --grset "${OUTDIR}/rdata/GenomicRatioSet_sample_filtered.RData" \
  --flagged "${OUTDIR}/flag/probe/snp_containing_SBE_CpG.csv" \
  --cutoff 0.1 \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/05_snp_containing_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"

### STEP 6 : Probe QC overview only ###
: > "${OUTDIR}/flag/all_flagged_probes_list.csv"
for file in "${OUTDIR}/flag/probe/"*.csv; do
    filename_with_ext=$(basename "$file")
    filename_without_ext="${filename_with_ext%.csv}"
    echo "${filename_without_ext},${file}" >> "${OUTDIR}/flag/all_flagged_probes_list.csv"
done

Rscript "${SCRIPTDIR}/06_probe_qc_overview.R" \
  --flaggedcombined "${OUTDIR}/flag/all_flagged_probes.csv" \
  --flaggedlist "${OUTDIR}/flag/all_flagged_probes_list.csv" \
  --pdf_output "${OUTDIR}/pdf/06_probe_qc_overview.pdf" \
  --min_flag_overlap 2 \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/06_probe_qc_overview.$(date +%Y-%m-%d_%H-%M-%S).log"

echo "======================================================"
echo "JOB 2 COMPLETE"
echo "Now manually inspect the probe QC results and create:"
echo "  ${OUTDIR}/flag/excluded_probes.csv"
echo "Then submit Job 3."
echo "======================================================"