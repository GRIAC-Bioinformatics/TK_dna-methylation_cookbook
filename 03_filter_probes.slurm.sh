#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --job-name=dna.methylation.final_filter
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

mkdir -p "${OUTDIR}/pdf" "${OUTDIR}/rdata" "${LOGDIR}"

if [ ! -f "${OUTDIR}/flag/excluded_probes.csv" ]; then
  echo "ERROR: ${OUTDIR}/flag/excluded_probes.csv does not exist."
  echo "Create it manually after reviewing Job 2 outputs."
  exit 1
fi

### STEP 6 : Filter flagged probes ###
Rscript "${SCRIPTDIR}/06_filter_flagged_probes.R" \
  --rgset "${OUTDIR}/rdata/RGChannelSet_sample_filtered.RData" \
  --mset "${OUTDIR}/rdata/MethylChannelSet_sample_filtered.RData" \
  --grset "${OUTDIR}/rdata/GenomicRatioSet_sample_filtered.RData" \
  --rgsetext "${OUTDIR}/rdata/RGChannelSetExtended_sample_filtered.RData" \
  --flagged "${OUTDIR}/flag/excluded_probes.csv" \
  --base_suffix "_probe_filtered" \
  --pdf_output "${OUTDIR}/pdf/06_filter_flagged_probes.pdf" \
  2>&1 | ts '[%Y-%m-%d %H:%M:%S]' >> "${LOGDIR}/06_filter_flagged_probes.$(date +%Y-%m-%d_%H-%M-%S).log"

echo "======================================================"
echo "JOB 3 COMPLETE"
echo "Final probe-filtered objects have been created."
echo "======================================================"