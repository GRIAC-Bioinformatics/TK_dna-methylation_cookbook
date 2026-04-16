## 🧬 DNA Methylation Array QC Cookbook

[![R](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![Illumina Arrays](https://img.shields.io/badge/platform-Illumina%20Methylation%20Arrays-orange.svg)](https://www.illumina.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A reproducible **Quality Control (QC) pipeline** for Illumina DNA methylation arrays (including **450K**, **EPIC**, and **EPICv2**).

This cookbook provides modular R scripts and a **human-in-the-loop SLURM workflow** to preprocess, filter, and QC your DNA methylation array data at both the **sample** and **probe** levels.

---

## 📂 Repository Structure

```text
├── config.sh
├── LICENSE
├── README.md
├── 01_run_until_sample_review.slurm.sh
├── 02_run_until_probe_review.slurm.sh
├── 03_filter_probes.slurm.sh
├── data
│   ├── EPICV2_probes_950K_CrossHybridization.csv
│   ├── idat
│   │   ├── 9374343009
│   │   │   ├── GSM8533627_9374343009_R04C02_Grn.idat
│   │   │   ├── GSM8533627_9374343009_R04C02_Red.idat
│   │   │   ├── GSM8533629_9374343009_R06C02_Grn.idat
│   │   │   └── GSM8533629_9374343009_R06C02_Red.idat
│   │   └── readme.sh
│   ├── manifest.annotation.key.csv
│   └── metadatasheet.csv
├── env
│   └── dna-methylation.yml
└── scripts
    ├── 00_check_my_input.R
    ├── 01_create_rgset_object.R
    ├── 02_meta_data_overview.R
    ├── 03_intensity_flagged_samples.R
    ├── 03_sample_dependent_controls_bisulfite_conversion.R
    ├── 03_sample_dependent_controls_negative.R
    ├── 03_sample_dependent_controls_nonpolymeric.R
    ├── 03_sample_dependent_controls_specificity.R
    ├── 03_sample_detection_p_value.R
    ├── 03_sample_independent_controls_extension.R
    ├── 03_sample_independent_controls_hybridization.R
    ├── 03_sample_independent_controls_staining.R
    ├── 03_sample_independent_controls_target_removal.R
    ├── 03_sample_sex_mismatch.R
    ├── 04_filter_flagged_samples.R
    ├── 04_sample_qc_overview.R
    ├── 05_detection_pvalue_flagged_probes.R
    ├── 05_high_intensity_flagged_probes.R
    ├── 05_low_beadcount_flagged_probes.R
    ├── 05_sex_chromosome_probes.R
    ├── 05_snp_containing_flagged_probes.R
    ├── 06_filter_flagged_probes.R
    ├── 06_probe_qc_overview.R
    └── helper_functions
        ├── convert_epic_samplesheet.R
        └── merge_excluded_probes.R
```

- **scripts/** → All R scripts for Illumina preprocessing and QC.
- **scripts/helper_functions** → Helper functions to create imput files.
- **Three SLURM entry points** are provided to match the two required manual inspection steps in the workflow.
- **data/** → Example dataset with:
  - Example 450K IDAT files
  - `metadatasheet.csv` (sample metadata)
  - `manifest.annotation.key.csv` (maps to Illumina manifests and annotation files for 450K/EPIC/EPICv2)

---

## 🚀 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/GRIAC-Bioinformatics/TK_dna-methylation_cookbook.git
cd TK_dna-methylation_cookbook
```

### 2. Prepare environment

We recommend using the provided YAML environment file to ensure reproducibility.

```bash
conda env create -f env/dna-methylation.yml
conda activate dna-methylation
```

On HPC:

```bash
ml Anaconda3
export CONDA_PKGS_DIRS=/groups/${choose_one_of_your_groups}/tmp02/some/example/folder/conda/packages/
conda info
conda env create -f env/dna-methylation.yml --prefix /groups/${choose_one_of_your_groups}/tmp02/some/example/folder/conda/environment/
```

### 3. Configure

Edit `config.sh` to specify:

- Platform (`450K`, `EPIC`, `EPICv2`)
- Assembly (`hg19` or `hg38`) — **EPICv2 is always hg38**
- Paths to your IDATs, metadata, manifest key (`data/manifest.annotation.key.csv`), output directory, and logging directory

### 4. Prepare metadata

Create a metadata file similar to [data/metadatasheet.csv](data/metadatasheet.csv). Make sure that the column headers match.

For EPIC metadata, a conversion script is available [here](scripts/helper_functions/convert_epic_samplesheet.R).

If your file looks similar to:

```text
[Header]
Investigator Name,,XXXXXXXX
Project Name,,XXXXXXXX
Experiment Name,,
Date,,XX_XX_XXXX
[Data]
Sample_Well,Sample_Name,Sample_Group,Sample_Plate,Sentrix_ID,Sentrix_Position,Project,Gender,Pool_ID
A1,L101817420,,OV2300861CDNA03,209548820026,R01C01,S2025-71 EPIC,F,
A1,FS57750722,,OV2300861CDNA04,209536060010,R01C01,S2025-71 EPIC,F,
```

Usage example:

```bash
chmod +x convert_epic_samplesheet.R
./convert_epic_samplesheet.R \
  --input /path/to/metadata.csv \
  --output /path/to/output/metadatasheet.csv \
  --basename-prefix /path/to/idats \
  --check-files \
  --stop-on-duplicate
```

### 5. Run pipeline

The workflow is intentionally split into **3 SLURM jobs** because there are **2 manual review checkpoints**:

1. **Job 1** runs input checks, object creation, metadata overview, sample-level QC, and the sample QC overview.
2. **Manual checkpoint 1**: review sample QC outputs and create `excluded_samples.csv`.
3. **Job 2** filters samples, runs probe-level QC, and creates the probe QC overview.
4. **Manual checkpoint 2**: review probe QC outputs and create `excluded_probes.csv`.
5. **Job 3** filters probes and creates the final probe-filtered objects.

#### HPC / SLURM usage

```bash
sbatch 01_run_until_sample_review.slurm.sh
```

After manually creating `flag/excluded_samples.csv`:

```bash
sbatch 02_run_until_probe_review.slurm.sh
```

After manually creating `flag/excluded_probes.csv`:

```bash
sbatch 03_filter_probes.slurm.sh
```

#### Running locally (without SLURM)

Remove or adapt the `#SBATCH` lines and execute each script in the same order:

```bash
bash 01_run_until_sample_review.slurm.sh
bash 02_run_until_probe_review.slurm.sh
bash 03_filter_probes.slurm.sh
```

---

## 🧩 Workflow Logic

The scripts are modular and ordered numerically. You can run them individually or in blocks, but the recommended workflow is the **3-job split** because filtering requires manual curation at two stages.

### 🧩 Workflow Steps

| Step | Purpose | Notes |
|------|---------|-------|
| **00_** | Input checks | Validate raw IDATs and metadata consistency |
| **01_** | Create RGSet | Build raw RGSet object (starting point for downstream QC) |
| **02_** | Metadata overview | Summarize and visualize sample sheet information |
| **03_** | Sample-level QC | Evaluate control probes (bisulfite conversion, staining, sex check, etc.) ✅ Can be run in parallel |
| **04_** | Filter samples | Exclude low-quality samples ⚠️ Requires **manual review** |
| **05_** | Probe-level QC | Identify problematic probes (SNPs, low beadcount, detection p-value, intensity) ✅ Can be run in parallel |
| **06_** | Filter probes | Exclude low-quality probes ⚠️ Requires **manual review** |

### Recommended execution model

#### **Job 1 — run until sample review**
Runs:
- `00_check_my_input.R`
- `01_create_rgset_object.R`
- `02_meta_data_overview.R`
- all `03_*` sample QC scripts
- `04_sample_qc_overview.R`

Then manually inspect:
- `pdf/04_sample_qc_overview.pdf`
- `flag/all_flagged_samples.csv`
- `flag/all_flagged_samples_list.csv`

Create:
- `flag/excluded_samples.csv`

#### **Job 2 — run until probe review**
Runs:
- `04_filter_flagged_samples.R`
- all `05_*` probe QC scripts
- `06_probe_qc_overview.R`

Then manually inspect:
- `pdf/06_probe_qc_overview.pdf`
- `flag/all_flagged_probes.csv`
- `flag/all_flagged_probes_list.csv`

Create:
- `flag/excluded_probes.csv`

#### **Job 3 — final probe filtering**
Runs:
- `06_filter_flagged_probes.R`

Output:
- final sample- and probe-filtered `.RData` objects for downstream analysis

---

## ⚠️ Important

- For `06_filter_flagged_probes.R` and **EPICv2 / 950K**, also filter out potential cross-hybridization probes, see: [VB_Illumina_DNAmethylation_450k_850k_probe_filtering_cookbook](https://github.com/GRIAC-Bioinformatics/VB_Illumina_DNAmethylation_450k_850k_probe_filtering_cookbook?tab=readme-ov-file#epicv2950k-array)
- After removal of bad samples and probes, the data should be normalized between samples. See: [VB_nf_dna_methylation_sesame_pipeline](https://github.com/GRIAC-Bioinformatics/VB_nf_dna_methylation_sesame_pipeline)

---

## 📊 Outputs

The pipeline generates:

- QC reports and plots for each step
- Logs of SLURM/bash runs (`logs/`)
- Intermediate `.RData` objects after raw import and sample filtering
- Final sample- and probe-filtered `.RData` objects for downstream normalization and analysis

---

## 📚 References

- Agata Smialowska DNA Methylation: Array Workflow: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html#normalization
- Jovana Maksimovic pipeline: https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#quality-control
- FOXO pipeline: https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/docs/introduction/introduction.html#introduction
- Normalization benchmark paper: https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-023-01459-z
- Illumina - Leveraging Hidden Information to Correct for Background Fluorescence with Kim Siegmund: https://www.youtube.com/watch?v=_IWhwXnAAls
- Illumina manifest (product) files: https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html

---

## 🙋 Contact

Maintainer: Vartika Bisht  
📧 v.bisht@umcg.nl

✨ Happy QC’ing! Your methylation data deserves it ✨
