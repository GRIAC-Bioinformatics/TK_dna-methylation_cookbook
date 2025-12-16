# 🧬 DNA Methylation Array QC Cookbook

[![R](https://img.shields.io/badge/R-%23276DC3.svg?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![Illumina Arrays](https://img.shields.io/badge/Illumina-450K%20|%20EPIC%20|%20EPICv2-orange)](https://www.illumina.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A reproducible **Quality Control (QC) pipeline** for Illumina DNA methylation arrays  
(including **450K**, **EPIC**, and **EPICv2**).  

This cookbook provides modular R scripts and a bash/SLURM workflow to help you preprocess, filter, and QC your DNA methylation array data at both the **sample** and **probe** levels.

---

## 📂 Repository Structure
```
├── config.sh
├── LICENSE
├── main.slurm.sh
├── README.md
├── data
│   ├── EPICV2_probes_950K_CrossHybridization.csv
│   ├── idat
│   │   ├── 9374343009
│   │   │   ├── GSM8533627_9374343009_R04C02_Grn.idat
│   │   │   ├── GSM8533627_9374343009_R04C02_Red.idat
│   │   │   ├── GSM8533629_9374343009_R06C02_Grn.idat
│   │   │   └── GSM8533629_9374343009_R06C02_Red.idat
│   │   └── readme.sh
│   ├── manifest.annotation.key.csv
│   └── metadatasheet.csv
├── env
│   └── dna-methylation.yml
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
    └── 06_probe_qc_overview.R

```

<details>
<summary>📝 About the directories</summary>

- **`scripts/`** → All R scripts for Illumina preprocessing and QC.  
  Each script is **modular**, with comments describing input, output, and QC checks.

- **`data/`** → Example dataset with:
  - Example 450K IDAT files
  - `metadatasheet.csv` (sample metadata)
  - `manifest.annotation.key.csv` (maps to Illumina manifests and annotation files for 450K/EPIC/EPICv2)

- **`main.slurm.sh`** → Entry point for running the full pipeline.  
  Includes SLURM job specs (remove/adapt if not using HPC).
</details>

---

## 🚀 Getting Started

### 1. Clone the repository
```
git clone https://github.com/GRIAC-Bioinformatics/TK_dna-methylation_cookbook.git
cd TK_dna-methylation_cookbook
```

### 2. Prepare environment

We recommend using the provided YAML environment file to ensure reproducibility.

```
conda env create -f dna-methylation.yml
conda activate dna-methylation
```

### 3. Configure

Edit config.sh to specify:

- Platform (450K, EPIC, EPICv2)
- Assembly (hg19 or hg38)
- Paths to your IDATs, metadata, manifest key ( data/manifest.annotation.key.csv ), output and logging directory

### 4. Run pipeline

If running on an HPC with SLURM:
```
sbatch main.slurm.sh
```

If running locally (without SLURM), remove the slurm specifications and simply execute:
```
bash main.slurm.sh
```
## 🧩 Workflow Logic

The scripts are modular and ordered numerically. You can run them individually or in blocks:

### 🧩 Workflow Steps

| Step  | Purpose            | Notes                                                                 |
|-------|--------------------|----------------------------------------------------------------------|
| **00_** | Input checks       | Validate raw IDATs and metadata consistency                         |
| **01_** | Create RGSet       | Build raw RGSet object (starting point for downstream QC)           |
| **02_** | Metadata overview  | Summarize and visualize sample sheet information                    |
| **03_** | Sample-level QC    | Evaluate control probes (bisulfite conversion, staining, sex check, etc.) <br> ✅ Can be run in parallel |
| **04_** | Filter samples     | Exclude low-quality samples <br> ⚠️ Requires **manual review**       |
| **05_** | Probe-level QC     | Identify problematic probes (SNPs, low beadcount, detection p-value, intensity) <br> ✅ Can be run in parallel |
| **06_** | Filter probes      | Exclude low-quality probes <br> ⚠️ Requires **manual review**        |


## ⚠️ Important:
04_filter_flagged_samples.R and 06_filter_flagged_probes.R require manual curation of flagged samples/probes before filtering. This safeguard prevents accidental data loss.

## 📊 Outputs

- QC reports and plots for each step
- Logs of SLURM/bash runs (logs/)
- Filtered RGSet objects for downstream analysis

## 📚 References
 - Agata Smialowska DNA Methylation: Array Workflow: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html#normalization
 - Jovana Maksimovic pipeline: https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#quality-control
 - FOXO pipeline: https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/docs/introduction/introduction.html#introduction
 - Noramalization benchmark paper: https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-023-01459-z
 - Illumina - Leveraging Hidden Information to Correct for Background Fluorescence with Kim Siegmund: https://www.youtube.com/watch?v=_IWhwXnAAls
 - Illumina manifest (product) files: https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html

## 🙋 Contact

Maintainer: Vartika Bisht
📧 v.bisht@umcg.nl

---

✨ Happy QC’ing! Your methylation data deserves it ✨
