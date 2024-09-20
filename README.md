# Methylation Array QC Pipeline

## Overview

This repository contains scripts for performing quality control (QC) steps on DNA methylation array data (e.g., Illumina 450K or EPIC or EPIC V2 array). The pipeline includes key steps such as data preprocessing, sample level QC and probe-level QC.

## Repository Structure

- **File 1: `01_data_preparation.R`**
  - **Purpose:** Preprocesses and filters the raw methylation data, preparing it for quality control.
  - **Key Steps:**
    - Load raw IDAT files
    - Create minfi objects: RGChannelSet, MethylSet, RatioSet
    - Upload annotation 

- **File 2: `02_sample_QC.R`**
  - **Purpose:** Conducts initial sample-level quality control.
  - **Key Steps:**
    - Check methylated/unmethylated intensities
    - Check detection P-values
    - Check sex concordance
    - Check technical variation
    - Check bisulfite conversion step
    - Check beta-values disitribution 

- **File 3: `03_probe_QC.R`**
  - **Purpose:** Conducts probe-level quality control.
  - **Key Steps:**
    - Check control probes
    - Check detection p-values
    - Check cross-reactive probes
    - Check SNP-containing probes
    - Check XY-related probes
    - Check high-intensity probes
    - Check probes with low bead counts

## Usage

### Prerequisites

- R (version 4.0 or later recommended)
- Required R packages: `minfi`, `tidyverse`, `ENmix`, `wateRmelon`
- Depending on the array version you need to use different manifest and annotation files.
- Depending on the array version you might need to use different list of cross-reactive probes

### Running the Scripts

Depending on the size of the dataset you might need to run the scripts on HPC using batch jobs. For more info about how to do this on nibbler visit this link: https://docs.gcc.rug.nl/nibbler/analysis/

## References: 

 - Agata Smialowska DNA Methylation: Array Workflow: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html#normalization
 - Jovana Maksimovic pipeline: https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#quality-control
 - FOXO pipeline: https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/docs/introduction/introduction.html#introduction
 - Noramalization benchmark paper: https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-023-01459-z
 - Illumina - Leveraging Hidden Information to Correct for Background Fluorescence with Kim Siegmund: https://www.youtube.com/watch?v=_IWhwXnAAls


## Contact

For questions or issues, please contact [Tatiana Karp](t.karp@rug.nl), [Martin Banchero](m.banchero@umcg.nl), [Maaike de Vries](m.de.vries04@umcg.nl).


