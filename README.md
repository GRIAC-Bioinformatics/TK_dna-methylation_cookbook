# Methylation Array QC Pipeline

## Overview

This repository contains scripts for performing quality control (QC) steps on DNA methylation array data (e.g., Illumina 450K or EPIC array). The pipeline includes key steps such as data preprocessing, sample level QC and probe-level QC.

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

### Running the Scripts

Depending on the size of the dataset you might need to run the scripts on HPC usong batch jobs. For more info visit this [link] https://docs.gcc.rug.nl/nibbler/analysis/

## Contact

For questions or issues, please contact [Tatiana Karp](t.karp@rug.nl), [Martin Banchero] (m.banchero@umcg.nl), [Maaike de Vries] (m.de.vries04@umcg.nl).
