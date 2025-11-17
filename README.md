# Advancing Spacecraft Microbial Monitoring through Absolute Quantification and Shotgun Metagenomic Diversity Profiling


**Note:** This repository is undergoing active maintenance as of 2025-11-17.  
A full, stable version will be available before December.


## Overview
This repository contains the R scripts, processed data, and reproducible workflows associated with the manuscript:

**"Advancing Spacecraft Microbial Monitoring through Absolute Quantification and Shotgun Metagenomic Diversity Profiling"**  
**Author:** Jeremy Kahsen (Jeremy_Kahsen@Rush.edu)  
**Manuscript status:** *[Under review / In preparation / Published]*  
**Affiliations:** *[Insert institution(s)]*

This project integrates:

- Digital PCR (dPCR) absolute quantification  
- Shotgun metagenomic diversity profiling  
- Amplicon-based microbiome analyses  
- Low-biomass contamination-aware workflows  
- Multi-environment spacecraft-associated samples (PPR, ILMAH, SAF, etc.)

The repository provides the complete computational workflow used to generate manuscript figures and tables.

---

## Repository Structure

```
/                                  # Root project directory
├── helper/                        # Shared helper functions
│   └── helperJ.R
│
├── PPR/                           # Dataset 1 (Planetary Protection Research)
│   ├── data_tables/               # Raw data for PPR
│   ├── output_data/               # Generated intermediate data
│   ├── output_plot/               # Publication-ready figures & tables
│   ├── scripts/                   # One R script per figure/table
│   └── _globalStuff.R             # PPR-specific global settings
│
├── JPL/                           # Dataset 2 (JPL dataset)
│   ├── data_tables/               # Raw data for JPL
│   ├── output_data/               # Generated intermediate data
│   ├── output_plot/               # Publication-ready figures & tables
│   ├── scripts/                   # One R script per figure/table
│   └── _globalStuff.R             # JPL-specific global settings
│
├── README.md                      # This file
└── Advancing_Spacecraft_Microbial_Monitoring.Rproj
```

---

## Reproducibility

### Requirements
- **R version:** ≥ 4.3  
- **Packages:**  
  `tidyverse`, `ggplot2`, `patchwork`, `vegan`, `compositions`,  
  `ggtext`, `ggh4x`, `effectsize`, plus those automatically loaded via helper scripts.

---

### Load shared helper functions

```r
source("helper/helperJ.R")
```

### Load dataset-specific global settings

**For PPR:**
```r
source("PPR/_globalStuff.R")
```

**For JPL:**
```r
source("JPL/_globalStuff.R")
```

### Run figure/table scripts

**Example:**
```r
source("PPR/scripts/Figure_5_Microbial_Load.R")
```

**Or for JPL:**
```r
source("JPL/scripts/Figure_3_Shotgun_Diversity.R")
```

### Outputs are saved automatically to:

```
PPR/output_data/
PPR/output_plot/

JPL/output_data/
JPL/output_plot/
```

---

## Placeholder Sections (Complete during manuscript finalization)

### Study Summary
*[Insert summary of sampling design, environments, and goals]*

### Methods Summary
*[Insert description of dPCR processing, shotgun workflows, QC steps, and statistical methods]*

### Citation (update once published)
*Kahsen J. et al., Year, Journal, DOI*

### Acknowledgments
*[Insert collaborators, funding, institutional support]*

### License
*[Insert license: MIT, CC-BY, or institution-specific requirements]*

---

## Contact
For questions or collaboration:

📧 **Jeremy_Kahsen@Rush.edu**

---

## Notes
- The repository is actively maintained.  
- Each figure/table corresponds to a dedicated script located in the appropriate dataset folder.  
- Shared helper functions are in `helper/`.  
- Global settings are dataset-specific and live inside `PPR/` and `JPL/`.
