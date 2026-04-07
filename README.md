# Advancing Spacecraft Microbial Monitoring through Absolute Quantification and Shotgun Metagenomic Diversity Profiling

This repository contains the **R scripts and analysis-level data tables used to generate manuscript figure outputs** for:

> **"Advancing Spacecraft Microbial Monitoring through Absolute Quantification and Shotgun Metagenomic Diversity Profiling"**

---

## Overview

This is a **figure-reproducibility repository**. It contains the code and input data needed to regenerate manuscript figure outputs from analysis-level abundance tables. Raw sequencing files are not included.

Preprocessing scripts convert merged feature tables into standardized working matrices used by downstream plotting scripts. Because preprocessing has already been completed in this repository copy, **re-run the preprocessing scripts only if the source feature tables or metadata change**.

Two datasets are included, each with its own subdirectory, scripts, data tables, and output folder:

| Dataset | Figures |
|---------|---------|
| `PPR/` | Figures 2-5 |
| `KSC/` | Figure 6 |

## Data

This repository includes analysis-level data tables used by the figure scripts.

Files in each `data_tables/original/` folder are source analysis tables used as preprocessing inputs. They are not raw sequencing files.

---

## Repository structure

```text
Advancing_Spacecraft_Microbial_Monitoring/
|-- helperJ.R
|-- _globalStuff.R
|
|-- PPR/
|   |-- _globalStuff.R
|   |-- 1_feature_tables_to_working_matrix.R
|   |-- Figure_2_run first_percent_relative abundance.R
|   |-- Figure_2_run_second_Neg_Control.R
|   |-- Figure_3_run first_community_line_raw_rel.R
|   |-- Figure_3_run second_LOD.R
|   |-- Figure_4_run_first_percent_relative_abundance_separate3.R
|   |-- Figure_4_run_second__Release_Recovery.R
|   |-- Figure_5_JPL_SAF.R
|   |-- Figure_5_JPL_SAF - dqPCR.R
|   |-- R_package_citations.txt
|   |
|   |-- data_tables/
|   |   |-- meta.csv
|   |   |-- matrix_names.csv
|   |   |-- sample_reads.csv
|   |   |
|   |   |-- original/
|   |   |   |-- Merged_Phylum_raw_abund.csv
|   |   |   |-- Merged_Family_raw_abund.csv
|   |   |   |-- Merged_Genus_raw_abund.csv
|   |   |   \-- Merged_Species_raw_abund.csv
|   |   |
|   |   \-- working_matrix/
|   |       |-- working_matrix_Phylum.csv
|   |       |-- working_matrix_Family.csv
|   |       |-- working_matrix_Genus.csv
|   |       \-- working_matrix_Species.csv
|   |
|
|-- KSC/
|   |-- _globalStuff.R
|   |-- 1_feature_tables_to_working_matrix.R
|   |-- Figure_6_KSC_PHSF.R
|   |
|   |-- data_tables/
|   |   |-- meta.csv
|   |   |-- matrix_names.csv
|   |   |-- sample_reads.csv
|   |   |
|   |   |-- original/
|   |   |   |-- Merged_Phylum_raw_abund.csv
|   |   |   |-- Merged_Family_raw_abund.csv
|   |   |   |-- Merged_Genus_raw_abund.csv
|   |   |   \-- Merged_Species_raw_abund.csv
|   |   |
|   |   \-- working_matrix/
|   |       |-- working_matrix_Phylum.csv
|   |       |-- working_matrix_Family.csv
|   |       |-- working_matrix_Genus.csv
|   |       \-- working_matrix_Species.csv
|   |
\-- LICENSE
```

## File descriptions

### Root-level files

| File | Description |
|------|-------------|
| `helperJ.R` | Shared helper functions used across repository scripts |
| `_globalStuff.R` | Root-level setup script used by repository-wide workflows |

### PPR scripts

| File | Description |
|------|-------------|
| `PPR/_globalStuff.R` | PPR-specific setup: metadata loading, palettes, and shared plotting helpers |
| `PPR/1_feature_tables_to_working_matrix.R` | Converts merged feature tables into standardized working matrices and updates `matrix_names.csv` and `sample_reads.csv` |
| `PPR/Figure_2_run first_percent_relative abundance.R` | Prepares abundance objects used by the Figure 2 plotting step |
| `PPR/Figure_2_run_second_Neg_Control.R` | Writes Figure 2 negative-control outputs |
| `PPR/Figure_3_run first_community_line_raw_rel.R` | Writes Figure 3 community-profile outputs |
| `PPR/Figure_3_run second_LOD.R` | Writes Figure 3 limit-of-detection outputs |
| `PPR/Figure_4_run_first_percent_relative_abundance_separate3.R` | Prepares intermediate objects used by the Figure 4 plotting step |
| `PPR/Figure_4_run_second__Release_Recovery.R` | Writes Figure 4 release and recovery outputs |
| `PPR/Figure_5_JPL_SAF.R` | Writes the Figure 5A SAF plot outputs |
| `PPR/Figure_5_JPL_SAF - dqPCR.R` | Writes the Figure 5 dPCR/qPCR outputs |

### KSC scripts

| File | Description |
|------|-------------|
| `KSC/_globalStuff.R` | KSC-specific setup: metadata loading, palettes, and shared plotting helpers |
| `KSC/1_feature_tables_to_working_matrix.R` | Converts merged feature tables into standardized working matrices and updates `matrix_names.csv` and `sample_reads.csv` |
| `KSC/Figure_6_KSC_PHSF.R` | Writes the Figure 6 KSC-PHSF plot outputs |

---

## Required input files

Each dataset requires the following files in `data_tables/original/` and `data_tables/`:

| File | Description |
|------|-------------|
| `Merged_Phylum_raw_abund.csv` | Merged phylum-level raw abundance table |
| `Merged_Family_raw_abund.csv` | Merged family-level raw abundance table |
| `Merged_Genus_raw_abund.csv` | Merged genus-level raw abundance table |
| `Merged_Species_raw_abund.csv` | Merged species-level raw abundance table |
| `meta.csv` | Sample metadata |
| `matrix_names.csv` | Matrix lookup table updated by preprocessing |
| `sample_reads.csv` | Per-sample read counts updated by preprocessing |

> Files in `data_tables/original/` are analysis-level inputs, not raw sequencing files. Files in `data_tables/working_matrix/` are the standardized matrices used by figure scripts.

---

## Running the analysis

### Step 1 - Preprocessing

Run these scripts only if the input tables or metadata have changed.

```r
source("PPR/1_feature_tables_to_working_matrix.R")
source("KSC/1_feature_tables_to_working_matrix.R")
```

### Step 2 - Generate figure outputs

Run the relevant figure scripts in order. Scripts labeled `run first` and `run second` should be executed sequentially because later steps depend on objects created earlier in the workflow.

**PPR figures**

```r
# Figure 2
source("PPR/Figure_2_run first_percent_relative abundance.R")
source("PPR/Figure_2_run_second_Neg_Control.R")

# Figure 3
source("PPR/Figure_3_run first_community_line_raw_rel.R")
source("PPR/Figure_3_run second_LOD.R")

# Figure 4
source("PPR/Figure_4_run_first_percent_relative_abundance_separate3.R")
source("PPR/Figure_4_run_second__Release_Recovery.R")

# Figure 5
source("PPR/Figure_5_JPL_SAF.R")
source("PPR/Figure_5_JPL_SAF - dqPCR.R")
```

**KSC figures**

```r
# Figure 6
source("KSC/Figure_6_KSC_PHSF.R")
```

Scripts write generated outputs into dataset-level output folders during runtime as needed.

---

## R package citations

A list of package citations is included in `PPR/R_package_citations.txt`.

## Contact

Jeremy Kahsen - [Jeremy_Kahsen@Rush.edu](mailto:Jeremy_Kahsen@Rush.edu)
