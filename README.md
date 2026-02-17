# Advancing Spacecraft Microbial Monitoring through Absolute Quantification and Shotgun Metagenomic Diversity Profiling

**Note:** This repository is undergoing active maintenance.

## Overview
This repository contains R scripts, processed data, and reproducible workflows used to generate manuscript figures and tables for spacecraft microbial monitoring analyses.

---

## Active Repository Structure

```text
/
├── PPR/                       # Primary dataset workflows and data
│   ├── _globalStuff.R
│   ├── 2_Feature_Tables_to_Rarefied_merged.R
│   └── data_tables/
│       ├── original/
│       └── raw_matrix/
│
├── JPL/                       # Secondary dataset workflows (renamed from Clipper)
│   ├── _globalStuff.R
│   ├── 2_Feature_Tables_to_Rarefied.R
│   └── data_tables/
│       ├── original/
│       └── raw_matrix/
│
├── helperJ.R                  # Shared helper functions
├── _globalStuff.R             # Root/global bootstrap script
├── actual figures/            # Canonical figure outputs to recreate
└── unneeded/                  # Archived legacy/duplicate/temp files
```

---

## Current Curation Rules

- Canonical manuscript figures are in `actual figures/`.
- Active feature-table preprocessing is script-driven from:
  - `PPR/2_Feature_Tables_to_Rarefied_merged.R`
  - `JPL/2_Feature_Tables_to_Rarefied.R`
- Each feature-table script sources its local `_globalStuff.R` (which sets working directory and loads `helperJ.R`).
- `look_at_counts`, rarefied output, and filtered output steps are intentionally removed from the active preprocessing path.
- The active preprocessing path now does only:
  1. create/update `data_tables/matrix_names.csv`
  2. transform `data_tables/original/Merged_*_raw_abund.csv` into `data_tables/raw_matrix/raw_matrix_*.csv`
- In `matrix_names.csv`, `file_name` and `file_path` now target the raw-matrix output path.
- Legacy JPL content and other non-essential files were moved to `unneeded/` for temporary retention.

---

## Required Input Files (per dataset)

Only these canonical merged taxonomic files are expected in each active dataset's `data_tables/original/` folder:

- `Merged_Phylum_raw_abund.csv`
- `Merged_Family_raw_abund.csv`
- `Merged_Genus_raw_abund.csv`
- `Merged_Species_raw_abund.csv`

Metadata files (e.g., `meta.csv`) are expected to already exist; metadata regeneration is not part of the active run path.

---

## Reproducibility Notes

- R version: >= 4.1 recommended.
- Core packages are loaded via each dataset `_globalStuff.R`.
- Run scripts non-interactively from the repository with Linux-style paths.

### Example

```r
source('PPR/2_Feature_Tables_to_Rarefied_merged.R')
source('JPL/2_Feature_Tables_to_Rarefied.R')
```

---

## Contact
Jeremy_Kahsen@Rush.edu
