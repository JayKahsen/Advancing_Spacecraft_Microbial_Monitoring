# unneeded/

This folder temporarily stores files removed from active workflows during cleanup.

## Why these files were moved

- Duplicate, exploratory, or copy scripts that are not part of the current figure-rebuild path.
- Office lock/temp artifacts (e.g., files prefixed with `~$`).
- Legacy metadata regeneration scripts (`1_1_create_meta.R`, `make_meta_data.R`) retained for reference because existing `meta` files are already available.
- Non-canonical presentation and ad hoc output assets.

## Feature-table input policy

For `2_Feature_Tables_to_Rarefied*.R`, only these canonical taxonomic inputs remain in active `data_tables/original/` directories:

- `Merged_Phylum_raw_abund.csv`
- `Merged_Family_raw_abund.csv`
- `Merged_Genus_raw_abund.csv`
- `Merged_Species_raw_abund.csv`

Everything else from those `original/` folders was moved here for safekeeping.

- Legacy JPL directory content moved to `unneeded/JPL_legacy_active/` after promoting the former `Clipper/` dataset to active `JPL/`.
