# data_processed — Intermediate data for figure reproduction

Minimum data needed to reproduce all manuscript figures without re-running
the full SDM pipeline (~3,500 CPU-hours on HPC).

## Contents

| Folder | Files | Size | Description |
|--------|-------|------|-------------|
| `occurrence/` | `{Country}_final_residents.csv` x 3 | ~62 MB | Cleaned resident bird occurrence records |
| `species_lists/` | `Supplementary_Table_S3_species_classification.csv` | 57 KB | Species list with habitat guild + trophic niche + migration status |
| `model_evaluation/` | `{Country}_all_species_summary.csv` x 3 | ~150 KB | Per-species ensemble SDM performance (mean TSS, AUC) |
| `psm_matched/` | `{Country}_psm_matched_complete.csv` x 3 | ~1 MB | PSM matched admin units with covariates, richness, scenarios, climate trends, and guild richness |
| `results_tables/` | `DR_ATT_results.csv`, `guild_DR_ATT_results.csv` | ~10 KB | Final doubly robust ATT estimates |

Total: ~63 MB

## Figure to data mapping

| Figure | Input data |
|--------|-----------|
| Fig 2e (total boxplot) | `psm_matched/*_psm_matched_complete.csv` → `prob_total` column |
| Fig 3 (causal effects) | `results_tables/DR_ATT_results.csv` |
| Fig 4d (guild boxplot) | `psm_matched/*_psm_matched_complete.csv` → `Woodland/Openhabitat/Wetland_change` columns |
| Fig 4e (guild ATT) | `results_tables/guild_DR_ATT_results.csv` |
| ED Fig 1e (pathway boxplot) | `psm_matched/*_psm_matched_complete.csv` → `prob_climate/vegetation/interaction` columns |
| ED Fig 2 (SDM performance) | `model_evaluation/*_all_species_summary.csv` |
| SI Fig S1 (PSM diagnostics) | `psm_matched/*_psm_matched_complete.csv` → covariates + treatment columns |

## Column structure of `*_psm_matched_complete.csv`

Each row = one administrative unit (PSM 1:1 matched).

| Column group | Columns | Description |
|-------------|---------|-------------|
| ID | `GID_*`, `NAME_*`, `treatment` | Admin unit ID + GGW treatment assignment |
| PSM covariates | `temp_mean`, `prec_mean`, `ndvi_mean`, ... (8 vars) | Pre-treatment conditions (2000-2006) |
| Richness S0/S3 | `prob_before`, `prob_after`, `prob_change`, ... | Expected species richness (binary/prob/TSS-weighted) |
| Scenario decomposition | `prob_climate`, `prob_vegetation`, `prob_interaction`, `prob_total` | Counterfactual effect decomposition |
| Climate trends | `bio*_trend` (2-6 per country) | Bioclimatic trend controls for DR ATT |
| Guild richness | `richness_prob_Woodland_S0/S3`, `Woodland_change`, ... | Guild-level richness (3 guilds x S0/S3/change) |
| Matching info | `distance`, `weights`, `subclass` | PSM diagnostics |

## Data sources

- Occurrence records: GBIF (https://www.gbif.org) + eBird (Cornell Lab of Ornithology)
- Migration status: AVONET (Tobias et al. 2022, *Ecology Letters*, doi:10.1111/ele.13898)
- Environmental variables: WorldClim, MODIS NDVI, SRTM, WorldPop, ESA CCI Land Cover
