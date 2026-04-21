# GGW Bird Biodiversity — Reproducibility Code

Code for **Wang et al. — "Causal effects of Africa's Great Green Wall on avian biodiversity across Nigeria, Senegal, and Ethiopia"** (University of Leeds).

## Repository Structure

```
Github_coding/
├── main.R                              # Master orchestration script
├── 00_setup.R                          # Paths, parameters, shared species lists
├── functions.R                         # Shared helper functions
├── manual_review_decisions.R           # Hand-curated species decisions (3 countries)
│
├── 01_data_acquisition/                # GBIF download + eBird merge (11 GGW countries)
│   ├── 01_download_gbif.R
│   └── 02_standardize_merge.R
│
├── 02_data_cleaning/                   # Coordinate cleaning + migration filtering
│   ├── 01_clean_classify.R
│   └── 02_apply_review.R
│
├── 03_species_distribution_models/     # HPC: biomod2 ensemble SDMs
│   ├── 00_hpc_config.R
│   ├── 02_biomod2_modeling.R
│   ├── 03_calculate_richness.R
│   ├── 04_extract_var_importance.R
│   └── submit_example.sh
│
├── 04_scenario_decomposition/          # HPC: counterfactual scenarios
│   ├── 01_scenario_projection.R
│   ├── 02_scenario_richness.R
│   └── 03_guild_richness.R
│
├── 05_causal_inference/                # PSM + Doubly Robust ATT
│   ├── README.md
│   ├── 01a_prepare_covariates_senegal.R
│   ├── 01b_prepare_covariates_ethiopia.R
│   ├── 02a_climate_trends_senegal_ethiopia.R
│   ├── 02b_full_pipeline_nigeria.R
│   ├── 03_psm_matching.R
│   ├── 04_scenario_extraction.R
│   ├── 05_doubly_robust_ATT.R
│   ├── 06_scenario_DR_ATT.R
│   ├── 07_guild_DR_ATT_nigeria.R
│   ├── 07_guild_DR_ATT_senegal.R
│   └── 07_guild_DR_ATT_ethiopia.R
│
├── 06_figures/                         # Publication figures (R + Python)
│   ├── Fig2ad_scenario_maps.R
│   ├── Fig2e_total_boxplot.py
│   ├── Fig3_causal_effects.py
│   ├── Fig4abc_guild_maps.R
│   ├── Fig4d_guild_boxplot.py
│   ├── Fig4_assembly.py
│   ├── Fig5_guild_ATT.py
│   ├── ED_Fig1e_pathway_boxplot.py
│   ├── ED_Fig2_SDM_performance.py
│   └── SI_FigS1_PSM_diagnostics.py
│
├── environment/
│   └── 01_env_variable_selection_VIF.R
│
└── data_processed/
```

## Pipeline

| Stage | Folder | Description |
|-------|--------|-------------|
| 1 | `01_data_acquisition/` | Download GBIF + standardise eBird for 11 GGW countries |
| 2 | `02_data_cleaning/` | 7-step coordinate cleaning + AVONET migration filtering |
| 3 | `03_species_distribution_models/` | Ensemble SDMs (GLM, GAM, RF, MAXNET, GBM) on HPC |
| 4 | `04_scenario_decomposition/` | 2x2 climate x vegetation counterfactual scenarios |
| 5 | `05_causal_inference/` | PSM matching + Doubly Robust ATT estimation |
| 6 | `06_figures/` | Publication-quality figures (NEE format) |

## Key Methods

- **Occurrence threshold**: >= 50 records (2016-2024) for SDM fitting
- **SDMs**: 5 algorithms, ensemble via weighted mean + committee averaging (TSS >= 0.4)
- **Causal inference**: PSM 1:1 nearest-neighbour + Doubly Robust ATT (controlling climate trends)
- **Scenarios**: S0 (baseline) / S1 (climate only) / S2 (vegetation only) / S3 (observed)

## Requirements

R: tidyverse, data.table, readxl, sf, terra, biomod2, maxnet, MatchIt, cobalt, weights
Python: pandas, numpy, matplotlib, scipy, scikit-learn, openpyxl, PyMuPDF

## Citation

> Wang, Y. et al. (in prep). Causal effects of Africa's Great Green Wall
> initiative on avian biodiversity across Nigeria, Senegal, and Ethiopia.
