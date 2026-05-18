# GGW Bird Biodiversity — Reproducibility Code

Code for:

> Wang, Y., Scott, C. E. & Dallimer, M. "Biodiversity co-benefits of large-scale ecosystem restoration: Causal evidence from avian responses to Africa's Great Green Wall." *Nature Communications* (under review, 2026). Manuscript ID: NCOMMS-26-038140-T

All analysis code is publicly available at: https://github.com/eeywang/ggw-bird-biodiversity

Processed datasets used to generate figures will be archived on Zenodo upon acceptance.

---

## System Requirements

### Software dependencies

**R v4.5.1** — required packages (with versions used in analysis):

| Package | Version | Purpose |
|---------|---------|---------|
| biomod2 | 4.2-4 | Ensemble species distribution modelling |
| maxnet | 0.1.4 | MAXENT implementation for biomod2 |
| MatchIt | 4.5 | Propensity score matching |
| cobalt | 5.x | PSM balance diagnostics |
| weights | 1.0.x | Weighted regression for doubly robust ATT |
| terra | 1.7+ | Raster data processing |
| sf | 1.0+ | Vector spatial data |
| exactextractr | 0.10+ | Zonal statistics extraction |
| tidyverse | 2.0+ | Data wrangling |
| data.table | 1.15+ | Fast data manipulation |
| spThin | 0.2+ | Spatial thinning of occurrence records |
| CoordinateCleaner | 3.0+ | Occurrence record quality control |
| rgbif | 3.7+ | GBIF data download |
| usdm | 2.1+ | Variance inflation factor screening |

**Python 3.11** — required packages:

| Package | Version | Purpose |
|---------|---------|---------|
| matplotlib | 3.8 | Figure generation |
| pandas | 2.x | Data manipulation |
| numpy | 1.x | Numerical computation |
| scipy | 1.x | Statistical functions |
| PyMuPDF (fitz) | 1.23+ | PDF figure assembly |

### Operating systems tested on

- Windows 10/11 (local scripts: causal inference, figures)
- Linux (Ubuntu 22.04) — University of Leeds Aire HPC cluster, SLURM scheduler (SDM and scenario decomposition steps)

### Non-standard hardware

The **full SDM pipeline** (Stage 3–4) requires access to a high-performance computing (HPC) cluster. A SLURM submission script is provided (`03_species_distribution_models/submit_example.sh`).

The **causal inference and figure generation** stages (Stages 5–6) can be run on a standard desktop or laptop computer.

---

## Installation Guide

### R environment

Install required R packages from CRAN:

```r
install.packages(c(
  "biomod2", "maxnet", "MatchIt", "cobalt", "weights",
  "terra", "sf", "exactextractr", "tidyverse", "data.table",
  "spThin", "CoordinateCleaner", "rgbif", "usdm"
))
```

**Note:** `biomod2 v4.2-4` was used in the analysis. Later versions may have different function signatures. To install the specific version:

```r
remotes::install_version("biomod2", version = "4.2-4")
```

### Python environment

```bash
pip install matplotlib==3.8 pandas numpy scipy PyMuPDF openpyxl
```

Or using conda:

```bash
conda create -n ggw-birds python=3.11
conda activate ggw-birds
pip install matplotlib==3.8 pandas numpy scipy PyMuPDF openpyxl
```

### Typical installation time

Approximately **10–20 minutes** on a standard desktop computer with a stable internet connection.

---

## Repository Structure

```
Github_coding/
├── main.R                              # Master orchestration script
├── 00_setup.R                          # Paths, parameters, shared species lists
├── functions.R                         # Shared helper functions
├── manual_review_decisions.R           # Hand-curated species decisions (3 countries)
│
├── 01_data_acquisition/                # GBIF download + eBird merge
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
│   └── submit_example.sh              # Example SLURM submission script
│
├── 04_scenario_decomposition/          # HPC: counterfactual scenarios (S0–S3)
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
└── data_processed/                     # Processed outputs (to be archived on Zenodo)
```

---

## Pipeline Overview

| Stage | Folder | Description | Platform |
|-------|--------|-------------|----------|
| 1 | `01_data_acquisition/` | Download GBIF + standardise eBird for 11 GGW countries | Desktop |
| 2 | `02_data_cleaning/` | 7-step coordinate cleaning + AVONET migration filtering | Desktop |
| 3 | `03_species_distribution_models/` | Ensemble SDMs (GLM, GAM, RF, MAXNET, GBM) | **HPC required** |
| 4 | `04_scenario_decomposition/` | 2×2 climate × vegetation counterfactual scenarios | **HPC required** |
| 5 | `05_causal_inference/` | PSM 1:1 matching + Doubly Robust ATT estimation | Desktop |
| 6 | `06_figures/` | Publication-quality figures | Desktop |

---

## Demo: Reproducing Causal Inference Results

The causal inference stage (Stage 5) can be run on a standard desktop computer using the processed outputs provided in `data_processed/` (to be made available via Zenodo upon acceptance). This provides the fastest route to reproducing the main results in the paper.

### Step-by-step

1. Clone this repository:
```bash
git clone https://github.com/eeywang/ggw-bird-biodiversity.git
cd ggw-bird-biodiversity
```

2. Open R and set your working directory to the repository root. Run `00_setup.R` to configure paths.

3. Run the causal inference pipeline for each country:
```r
source("05_causal_inference/03_psm_matching.R")          # PSM matching
source("05_causal_inference/05_doubly_robust_ATT.R")     # Overall ATT
source("05_causal_inference/06_scenario_DR_ATT.R")       # Pathway decomposition
source("05_causal_inference/07_guild_DR_ATT_nigeria.R")  # Guild-level (Nigeria)
source("05_causal_inference/07_guild_DR_ATT_senegal.R")  # Guild-level (Senegal)
source("05_causal_inference/07_guild_DR_ATT_ethiopia.R") # Guild-level (Ethiopia)
```

### Expected output

- PSM matched sample with covariate balance diagnostics (Love plots, SMD tables)
- Doubly robust ATT estimates for overall richness (Senegal: −3.30; Nigeria: +2.19; Ethiopia: −1.59)
- Pathway-specific ATTs (climate, vegetation, interaction components)
- Guild-specific ATT estimates (woodland, open-habitat, wetland)

These correspond to Figures 3 and 5 and Extended Data Table 1 in the manuscript.

### Expected run time (Stage 5 only, standard desktop)

Approximately **5–10 minutes** for all three countries.

### Full pipeline run time (Stages 3–4, HPC)

Approximately **24–48 hours per country** on a 16-core HPC node (University of Leeds Aire cluster, SLURM). Total HPC compute time for all three countries: ~100 core-hours.

---

## Key Methods Summary

- **Occurrence threshold**: ≥50 unique records after 5-km spatial thinning for SDM fitting
- **SDMs**: 5 algorithms (GLM, GAM, RF, MAXNET, GBM), ensemble via weighted mean + committee averaging, TSS ≥ 0.4 retention threshold; 3-fold cross-validation × 2 repetitions = 30 sub-models per species
- **Causal inference**: 1:1 nearest-neighbour PSM + Doubly Robust ATT (Bang & Robins 2005), controlling for climate trends (Δbio terms)
- **Scenarios**: S0 (baseline) / S1 (climate only) / S2 (vegetation only) / S3 (observed), crossed climate × vegetation conditions

---

## Data Availability

Bird occurrence records are publicly available from:
- **GBIF**: https://www.gbif.org (accessed via `rgbif` package)
- **eBird**: https://ebird.org (Cornell Lab of Ornithology, custom data request)

Environmental data sources:
- TerraClimate: https://www.climatologylab.org/terraclimate.html
- CHELSA v2.1: https://chelsa-climate.org
- MODIS NDVI/land cover: https://lpdaac.usgs.gov
- SRTM: https://srtm.csi.cgiar.org
- WorldPop: https://www.worldpop.org
- Copernicus Global Land Service: https://lcviewer.vito.be
- GGW boundaries: Pan-African Agency of the Great Green Wall (PAGGW)

Processed outputs (administrative-unit-level richness, matched samples, scenario values) will be archived on **Zenodo** upon acceptance.

---

## Citation

> Wang, Y., Scott, C. E. & Dallimer, M. Biodiversity co-benefits of large-scale ecosystem restoration: Causal evidence from avian responses to Africa's Great Green Wall. *Nature Communications* (under review, 2026).
