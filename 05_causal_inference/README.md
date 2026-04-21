# 05_causal_inference — PSM + Doubly Robust ATT

## Pipeline

```
01  Prepare covariates          01a (Senegal), 01b (Ethiopia), embedded in 02b (Nigeria)
02  Climate trends (Δbio)       02a (Senegal+Ethiopia), 02b (Nigeria complete pipeline)
03  PSM 1:1 matching            03 — all 3 countries (set COUNTRY at top)
04  Scenario extraction         04 — extract S1/S2 for matched units
05  Doubly Robust ATT           05 — richness change ATT (controlling Δbio)
06  Scenario DR ATT             06 — decomposition (climate/vegetation/interaction)
07  Guild-level DR ATT          07_* — one file per country
```

## How to Run

Set `COUNTRY <- "Nigeria"` (or `"Senegal"` / `"Ethiopia"`) at the top of each script, then run Steps 03 → 04 → 05 → 06 in order.

## Files

| File | Countries | Description |
|------|-----------|-------------|
| **Covariate preparation** | | |
| `01a_prepare_covariates_senegal.R` | SEN | Extract covariates at commune (adm4) level |
| `01b_prepare_covariates_ethiopia.R` | ETH | Extract covariates at woreda (adm3) level |
| **Climate trends** | | |
| `02a_climate_trends_senegal_ethiopia.R` | SEN+ETH | Δbio trend rasters (2000–2024) |
| `02b_full_pipeline_nigeria.R` | NGA | Complete: data download + climate trends + covariate extraction + PSM + DR |
| **Modular pipeline (Step 3–6)** | | |
| `03_psm_matching.R` | **ALL 3** | PSM 1:1 nearest-neighbour matching |
| `04_scenario_extraction.R` | **ALL 3** | Extract S1/S2 scenario richness for matched admin units |
| `05_doubly_robust_ATT.R` | **ALL 3** | DR ATT: richness change controlling Δbio climate trends |
| `06_scenario_DR_ATT.R` | **ALL 3** | DR ATT: scenario decomposition (climate + vegetation + interaction) |
| **Guild analysis** | | |
| `07_guild_DR_ATT_nigeria.R` | NGA | Guild-specific DR ATT (Woodland / Open-habitat / Wetland) |
| `07_guild_DR_ATT_senegal.R` | SEN | Guild-specific DR ATT |
| `07_guild_DR_ATT_ethiopia.R` | ETH | Guild-specific DR ATT |

## Country-Specific Parameters

### PSM Covariates

| Variable | Nigeria | Senegal | Ethiopia |
|----------|:---:|:---:|:---:|
| temp_mean | ✓ | ✓ | ✓ |
| prec_mean | ✓ | ✓ | ✓ |
| prec_cv | ✓ | | ✓ |
| elev_mean | | ✓ | |
| roughness_mean | ✓ | | |
| slope_mean | | ✓ | ✓ |
| ndvi_mean | ✓ | ✓ | ✓ |
| log_pop | ✓ | ✓ | ✓ |
| lulc_cropland | ✓ | ✓ | ✓ |
| lulc_diversity | ✓ | ✓ | ✓ |

### DR Climate Trend Controls

| Country | Δbio variables |
|---------|---------------|
| Nigeria | Δbio1, Δbio12 |
| Senegal | Δbio2, Δbio4, Δbio13, Δbio18 |
| Ethiopia | Δbio3, Δbio4, Δbio7, Δbio9, Δbio13, Δbio19 |

### Administrative Units

| Country | Level | ID column | GGW units |
|---------|-------|-----------|-----------|
| Nigeria | LGA (adm2) | GID_2 | 262 |
| Senegal | Commune (adm4) | GID_4 | 76 |
| Ethiopia | Woreda (adm3) | GID_3 | 87 |

## Note on Nigeria

Nigeria's covariate extraction code is embedded in `02b_full_pipeline_nigeria.R`
(步骤3 + 步骤5). This was the original monolithic development script; the modular
Step 03–06 approach was adopted later for Senegal/Ethiopia and then Nigeria was
added to those scripts for consistency.
