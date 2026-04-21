# 06_figures — Figure generation code

## Main Figures

| Figure | Script | Language | Description |
|--------|--------|----------|-------------|
| Fig 2a-d | `Fig2ad_scenario_maps.R` | R | Scenario decomposition spatial maps (S0-S3), 180mm |
| Fig 2e | `Fig2e_total_boxplot.py` | Python | Total richness change boxplot (GGW vs Non-GGW) |
| Fig 3 | `Fig3_causal_effects.py` | Python | Causal effects: ATT + dumbbell + pathway bars (3x2 grid) |
| Fig 4a-c | `Fig4abc_guild_maps.R` | R | Guild spatial change maps (Woodland / Open-habitat / Wetland) |
| Fig 4d | `Fig4d_guild_boxplot.py` | Python | Guild richness change boxplot |
| Fig 5 | `Fig5_guild_ATT.py` | Python | Guild-specific DR ATT dot plot |
| Fig 4 | `Fig4_assembly.py` | Python | PyMuPDF assembly (panels a-d into combined PDF) |

## Extended Data & Supplementary

| Figure | Script | Description |
|--------|--------|-------------|
| ED Fig 1e | `ED_Fig1e_pathway_boxplot.py` | Scenario pathway decomposition boxplot (3 facets) |
| ED Fig 2 | `ED_Fig2_SDM_performance.py` | SDM performance ridgeline (TSS + AUC) |
| SI Fig S1 | `SI_FigS1_PSM_diagnostics.py` | PSM matching diagnostics (Love plot + PS density) |

## Style Standards (NEE)

- Font: Helvetica / Arial, regular weight throughout
- Panel labels (a, b, c...): 8-9 pt, **bold** (only exception to regular weight)
- Axis labels: 7-8 pt | Tick labels: 6-7 pt | Minimum anywhere: 5 pt
- Full box borders, inward-facing ticks
- Two-column width: 180 mm | `pdf.fonttype = 42` for Illustrator editability
- Country order: Senegal, Nigeria, Ethiopia (west-to-east)
- Country colours: Senegal `#E91E63` (pink), Nigeria `#009688` (teal), Ethiopia `#FF9800` (orange)
