# =============================================================================
# Step 2: - PSM(50)

# : Step1{COUNTRY}_1to1_matched_data_50threshold.csv


#   1. 50S1/S2 
#   2. Step1 PSM
#   3. 


#   - {COUNTRY}_matched_with_scenarios_50threshold.csv (+)
# =============================================================================

library(sf)
library(terra)
library(dplyr)
library(exactextractr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║     Step 2:  - PSM (50)       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# - (Step1)
# =============================================================================

COUNTRY <- "Nigeria"  # ← : "Nigeria", "Senegal", "Ethiopia"

if(COUNTRY == "Senegal") {
  
  adm_file     <- "E:/11.17progress/study_area/SEN_adm4.shp"
  id_col       <- "GID_4"
  scenario_dir <- "E:/2026.1.8_biomod2/SEN_scenario_analysis/ggw_senegal/scenario_analysis_50threshold"
  matched_file <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold/Senegal_1to1_matched_data_50threshold.csv"
  output_dir   <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold"
  
} else if(COUNTRY == "Ethiopia") {
  
  adm_file     <- "E:/11.17progress/study_area/ETH_adm3.shp"
  id_col       <- "GID_3"
  scenario_dir <- "E:/2026.1.8_biomod2/ETH_scenario_analysis/ggw_ethiopia/scenario_analysis_50threshold"
  matched_file <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold/Ethiopia_1to1_matched_data_50threshold.csv"
  output_dir   <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold"

} else if(COUNTRY == "Nigeria") {
  
  id_col       <- "GID_2"
  matched_file <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold/Nigeria_1to1_matched_data_50threshold.csv"
  scenario_dir <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition_with_trend"
  
  output_dir   <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold"
}

cat(sprintf(": %s\n", COUNTRY))

# =============================================================================
# 2.1 
# =============================================================================

if(!file.exists(matched_file)) {
  stop(sprintf(": PSM: %s\n Step1_PSM_matching_50threshold.R", matched_file))
}
if(!dir.exists(scenario_dir)) {
  stop(sprintf(": : %s", scenario_dir))
}

# Step1
matched_data <- read.csv(matched_file)
cat(sprintf("  PSM: %d \n", nrow(matched_data)))
cat(sprintf("    Treatment: %d, Control: %d\n",
            sum(matched_data$treatment == 1), sum(matched_data$treatment == 0)))

# =============================================================================
# 2.2 50
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  2.2 50                                        \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")


adm <- st_read(adm_file, quiet = TRUE)
cat(sprintf("  : %d\n", nrow(adm)))


ref_rast <- rast(file.path(scenario_dir, "richness_S1_binary.tif"))
adm_proj <- st_transform(adm, crs(ref_rast))

# ----- S1/S2 -----
cat("  【 S1/S2】\n")

adm$S1_binary <- exact_extract(rast(file.path(scenario_dir, "richness_S1_binary.tif")), adm_proj, 'mean')
adm$S1_prob   <- exact_extract(rast(file.path(scenario_dir, "richness_S1_prob.tif")), adm_proj, 'mean')
adm$S1_tss    <- exact_extract(rast(file.path(scenario_dir, "richness_S1_tss.tif")), adm_proj, 'mean')
adm$S2_binary <- exact_extract(rast(file.path(scenario_dir, "richness_S2_binary.tif")), adm_proj, 'mean')
adm$S2_prob   <- exact_extract(rast(file.path(scenario_dir, "richness_S2_prob.tif")), adm_proj, 'mean')
adm$S2_tss    <- exact_extract(rast(file.path(scenario_dir, "richness_S2_tss.tif")), adm_proj, 'mean')
cat("    ✓ S1/S2\n")

# ----- -----
cat("  【】\n")

# Binary
adm$binary_climate     <- exact_extract(rast(file.path(scenario_dir, "effect_binary_climate.tif")), adm_proj, 'mean')
adm$binary_vegetation  <- exact_extract(rast(file.path(scenario_dir, "effect_binary_vegetation.tif")), adm_proj, 'mean')
adm$binary_interaction <- exact_extract(rast(file.path(scenario_dir, "effect_binary_interaction.tif")), adm_proj, 'mean')
adm$binary_total       <- exact_extract(rast(file.path(scenario_dir, "effect_binary_total.tif")), adm_proj, 'mean')

# Probability
adm$prob_climate     <- exact_extract(rast(file.path(scenario_dir, "effect_probability_climate.tif")), adm_proj, 'mean')
adm$prob_vegetation  <- exact_extract(rast(file.path(scenario_dir, "effect_probability_vegetation.tif")), adm_proj, 'mean')
adm$prob_interaction <- exact_extract(rast(file.path(scenario_dir, "effect_probability_interaction.tif")), adm_proj, 'mean')
adm$prob_total       <- exact_extract(rast(file.path(scenario_dir, "effect_probability_total.tif")), adm_proj, 'mean')

# TSS-weighted
adm$tss_climate     <- exact_extract(rast(file.path(scenario_dir, "effect_tss_weighted_climate.tif")), adm_proj, 'mean')
adm$tss_vegetation  <- exact_extract(rast(file.path(scenario_dir, "effect_tss_weighted_vegetation.tif")), adm_proj, 'mean')
adm$tss_interaction <- exact_extract(rast(file.path(scenario_dir, "effect_tss_weighted_interaction.tif")), adm_proj, 'mean')
adm$tss_total       <- exact_extract(rast(file.path(scenario_dir, "effect_tss_weighted_total.tif")), adm_proj, 'mean')
cat("    ✓ \n")

# =============================================================================
# 2.3 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  2.3 PSM                                       \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")


adm_df <- st_drop_geometry(adm)

scenario_cols <- c(id_col,
                   "S1_binary", "S1_prob", "S1_tss",
                   "S2_binary", "S2_prob", "S2_tss",
                   "binary_climate", "binary_vegetation", "binary_interaction", "binary_total",
                   "prob_climate", "prob_vegetation", "prob_interaction", "prob_total",
                   "tss_climate", "tss_vegetation", "tss_interaction", "tss_total")
scenario_subset <- adm_df[, intersect(scenario_cols, names(adm_df))]

# ID
scenario_subset <- scenario_subset[!duplicated(scenario_subset[[id_col]]), ]


matched_with_scenarios <- merge(matched_data, scenario_subset, by = id_col, all.x = TRUE)

cat(sprintf("  : %d \n", nrow(matched_with_scenarios)))


n_with_scenarios <- sum(!is.na(matched_with_scenarios$binary_climate))
cat(sprintf("  : %d / %d (%.1f%%)\n",
            n_with_scenarios, nrow(matched_with_scenarios),
            100 * n_with_scenarios / nrow(matched_with_scenarios)))

# =============================================================================
# 2.4 
# =============================================================================

output_file <- file.path(output_dir, paste0(COUNTRY, "_matched_with_scenarios_50threshold.csv"))
write.csv(matched_with_scenarios, output_file, row.names = FALSE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Step 2 ！                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  : %s\n", output_file))
cat(sprintf("  : %d  × %d \n", nrow(matched_with_scenarios), ncol(matched_with_scenarios)))


cat("\n  【 - Binary】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
treat <- matched_with_scenarios[matched_with_scenarios$treatment == 1, ]
ctrl  <- matched_with_scenarios[matched_with_scenarios$treatment == 0, ]

for(ef in c("binary_climate", "binary_vegetation", "binary_interaction", "binary_total")) {
  if(ef %in% names(matched_with_scenarios)) {
    t_mean <- mean(treat[[ef]], na.rm = TRUE)
    c_mean <- mean(ctrl[[ef]], na.rm = TRUE)
    cat(sprintf("  %-25s GGW=%.2f  Non-GGW=%.2f  Diff=%.2f\n", ef, t_mean, c_mean, t_mean - c_mean))
  }
}

cat(sprintf("\n:  Step3_doubly_robust_ATT_50threshold.R\n"))
