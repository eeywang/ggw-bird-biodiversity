# =============================================================================
# Step 1: PSM 1:1 (50)


#   1. + Treatment
#   2. 50(S0=before, S3=after)
#   3. PSM 1:1
#   4. 


#   - {COUNTRY}_1to1_matched_data_50threshold.csv  ()
#   - {COUNTRY}_50threshold_balance_loveplot.png    (Love plot)
#   - {COUNTRY}_50threshold_propensity_score.png    ()
# =============================================================================

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(exactextractr)
library(MatchIt)
library(cobalt)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║        Step 1: PSM 1:1  (50)                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# - 
# =============================================================================

COUNTRY <- "Nigeria"  # ← : "Nigeria", "Senegal", "Ethiopia"

if(COUNTRY == "Senegal") {
  
  adm_file   <- "E:/11.17progress/study_area/SEN_adm4.shp"
  id_col     <- "GID_4"
  name_col   <- "NAME_4"
  cov_file   <- "E:/2026.1.8_biomod2/SEN_covariates/Senegal_Adm4_covariates_clean.csv"
  richness_dir <- "E:/2026.1.8_biomod2/senegal_results/richness/richness_50threshold"
  
  psm_covars <- c("temp_mean", "prec_mean", "elev_mean", "slope_mean",
                  "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")
  
  ggw_ids <- c(
    "SEN.8.1.1.6_1", "SEN.8.1.1.7_1", "SEN.8.1.2.1_1", "SEN.8.1.2.2_1", "SEN.8.1.2.3_1",
    "SEN.8.1.2.4_1", "SEN.8.1.2.5_1", "SEN.8.1.2.6_1", "SEN.8.1.2.7_1", "SEN.8.1.2.8_1",
    "SEN.8.1.2.9_1", "SEN.8.1.2.10_1", "SEN.8.2.1.1_1", "SEN.8.2.2.1_1", "SEN.8.2.2.2_1",
    "SEN.8.2.2.3_1", "SEN.8.2.3.2_1", "SEN.8.2.3.5_1", "SEN.8.2.4.1_1", "SEN.8.2.4.2_1",
    "SEN.8.2.4.3_1", "SEN.8.2.4.4_1", "SEN.8.3.1.1_1", "SEN.8.3.1.2_1", "SEN.8.3.1.3_1",
    "SEN.8.3.1.4_1", "SEN.8.3.2.1_1", "SEN.8.3.2.2_1", "SEN.8.3.2.3_1", "SEN.8.3.2.4_1",
    "SEN.8.3.3.1_1", "SEN.8.3.3.2_1", "SEN.8.3.3.3_1", "SEN.8.3.3.4_1", "SEN.8.3.4.1_1",
    "SEN.8.3.4.2_1", "SEN.8.3.4.3_1", "SEN.9.1.1.1_1", "SEN.9.1.1.2_1", "SEN.9.1.1.3_1",
    "SEN.9.1.2.1_1", "SEN.9.1.2.2_1", "SEN.9.2.1.1_1", "SEN.9.2.1.2_1", "SEN.9.2.1.3_1",
    "SEN.9.2.2.1_1", "SEN.9.2.2.2_1", "SEN.9.2.2.3_1", "SEN.9.3.1.1_1", "SEN.9.3.1.2_1",
    "SEN.9.3.1.3_1", "SEN.10.1.1.1_1", "SEN.10.1.1.2_1", "SEN.10.1.2.1_1", "SEN.10.1.2.2_1",
    "SEN.10.1.2.3_1", "SEN.10.2.1.1_1", "SEN.10.2.1.2_1", "SEN.10.2.1.3_1", "SEN.10.2.2.1_1",
    "SEN.10.2.2.2_1", "SEN.10.2.2.3_1", "SEN.10.2.3.1_1", "SEN.10.2.3.2_1", "SEN.10.2.4.1_1",
    "SEN.10.2.4.2_1", "SEN.10.3.1.1_1", "SEN.10.3.1.2_1", "SEN.10.3.1.3_1", "SEN.12.1.1.1_1",
    "SEN.12.1.1.2_1", "SEN.12.1.2.1_1", "SEN.12.1.2.2_1", "SEN.12.1.3.1_1", "SEN.12.1.3.2_1",
    "SEN.12.1.3.3_1", "SEN.12.2.2.1_1", "SEN.12.2.2.3_1", "SEN.12.2.2.4_1", "SEN.12.2.4.2_1"
  )
  
  output_dir <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold"
  
} else if(COUNTRY == "Ethiopia") {
  
  adm_file   <- "E:/11.17progress/study_area/ETH_adm3.shp"
  id_col     <- "GID_3"
  name_col   <- "NAME_3"
  cov_file   <- "E:/2026.1.8_biomod2/ETH_covariates/Ethiopia_Adm3_covariates_updated.csv"
  richness_dir <- "E:/2026.1.8_biomod2/ethiopia_results/richness/richness_50threshold"
  
  psm_covars <- c("prec_mean", "prec_cv", "temp_mean", "slope_mean",
                  "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")
  
  ggw_ids <- c(
    "ETH.2.1.1_1", "ETH.2.1.2_1", "ETH.2.1.3_1", "ETH.2.1.4_1", "ETH.2.1.5_1",
    "ETH.2.1.6_1", "ETH.2.1.7_1", "ETH.2.1.8_1", "ETH.2.2.1_1", "ETH.2.2.2_1",
    "ETH.2.2.3_1", "ETH.2.2.4_1", "ETH.2.2.6_1", "ETH.2.2.8_1", "ETH.2.4.1_1",
    "ETH.2.4.2_1", "ETH.2.4.3_1", "ETH.2.4.4_1", "ETH.2.4.5_1", "ETH.3.1.3_1",
    "ETH.3.1.7_1", "ETH.3.4.3_1", "ETH.3.4.4_1", "ETH.3.4.6_1", "ETH.3.4.7_1",
    "ETH.3.5.2_1", "ETH.3.5.18_1", "ETH.3.6.1_1", "ETH.3.6.3_1", "ETH.3.6.10_1",
    "ETH.3.9.2_1", "ETH.3.10.1_1", "ETH.3.10.2_1", "ETH.3.10.3_1", "ETH.3.10.4_1",
    "ETH.3.10.5_1", "ETH.3.10.6_1", "ETH.3.10.7_1", "ETH.3.10.8_1", "ETH.3.10.9_1",
    "ETH.3.10.10_1", "ETH.3.10.11_1", "ETH.3.10.12_1", "ETH.3.10.13_1", "ETH.3.10.14_1",
    "ETH.3.10.15_1", "ETH.3.10.16_1", "ETH.3.10.17_1", "ETH.3.10.18_1", "ETH.3.10.19_1",
    "ETH.3.10.20_1", "ETH.3.11.1_1", "ETH.3.11.3_1", "ETH.3.11.4_1", "ETH.3.11.5_1",
    "ETH.3.11.6_1", "ETH.3.11.7_1", "ETH.3.11.8_1", "ETH.3.11.9_1", "ETH.3.11.10_1",
    "ETH.3.12.1_1", "ETH.3.12.2_1", "ETH.3.12.3_1", "ETH.3.12.4_1", "ETH.3.12.5_1",
    "ETH.3.12.6_1", "ETH.4.3.2_1", "ETH.4.3.4_1", "ETH.4.3.5_1", "ETH.11.1.1_1",
    "ETH.11.1.2_1", "ETH.11.1.3_1", "ETH.11.1.4_1", "ETH.11.1.5_1", "ETH.11.1.6_1",
    "ETH.11.1.7_1", "ETH.11.1.8_1", "ETH.11.2.2_1", "ETH.11.2.3_1", "ETH.11.2.6_1",
    "ETH.11.2.9_1", "ETH.11.3.1_1", "ETH.11.3.2_1", "ETH.11.3.3_1", "ETH.11.5.1_1",
    "ETH.11.5.2_1", "ETH.11.5.4_1", "ETH.11.5.6_1"
  )
  
  output_dir <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold"

} else if(COUNTRY == "Nigeria") {
  
  adm_file   <- "E:/11.17progress/study_area/Nigeria_adm2.shp"
  id_col     <- "GID_2"
  name_col   <- "NAME_2"
  cov_file   <- "E:/2026.1.8_biomod2/NGA_covariates/Nigeria_LGA_covariates_complete.csv"
  richness_dir <- "E:/12.1progress_biomod2/result/richness"
  
  psm_covars <- c("temp_mean", "prec_mean", "prec_cv", "roughness_mean",
                  "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")
  
  ggw_ids <- c(
    "NGA.2.1_1", "NGA.2.2_1", "NGA.2.3_1", "NGA.2.4_1", "NGA.2.5_1", "NGA.2.6_1",
  "NGA.2.7_1", "NGA.2.8_1", "NGA.2.9_1", "NGA.2.10_1", "NGA.2.11_1", "NGA.2.12_1",
  "NGA.2.13_1", "NGA.2.14_1", "NGA.2.15_1", "NGA.2.16_1", "NGA.2.17_1", "NGA.2.18_1",
  "NGA.2.19_1", "NGA.2.20_1", "NGA.2.21_1", "NGA.5.1_1", "NGA.5.2_1", "NGA.5.3_1",
  "NGA.5.4_1", "NGA.5.5_1", "NGA.5.6_1", "NGA.5.7_1", "NGA.5.8_1", "NGA.5.9_1",
  "NGA.5.10_1", "NGA.5.11_1", "NGA.5.12_1", "NGA.5.13_1", "NGA.5.14_1", "NGA.5.15_1",
  "NGA.5.16_1", "NGA.5.17_1", "NGA.5.18_1", "NGA.5.19_1", "NGA.5.20_1", "NGA.8.1_1",
  "NGA.8.2_1", "NGA.8.3_1", "NGA.8.4_1", "NGA.8.5_1", "NGA.8.6_1", "NGA.8.7_1",
  "NGA.8.8_1", "NGA.8.9_1", "NGA.8.10_1", "NGA.8.11_1", "NGA.8.12_1", "NGA.8.13_1",
  "NGA.8.14_1", "NGA.8.15_1", "NGA.8.16_1", "NGA.8.17_1", "NGA.8.18_1", "NGA.8.19_1",
  "NGA.8.20_1", "NGA.8.21_1", "NGA.8.22_1", "NGA.8.23_1", "NGA.8.24_1", "NGA.8.25_1",
  "NGA.8.26_1", "NGA.8.27_1", "NGA.8.28_1", "NGA.16.1_1", "NGA.16.2_1", "NGA.16.3_1",
  "NGA.16.4_1", "NGA.16.5_1", "NGA.16.6_1", "NGA.16.7_1", "NGA.16.8_1", "NGA.16.9_1",
  "NGA.16.10_1", "NGA.16.11_1", "NGA.18.1_1", "NGA.18.2_1", "NGA.18.3_1", "NGA.18.4_1",
  "NGA.18.5_1", "NGA.18.6_1", "NGA.18.7_1", "NGA.18.8_1", "NGA.18.9_1", "NGA.18.10_1",
  "NGA.18.11_1", "NGA.18.12_1", "NGA.18.13_1", "NGA.18.14_1", "NGA.18.15_1", "NGA.18.16_1",
  "NGA.18.17_1", "NGA.18.18_1", "NGA.18.19_1", "NGA.18.20_1", "NGA.18.21_1", "NGA.18.22_1",
  "NGA.18.23_1", "NGA.18.24_1", "NGA.18.25_1", "NGA.18.26_1", "NGA.18.27_1", "NGA.19.5_1",
  "NGA.19.17_1", "NGA.20.1_1", "NGA.20.2_1", "NGA.20.3_1", "NGA.20.4_1", "NGA.20.5_1",
  "NGA.20.6_1", "NGA.20.7_1", "NGA.20.8_1", "NGA.20.9_1", "NGA.20.10_1", "NGA.20.11_1",
  "NGA.20.12_1", "NGA.20.13_1", "NGA.20.14_1", "NGA.20.15_1", "NGA.20.16_1", "NGA.20.17_1",
  "NGA.20.18_1", "NGA.20.19_1", "NGA.20.20_1", "NGA.20.21_1", "NGA.20.22_1", "NGA.20.23_1",
  "NGA.20.24_1", "NGA.20.25_1", "NGA.20.26_1", "NGA.20.27_1", "NGA.20.28_1", "NGA.20.29_1",
  "NGA.20.30_1", "NGA.20.31_1", "NGA.20.32_1", "NGA.20.33_1", "NGA.20.34_1", "NGA.20.35_1",
  "NGA.20.36_1", "NGA.20.37_1", "NGA.20.38_1", "NGA.20.39_1", "NGA.20.40_1", "NGA.20.41_1",
  "NGA.20.42_1", "NGA.20.43_1", "NGA.20.44_1", "NGA.21.1_1", "NGA.21.2_1", "NGA.21.3_1",
  "NGA.21.4_1", "NGA.21.5_1", "NGA.21.6_1", "NGA.21.7_1", "NGA.21.8_1", "NGA.21.9_1",
  "NGA.21.10_1", "NGA.21.11_1", "NGA.21.12_1", "NGA.21.13_1", "NGA.21.14_1", "NGA.21.15_1",
  "NGA.21.16_1", "NGA.21.17_1", "NGA.21.18_1", "NGA.21.19_1", "NGA.21.20_1", "NGA.21.21_1",
  "NGA.21.22_1", "NGA.21.23_1", "NGA.21.24_1", "NGA.21.25_1", "NGA.21.26_1", "NGA.21.27_1",
  "NGA.21.28_1", "NGA.21.29_1", "NGA.21.30_1", "NGA.21.31_1", "NGA.21.32_1", "NGA.21.33_1",
  "NGA.21.34_1", "NGA.22.1_1", "NGA.22.2_1", "NGA.22.3_1", "NGA.22.4_1", "NGA.22.5_1",
  "NGA.22.6_1", "NGA.22.7_1", "NGA.22.8_1", "NGA.22.9_1", "NGA.22.10_1", "NGA.22.11_1",
  "NGA.22.12_1", "NGA.22.13_1", "NGA.22.14_1", "NGA.22.15_1", "NGA.22.16_1", "NGA.22.17_1",
  "NGA.22.18_1", "NGA.22.19_1", "NGA.22.20_1", "NGA.22.21_1", "NGA.27.2_1", "NGA.27.15_1",
  "NGA.27.21_1", "NGA.32.4_1", "NGA.32.7_1", "NGA.34.1_1", "NGA.34.2_1", "NGA.34.3_1",
  "NGA.34.4_1", "NGA.34.5_1", "NGA.34.6_1", "NGA.34.7_1", "NGA.34.8_1", "NGA.34.9_1",
  "NGA.34.10_1", "NGA.34.11_1", "NGA.34.12_1", "NGA.34.13_1", "NGA.34.14_1", "NGA.34.15_1",
  "NGA.34.16_1", "NGA.34.17_1", "NGA.34.18_1", "NGA.34.19_1", "NGA.34.20_1", "NGA.34.21_1",
  "NGA.34.22_1", "NGA.34.23_1", "NGA.35.2_1", "NGA.35.4_1", "NGA.35.8_1", "NGA.36.1_1",
  "NGA.36.2_1", "NGA.36.3_1", "NGA.36.4_1", "NGA.36.5_1", "NGA.36.6_1", "NGA.36.7_1",
  "NGA.36.8_1", "NGA.36.9_1", "NGA.36.10_1", "NGA.36.11_1", "NGA.36.12_1", "NGA.36.13_1",
  "NGA.36.14_1", "NGA.36.15_1", "NGA.36.16_1", "NGA.36.17_1", "NGA.37.1_1", "NGA.37.2_1",
  "NGA.37.3_1", "NGA.37.4_1", "NGA.37.5_1", "NGA.37.6_1", "NGA.37.7_1", "NGA.37.8_1",
  "NGA.37.9_1", "NGA.37.10_1", "NGA.37.11_1", "NGA.37.12_1", "NGA.37.13_1", "NGA.37.14_1")
  
  output_dir <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold"

}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf(": %s\n", COUNTRY))
cat(sprintf("ID: %s\n", id_col))
cat(sprintf("GGW: %d\n", length(ggw_ids)))

# =============================================================================
# 1.1 + Treatment
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  1.1  + Treatment                                  \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cov_data <- read.csv(cov_file)
cat(sprintf("  : %d \n", nrow(cov_data)))

# TreatmentArcGIS
cov_data$treatment_old <- cov_data$treatment
cov_data$treatment <- ifelse(cov_data[[id_col]] %in% ggw_ids, 1, 0)

cat("  Treatment:\n")
cat(sprintf("     - GGW: %d\n", sum(cov_data$treatment_old == 1)))
cat(sprintf("     - GGW: %d\n", sum(cov_data$treatment == 1)))
cat(sprintf("     - GGW: %d\n", sum(cov_data$treatment == 0)))

# =============================================================================
# 1.2 50
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  1.2 50 (S0/S3)                                 \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

if(!dir.exists(richness_dir)) stop(sprintf(": : %s", richness_dir))

adm <- st_read(adm_file, quiet = TRUE)
cat(sprintf("  : %d\n", nrow(adm)))

ref_rast <- rast(file.path(richness_dir, "richness_binary_after.tif"))
adm_proj <- st_transform(adm, crs(ref_rast))

# S0 (before) S3 (after) 
adm$bin_before  <- exact_extract(rast(file.path(richness_dir, "richness_binary_before.tif")), adm_proj, 'mean')
adm$bin_after   <- exact_extract(rast(file.path(richness_dir, "richness_binary_after.tif")), adm_proj, 'mean')
adm$prob_before <- exact_extract(rast(file.path(richness_dir, "richness_prob_before.tif")), adm_proj, 'mean')
adm$prob_after  <- exact_extract(rast(file.path(richness_dir, "richness_prob_after.tif")), adm_proj, 'mean')
adm$tss_before  <- exact_extract(rast(file.path(richness_dir, "richness_tss_before.tif")), adm_proj, 'mean')
adm$tss_after   <- exact_extract(rast(file.path(richness_dir, "richness_tss_after.tif")), adm_proj, 'mean')
cat("  ✓ S0/S3\n")


adm$bin_change  <- adm$bin_after  - adm$bin_before
adm$prob_change <- adm$prob_after - adm$prob_before
adm$tss_change  <- adm$tss_after  - adm$tss_before


adm_df <- st_drop_geometry(adm)

# =============================================================================
# 1.3 + 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  1.3                                           \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

richness_cols <- c(id_col, name_col,
                   "bin_before", "bin_after", "bin_change",
                   "prob_before", "prob_after", "prob_change",
                   "tss_before", "tss_after", "tss_change")
adm_subset <- adm_df[, intersect(richness_cols, names(adm_df))]

full_data <- merge(cov_data, adm_subset, by = id_col, all.x = TRUE)
cat(sprintf("  : %d \n", nrow(full_data)))
cat(sprintf("    GGW (Treatment=1): %d\n", sum(full_data$treatment == 1)))
cat(sprintf("    GGW (Treatment=0): %d\n", sum(full_data$treatment == 0)))

# =============================================================================
# 1.4 PSM 1:1 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  1.4 PSM 1:1                                               \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("  PSM (%d):\n", length(psm_covars)))
for(v in psm_covars) cat(sprintf("    - %s\n", v))

outcome_vars <- c("bin_change", "prob_change", "tss_change",
                  "bin_before", "bin_after", "prob_before", "prob_after",
                  "tss_before", "tss_after")

match_vars <- c("treatment", id_col, psm_covars, outcome_vars)
match_data <- full_data[, intersect(match_vars, names(full_data))]
match_data <- na.omit(match_data)

cat(sprintf("\n  :\n"))
cat(sprintf("    Treatment (GGW): %d\n", sum(match_data$treatment == 1)))
cat(sprintf("    Control (GGW): %d\n", sum(match_data$treatment == 0)))

psm_formula <- as.formula(paste("treatment ~", paste(psm_covars, collapse = " + ")))
cat("\n  PSM:\n    ")
print(psm_formula)

cat("\n  1:1...\n")
set.seed(123)
psm_nn <- matchit(
  psm_formula,
  data = match_data,
  method = "nearest",
  distance = "glm",
  link = "logit",
  ratio = 1,
  replace = FALSE,
  estimand = "ATT"
)

cat("\n  ===  ===\n")
print(summary(psm_nn))

matched_data <- match.data(psm_nn)
cat(sprintf("\n  :\n"))
cat(sprintf("    Treatment: %d\n", sum(matched_data$treatment == 1)))
cat(sprintf("    Control: %d\n", sum(matched_data$treatment == 0)))
cat(sprintf("    : %d\n", nrow(matched_data)))

# =============================================================================
# 1.5 
# =============================================================================

cat("\n  【】\n")
bal_tab <- bal.tab(psm_nn, un = TRUE, thresholds = c(m = 0.1))
print(bal_tab)

# Love plot
png(file.path(output_dir, paste0(COUNTRY, "_50threshold_balance_loveplot.png")),
    width = 2400, height = 1600, res = 300)
love.plot(psm_nn,
          binary = "std",
          thresholds = c(m = 0.1),
          var.order = "unadjusted",
          title = sprintf("%s: Covariate Balance (50-threshold, 1:1 Matching)", COUNTRY),
          colors = c("#E74C3C", "#27AE60"))
dev.off()
cat("  ✓ Love plot\n")


png(file.path(output_dir, paste0(COUNTRY, "_50threshold_propensity_score.png")),
    width = 2000, height = 1200, res = 300)
plot(psm_nn, type = "jitter", interactive = FALSE)
dev.off()
cat("  ✓ \n")

# =============================================================================
# 1.6 
# =============================================================================

matched_file <- file.path(output_dir, paste0(COUNTRY, "_1to1_matched_data_50threshold.csv"))
write.csv(matched_data, matched_file, row.names = FALSE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Step 1 ！                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  : %s\n", matched_file))
cat(sprintf("  : %d (Treatment: %d, Control: %d)\n",
            nrow(matched_data), sum(matched_data$treatment == 1), sum(matched_data$treatment == 0)))
cat(sprintf("\n:  Step2_scenario_decomposition_50threshold.R\n"))
-