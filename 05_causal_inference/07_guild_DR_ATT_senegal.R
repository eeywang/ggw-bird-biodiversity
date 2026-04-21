# =============================================================================
# Senegal Guild-specific ATT

# : PSM + WLS
# : PSM

# Senegal: bio2_trend, bio4_trend, bio13_trend, bio18_trend


#   - guild richness(S0, S3) × 3 guilds
#   - Senegalshapefile (adm4)
#   - PSM 1:1 matched data
#   - 


#   - Senegal_guild_richness_by_admin.csv
#   - Senegal_guild_DR_ATT_results.csv
#   - Senegal_matched_with_guild_richness.csv
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(weights)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║     Senegal Guild-specific ATT                           ║\n")
cat("║     : PSM + WLS (Δbio2, Δbio4, Δbio13, Δbio18)             ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. 
# =============================================================================

COUNTRY <- "Senegal"


guild_raster_dir <- "D:/OneDrive - University of Leeds/biodiversity_manscript/results/Senegal/guild_richness"
admin_shp <- "E:/11.17progress/study_area/SEN_adm4.shp"
psm_file <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold/Senegal_1to1_matched_data_50threshold.csv"
climate_file <- "E:/2026.1.8_biomod2/senegal_results/climate_trends/Senegal_PSM_with_climate_trends.csv"


output_dir <- "D:/OneDrive - University of Leeds/biodiversity_manscript/results/Senegal/guild_richness"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Senegal
climate_vars <- c("bio2_trend", "bio4_trend", "bio13_trend", "bio18_trend")
id_col <- "GID_4"

cat(sprintf(": %s\n", COUNTRY))
cat(sprintf(": %s\n", paste(climate_vars, collapse = ", ")))
cat(sprintf("ID: %s\n", id_col))

# =============================================================================
# 2. 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  2.1                                                         \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# --- 2.1 ---
cat("  shapefile...\n")
admin <- st_read(admin_shp, quiet = TRUE)
cat(sprintf("    %d \n", nrow(admin)))

# --- 2.2 guild richness---
cat("\n  guild richness...\n")

guilds <- c("Woodland", "Openhabitat", "Wetland")
periods <- c("S0", "S3")

guild_rasters <- list()

for(guild in guilds) {
  for(period in periods) {
    raster_name <- paste0("richness_prob_", guild, "_", period)
    raster_path <- file.path(guild_raster_dir, paste0(raster_name, ".tif"))
    
    if(file.exists(raster_path)) {
      guild_rasters[[raster_name]] <- rast(raster_path)
      cat(sprintf("    ✓ %s\n", basename(raster_path)))
    } else {
      cat(sprintf("    ✗ : %s\n", basename(raster_path)))
    }
  }
}

cat(sprintf("\n     %d \n", length(guild_rasters)))

# --- 2.3 PSM matched data ---
cat("\n  PSM matched data...\n")
matched_data <- read.csv(psm_file)
cat(sprintf("    %d \n", nrow(matched_data)))
cat(sprintf("    Treatment (GGW): %d\n", sum(matched_data$treatment == 1)))
cat(sprintf("    Control (Non-GGW): %d\n", sum(matched_data$treatment == 0)))

# --- 2.4 ---
cat("\n  ...\n")

if(!file.exists(climate_file)) {
  stop(sprintf(": : %s", climate_file))
}

climate_data <- read.csv(climate_file)
climate_subset <- climate_data %>%
  dplyr::select(any_of(c(id_col, climate_vars))) %>%
  distinct()

matched_data <- matched_data %>%
  left_join(climate_subset, by = id_col)


available_climate <- intersect(climate_vars, names(matched_data))
cat(sprintf("    : %s\n", paste(available_climate, collapse = ", ")))

for(cv in available_climate) {
  na_count <- sum(is.na(matched_data[[cv]]))
  cat(sprintf("      %s: %d  (%d NA)\n", cv, nrow(matched_data) - na_count, na_count))
}

if(length(available_climate) == 0) {
  stop(": ！")
}

# =============================================================================
# 3. guild richness
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  3. guild richness                                     \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

admin_vect <- vect(admin)
guild_by_admin <- data.frame(ID = admin[[id_col]])
names(guild_by_admin)[1] <- id_col

for(raster_name in names(guild_rasters)) {
  cat(sprintf("  : %s\n", raster_name))
  r <- guild_rasters[[raster_name]]
  extracted <- terra::extract(r, admin_vect, fun = mean, na.rm = TRUE)
  guild_by_admin[[raster_name]] <- extracted[, 2]
}

# guildchange
for(guild in guilds) {
  s3_col <- paste0("richness_prob_", guild, "_S3")
  s0_col <- paste0("richness_prob_", guild, "_S0")
  change_col <- paste0(guild, "_change")
  
  if(s3_col %in% names(guild_by_admin) && s0_col %in% names(guild_by_admin)) {
    guild_by_admin[[change_col]] <- guild_by_admin[[s3_col]] - guild_by_admin[[s0_col]]
  }
}

cat(sprintf("\n  : %d \n", nrow(guild_by_admin)))


guild_admin_file <- file.path(output_dir, "Senegal_guild_richness_by_admin.csv")
write.csv(guild_by_admin, guild_admin_file, row.names = FALSE)
cat(sprintf("  ✓ : %s\n", guild_admin_file))

# =============================================================================
# 4. matched data
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  4. guild richnessPSM matched data                             \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

matched_guild <- matched_data %>%
  left_join(guild_by_admin, by = id_col)

cat(sprintf("  : %d\n", nrow(matched_guild)))

for(guild in guilds) {
  change_col <- paste0(guild, "_change")
  if(change_col %in% names(matched_guild)) {
    na_count <- sum(is.na(matched_guild[[change_col]]))
    cat(sprintf("    %s: %d  (%d NA)\n", 
                guild, nrow(matched_guild) - na_count, na_count))
  }
}

# =============================================================================
# 5. PSM
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  5.                         \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

calc_smd <- function(data, var_name) {
  treat <- data %>% filter(treatment == 1)
  ctrl <- data %>% filter(treatment == 0)
  
  mean_t <- mean(treat[[var_name]], na.rm = TRUE)
  mean_c <- mean(ctrl[[var_name]], na.rm = TRUE)
  sd_t <- sd(treat[[var_name]], na.rm = TRUE)
  sd_c <- sd(ctrl[[var_name]], na.rm = TRUE)
  
  pooled_sd <- sqrt((sd_t^2 + sd_c^2) / 2)
  smd <- (mean_t - mean_c) / pooled_sd
  
  return(data.frame(
    Variable = var_name,
    Mean_GGW = mean_t,
    Mean_NonGGW = mean_c,
    SMD = smd,
    Balanced = ifelse(abs(smd) < 0.25, "✓", "✗ ")
  ))
}

cat("  【PSM】\n")
cat(sprintf("  %-15s %12s %12s %8s %12s\n", "Variable", "GGW", "Non-GGW", "SMD", "Status"))
cat("  ─────────────────────────────────────────────────────────────────\n")

balance_climate <- do.call(rbind, lapply(available_climate, function(v) calc_smd(matched_guild, v)))

for(i in 1:nrow(balance_climate)) {
  cat(sprintf("  %-15s %12.5f %12.5f %8.3f %12s\n",
              balance_climate$Variable[i],
              balance_climate$Mean_GGW[i],
              balance_climate$Mean_NonGGW[i],
              balance_climate$SMD[i],
              balance_climate$Balanced[i]))
}

cat("\n  → SMD > 0.25PSM\n")

# =============================================================================
# 6. Guild-specific ATT
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║         6. Guild-specific ATT                            ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

# Guild
guild_labels <- c(
  "Woodland" = "Woodland species",
  "Openhabitat" = "Open-habitat species",
  "Wetland" = "Wetland species"
)


att_results <- data.frame()

for(guild in guilds) {
  
  change_col <- paste0(guild, "_change")
  
  if(!change_col %in% names(matched_guild)) {
    cat(sprintf("  ⚠ %s: \n", guild))
    next
  }
  
  cat(sprintf("【%s】\n", guild_labels[guild]))
  cat("─────────────────────────────────────────────────────────────────────────\n")
  
  treat <- matched_guild %>% filter(treatment == 1)
  ctrl <- matched_guild %>% filter(treatment == 0)
  
  # ===== 1: PSM (Weighted t-test) =====
  wt_result <- wtd.t.test(
    x = treat[[change_col]], y = ctrl[[change_col]],
    weight = treat$weights, weighty = ctrl$weights,
    samedata = FALSE
  )
  
  simple_att <- wt_result$additional["Difference"]
  simple_se <- wt_result$additional["Std. Err"]
  simple_p <- wt_result$coefficients["p.value"]
  
  att_results <- rbind(att_results, data.frame(
    Guild = guild_labels[guild],
    Method = "Simple PSM ()",
    ATT = simple_att,
    SE = simple_se,
    CI_Lower = simple_att - 1.96 * simple_se,
    CI_Upper = simple_att + 1.96 * simple_se,
    p_value = simple_p,
    Sig = add_sig(simple_p)
  ))
  
  # ===== 2: (PSM + WLS) - =====
  formula_dr <- as.formula(paste(change_col, "~ treatment +", 
                                 paste(available_climate, collapse = " + ")))
  
  # NA
  cols_needed <- c("treatment", change_col, available_climate, "weights")
  data_clean <- matched_guild[complete.cases(matched_guild[, cols_needed]), ]
  
  model_dr <- lm(formula_dr, data = data_clean, weights = weights)
  coef_dr <- summary(model_dr)$coefficients
  
  dr_att <- coef_dr["treatment", "Estimate"]
  dr_se <- coef_dr["treatment", "Std. Error"]
  dr_p <- coef_dr["treatment", "Pr(>|t|)"]
  
  att_results <- rbind(att_results, data.frame(
    Guild = guild_labels[guild],
    Method = "Doubly Robust ()",
    ATT = dr_att,
    SE = dr_se,
    CI_Lower = dr_att - 1.96 * dr_se,
    CI_Upper = dr_att + 1.96 * dr_se,
    p_value = dr_p,
    Sig = add_sig(dr_p)
  ))
  

  pct_change <- ifelse(abs(simple_att) > 0.001, 
                       (dr_att - simple_att) / abs(simple_att) * 100, NA)
  

  cat(sprintf("  PSM:     ATT = %7.3f (SE = %.3f), 95%% CI [%.3f, %.3f], p %s %s\n",
              simple_att, simple_se, 
              simple_att - 1.96*simple_se, simple_att + 1.96*simple_se,
              ifelse(simple_p < 0.001, "< 0.001", sprintf("= %.4f", simple_p)),
              add_sig(simple_p)))
  
  cat(sprintf("  :    ATT = %7.3f (SE = %.3f), 95%% CI [%.3f, %.3f], p %s %s\n",
              dr_att, dr_se, 
              dr_att - 1.96*dr_se, dr_att + 1.96*dr_se,
              ifelse(dr_p < 0.001, "< 0.001", sprintf("= %.4f", dr_p)),
              add_sig(dr_p)))
  
  cat(sprintf("  → ATT: %.1f%%\n", pct_change))
  cat("\n")
}

# =============================================================================
# 7. 
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                    ATT                               ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")


dr_results <- att_results %>% filter(Method == "Doubly Robust ()")

cat("【Guild-specific ATT ()】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %10s %10s %18s %6s\n", "Guild", "ATT", "SE", "95% CI", "Sig"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(dr_results)) {
  r <- dr_results[i, ]
  ci_str <- sprintf("[%.3f, %.3f]", r$CI_Lower, r$CI_Upper)
  cat(sprintf("%-25s %10.3f %10.3f %18s %6s\n",
              r$Guild, r$ATT, r$SE, ci_str, r$Sig))
}
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n")

# =============================================================================
# 8. Total Richness
# =============================================================================

cat("\n")
cat("【Total Richness ATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

if("prob_change" %in% names(matched_guild)) {
  
  treat <- matched_guild %>% filter(treatment == 1)
  ctrl <- matched_guild %>% filter(treatment == 0)
  
  # Total richness ATT
  formula_total <- as.formula(paste("prob_change ~ treatment +", 
                                    paste(available_climate, collapse = " + ")))
  cols_needed <- c("treatment", "prob_change", available_climate, "weights")
  data_clean <- matched_guild[complete.cases(matched_guild[, cols_needed]), ]
  
  model_total <- lm(formula_total, data = data_clean, weights = weights)
  coef_total <- summary(model_total)$coefficients
  
  total_att <- coef_total["treatment", "Estimate"]
  total_se <- coef_total["treatment", "Std. Error"]
  total_p <- coef_total["treatment", "Pr(>|t|)"]
  
  cat(sprintf("\n  Total Richness (): ATT = %.3f (SE = %.3f) %s\n", 
              total_att, total_se, add_sig(total_p)))
  
  cat("\n  Guild:\n")
  guild_sum <- 0
  for(i in 1:nrow(dr_results)) {
    r <- dr_results[i, ]
    pct <- ifelse(abs(total_att) > 0.001, r$ATT / total_att * 100, NA)
    guild_sum <- guild_sum + r$ATT
    cat(sprintf("    %-22s: %7.3f (%6.1f%% of total) %s\n",
                r$Guild, r$ATT, pct, r$Sig))
  }
  
  cat(sprintf("\n  Guild ATT: %.3f (Total ATT: %.3f)\n", guild_sum, total_att))
  

  cat("\n  【】\n")
  
  sig_positive <- dr_results %>% filter(p_value < 0.05 & ATT > 0)
  sig_negative <- dr_results %>% filter(p_value < 0.05 & ATT < 0)
  not_sig <- dr_results %>% filter(p_value >= 0.05)
  
  if(nrow(sig_positive) > 0) {
    cat(sprintf("  ✓ : %s\n", paste(sig_positive$Guild, collapse = ", ")))
  }
  if(nrow(sig_negative) > 0) {
    cat(sprintf("  ✗ : %s\n", paste(sig_negative$Guild, collapse = ", ")))
  }
  if(nrow(not_sig) > 0) {
    cat(sprintf("  ○ : %s\n", paste(not_sig$Guild, collapse = ", ")))
  }
}

# =============================================================================
# 9. 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  9.                                                          \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ATT
att_results$Country <- COUNTRY
att_file <- file.path(output_dir, "Senegal_guild_DR_ATT_results.csv")
write.csv(att_results, att_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", att_file))


balance_climate$Country <- COUNTRY
balance_file <- file.path(output_dir, "Senegal_guild_climate_balance.csv")
write.csv(balance_climate, balance_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", balance_file))


matched_guild_file <- file.path(output_dir, "Senegal_matched_with_guild_richness.csv")
write.csv(matched_guild, matched_guild_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", matched_guild_file))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Senegal Guild！                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")
