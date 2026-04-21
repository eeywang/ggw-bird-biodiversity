# =============================================================================
# Step 3: - ATT (Δbio)

# : 
#   - Step1{COUNTRY}_1to1_matched_data_50threshold.csv
#   - 01_calculate_climate_trends.R


#   1. PSM+ (Δbio)
#   2. ATT (Weighted t-test, WLS, Mixed Effects)
#   3. ATT ()
#   4. ATTATT

# :
#   :   Δbio2, Δbio4, Δbio13, Δbio18
#   : Δbio3, Δbio4, Δbio7, Δbio9, Δbio13, Δbio19


#   - {COUNTRY}_ATT_3methods_50threshold.csv     (ATT)
#   - {COUNTRY}_DR_richness_ATT_50threshold.csv   (ATT)
#   - {COUNTRY}_50threshold_Forest_ATT.png/pdf    ()
#   - {COUNTRY}_50threshold_Violin_Change.png/pdf ()
#   - {COUNTRY}_50threshold_Boxplot_Temporal.png/pdf ()
# =============================================================================

library(dplyr)
library(tidyr)
library(weights)
library(survey)
library(lme4)
library(lmerTest)
library(ggplot2)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Step 3:  - ATT (Δbio)         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================

# =============================================================================

COUNTRY <- "Ethiopia"  # ← 

if(COUNTRY == "Senegal") {
  
  matched_file <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold/Senegal_1to1_matched_data_50threshold.csv"
  climate_file <- "E:/2026.1.8_biomod2/senegal_results/climate_trends/Senegal_PSM_with_climate_trends.csv"
  climate_vars <- c("bio2_trend", "bio4_trend", "bio13_trend", "bio18_trend")
  id_col       <- "GID_4"
  output_dir   <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold"
  
} else if(COUNTRY == "Ethiopia") {
  
  matched_file <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold/Ethiopia_1to1_matched_data_50threshold.csv"
  climate_file <- "E:/2026.1.8_biomod2/ethiopia_results/climate_trends/Ethiopia_PSM_with_climate_trends.csv"
  climate_vars <- c("bio3_trend", "bio4_trend", "bio7_trend", "bio9_trend", "bio13_trend", "bio19_trend")
  id_col       <- "GID_3"
  output_dir   <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold"

} else if(COUNTRY == "Nigeria") {
  
  matched_file <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold/Nigeria_1to1_matched_data_50threshold.csv"
  climate_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"
  climate_vars <- c("bio1_trend", "bio12_trend")
  id_col       <- "GID_2"
  output_dir   <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold"
}

cat(sprintf(": %s\n", COUNTRY))
cat(sprintf(": %s\n", paste(climate_vars, collapse = ", ")))

# =============================================================================
# 3.1 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  3.1 PSM +                                  \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# PSM
if(!file.exists(matched_file)) {
  stop(sprintf(": PSM: %s\n Step1", matched_file))
}
matched_data <- read.csv(matched_file)
cat(sprintf("  PSM: %d \n", nrow(matched_data)))


available_climate <- c()

if(file.exists(climate_file)) {
  climate_data <- read.csv(climate_file)
  climate_cols <- c(id_col, climate_vars)
  climate_subset <- climate_data[, intersect(climate_cols, names(climate_data))]
  climate_subset <- climate_subset[!duplicated(climate_subset[[id_col]]), ]
  
  matched_data <- merge(matched_data, climate_subset, by = id_col, all.x = TRUE)
  available_climate <- intersect(climate_vars, names(matched_data))
  
  cat(sprintf("  ✓ \n"))
  cat(sprintf("  : %s\n", paste(available_climate, collapse = ", ")))
  
  # NA
  for(cv in available_climate) {
    na_count <- sum(is.na(matched_data[[cv]]))
    cat(sprintf("    %s: %d  (%d NA)\n", cv, nrow(matched_data) - na_count, na_count))
  }
} else {
  cat("  ⚠ : ATT\n")
  cat(sprintf("  : %s\n", climate_file))
}

# =============================================================================
# 3.2 
# =============================================================================

add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

# =============================================================================
# 3.3 ATT - 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  3.3  ATT -                                        \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

treated <- matched_data[matched_data$treatment == 1, ]
control <- matched_data[matched_data$treatment == 0, ]

att_results_3methods <- data.frame()

richness_metrics <- c("bin_change", "prob_change", "tss_change")
metric_labels <- c("Binary-based", "Probability-based", "TSS-weighted")

# --- 1: Weighted t-test ---
cat("--- 1: Weighted t-test ---\n")

for(i in 1:length(richness_metrics)) {
  metric <- richness_metrics[i]
  label <- metric_labels[i]
  
  wt_result <- wtd.t.test(
    x = treated[[metric]], y = control[[metric]],
    weight = treated$weights, weighty = control$weights,
    samedata = FALSE
  )
  
  att <- wt_result$additional["Difference"]
  se <- wt_result$additional["Std. Err"]
  p_val <- wt_result$coefficients["p.value"]
  
  att_results_3methods <- rbind(att_results_3methods, data.frame(
    Richness_Metric = label, Method = "Weighted t-test",
    ATT = att, SE = se,
    CI_Lower = att - 1.96 * se, CI_Upper = att + 1.96 * se,
    p_value = p_val
  ))
  cat(sprintf("  %s: ATT = %.3f [%.3f, %.3f] p = %s\n",
              label, att, att - 1.96*se, att + 1.96*se, format(p_val, digits = 3)))
}

# --- 2: WLS + Covariates ---
cat("\n--- 2: WLS + Covariates ---\n")

design_matched <- svydesign(ids = ~1, weights = ~weights, data = matched_data)

for(i in 1:length(richness_metrics)) {
  metric <- richness_metrics[i]
  label <- metric_labels[i]
  
  wls_formula <- as.formula(paste(metric, "~ treatment + prec_mean + temp_mean + ndvi_mean + slope_mean"))
  wls_model <- svyglm(wls_formula, design = design_matched)
  
  coef_sum <- summary(wls_model)$coefficients
  att <- coef_sum["treatment", "Estimate"]
  se <- coef_sum["treatment", "Std. Error"]
  p_val <- coef_sum["treatment", "Pr(>|t|)"]
  
  att_results_3methods <- rbind(att_results_3methods, data.frame(
    Richness_Metric = label, Method = "WLS + Covariates",
    ATT = att, SE = se,
    CI_Lower = att - 1.96 * se, CI_Upper = att + 1.96 * se,
    p_value = p_val
  ))
  cat(sprintf("  %s: ATT = %.3f [%.3f, %.3f] p = %s\n",
              label, att, att - 1.96*se, att + 1.96*se, format(p_val, digits = 3)))
}

# --- 3: Mixed Effects Model ---
cat("\n--- 3: Mixed Effects Model ---\n")

for(i in 1:length(richness_metrics)) {
  metric <- richness_metrics[i]
  label <- metric_labels[i]
  
  lmer_formula <- as.formula(paste(metric, "~ treatment + (1|subclass)"))
  lmer_model <- lmer(lmer_formula, data = matched_data, weights = weights)
  
  coef_sum <- summary(lmer_model)$coefficients
  att <- coef_sum["treatment", "Estimate"]
  se <- coef_sum["treatment", "Std. Error"]
  p_val <- coef_sum["treatment", "Pr(>|t|)"]
  
  att_results_3methods <- rbind(att_results_3methods, data.frame(
    Richness_Metric = label, Method = "Mixed Effects",
    ATT = att, SE = se,
    CI_Lower = att - 1.96 * se, CI_Upper = att + 1.96 * se,
    p_value = p_val
  ))
  cat(sprintf("  %s: ATT = %.3f [%.3f, %.3f] p = %s\n",
              label, att, att - 1.96*se, att + 1.96*se, format(p_val, digits = 3)))
}

att_results_3methods$Significance <- add_sig(att_results_3methods$p_value)
rownames(att_results_3methods) <- NULL

# =============================================================================
# 3.4 ATT - Δbio
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  3.4 ATT - Δbio                                 \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")


calc_dr_att <- function(data, var_name, var_label, climate_vars) {
  treat <- data[data$treatment == 1, ]
  ctrl <- data[data$treatment == 0, ]
  simple_att <- mean(treat[[var_name]], na.rm = TRUE) - mean(ctrl[[var_name]], na.rm = TRUE)
  simple_t <- t.test(treat[[var_name]], ctrl[[var_name]])
  simple_p <- simple_t$p.value
  
  available <- intersect(climate_vars, names(data))
  
  if(length(available) == 0) {
    return(data.frame(Variable = var_label, Simple_ATT = simple_att, Simple_p = simple_p,
                      DR_ATT = NA, DR_SE = NA, DR_p = NA, Change_pct = NA))
  }
  
  cols_needed <- c("treatment", var_name, available)
  data_clean <- na.omit(data[, cols_needed])
  
  if(nrow(data_clean) < 20) {
    return(data.frame(Variable = var_label, Simple_ATT = simple_att, Simple_p = simple_p,
                      DR_ATT = NA, DR_SE = NA, DR_p = NA, Change_pct = NA))
  }
  
  formula_adj <- as.formula(paste(var_name, "~ treatment +", paste(available, collapse = " + ")))
  model_adj <- lm(formula_adj, data = data_clean)
  coef_adj <- summary(model_adj)$coefficients
  
  dr_att <- coef_adj["treatment", "Estimate"]
  dr_se <- coef_adj["treatment", "Std. Error"]
  dr_p <- coef_adj["treatment", "Pr(>|t|)"]
  change_pct <- ifelse(abs(simple_att) > 0.001, (dr_att - simple_att) / abs(simple_att) * 100, NA)
  
  data.frame(Variable = var_label, Simple_ATT = simple_att, Simple_p = simple_p,
             DR_ATT = dr_att, DR_SE = dr_se, DR_p = dr_p, Change_pct = change_pct)
}

# ATT
dr_richness <- data.frame()
for(i in 1:length(richness_metrics)) {
  result <- calc_dr_att(matched_data, richness_metrics[i], metric_labels[i], climate_vars)
  dr_richness <- rbind(dr_richness, result)
}

dr_richness$Simple_Sig <- add_sig(dr_richness$Simple_p)
dr_richness$DR_Sig <- ifelse(is.na(dr_richness$DR_p), "", add_sig(dr_richness$DR_p))

# =============================================================================
# 3.5 
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  %s ATT (50)                                ║\n", COUNTRY))
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("【 - ATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-22s %12s %15s %10s\n", "", "ATT", "ATT", "%"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(dr_richness)) {
  r <- dr_richness[i, ]
  dr_str <- ifelse(is.na(r$DR_ATT), "NA", sprintf("%.2f%s", r$DR_ATT, r$DR_Sig))
  pct_str <- ifelse(is.na(r$Change_pct), "NA", sprintf("%.1f%%", r$Change_pct))
  cat(sprintf("%-22s %10.2f%s %15s %10s\n",
              r$Variable, r$Simple_ATT, r$Simple_Sig, dr_str, pct_str))
}
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n")


cat("\n")
bin_dr <- dr_richness[dr_richness$Variable == "Binary-based", ]
if(nrow(bin_dr) > 0) {
  cat("【Binary - 】\n")
  cat(sprintf("  ATT: %.2f%s\n", bin_dr$Simple_ATT, bin_dr$Simple_Sig))
  if(!is.na(bin_dr$DR_ATT)) {
    cat(sprintf("  ATT: %.2f%s\n", bin_dr$DR_ATT, bin_dr$DR_Sig))
    cat(sprintf("  : %.1f%%\n", bin_dr$Change_pct))
    if(!is.na(bin_dr$DR_p) && bin_dr$DR_p < 0.05) {
      if(bin_dr$DR_ATT > 0) { cat("  → ✓ GGW\n")
      } else { cat("  → ⚠ GGW\n") }
    } else {
      cat("  → GGW\n")
    }
  }
}

# =============================================================================
# 3.6 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  3.6                                                         \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ATT
att_results_3methods$Country <- COUNTRY
att3_file <- file.path(output_dir, paste0(COUNTRY, "_ATT_3methods_50threshold.csv"))
write.csv(att_results_3methods, att3_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", att3_file))

# ATT
dr_richness$Country <- COUNTRY
dr_file <- file.path(output_dir, paste0(COUNTRY, "_DR_richness_ATT_50threshold.csv"))
write.csv(dr_richness, dr_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", dr_file))

# =============================================================================
# 3.7 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  3.7                                                           \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

color_ggw <- "#3B9AB2"
color_nonggw <- "#E78AC3"
method_colors <- c("WLS + Covariates" = "#E67E22", "Weighted t-test" = "#27AE60", "Mixed Effects" = "#8E7CC3")

# ---- 1: ----
cat("  ...\n")

att_results_3methods$Method <- factor(att_results_3methods$Method,
                                      levels = c("Mixed Effects", "Weighted t-test", "WLS + Covariates"))
att_results_3methods$Richness_Metric <- factor(att_results_3methods$Richness_Metric,
                                               levels = c("Binary-based", "Probability-based", "TSS-weighted"))

fig_forest <- ggplot(att_results_3methods, aes(x = ATT, y = Method, color = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, linewidth = 1) +
  geom_point(size = 3.5) +
  geom_text(aes(x = CI_Upper + max(abs(att_results_3methods$CI_Upper)) * 0.08, label = Significance),
            hjust = 0, size = 4, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Richness_Metric, scales = "free_x", nrow = 1) +
  scale_color_manual(values = method_colors) +
  labs(title = sprintf("%s: Treatment Effect on Bird Richness (50-threshold, 1:1 Matching)", COUNTRY),
       subtitle = sprintf("Matched: %d Treatment vs %d Control",
                          sum(matched_data$treatment == 1), sum(matched_data$treatment == 0)),
       x = "Average Treatment Effect (ATT)", y = NULL,
       caption = "*** p<0.001, ** p<0.01, * p<0.05") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10),
        plot.caption = element_text(hjust = 0.5, color = "gray50"),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(output_dir, paste0(COUNTRY, "_50threshold_Forest_ATT.png")), fig_forest,
       width = 14, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, paste0(COUNTRY, "_50threshold_Forest_ATT.pdf")), fig_forest,
       width = 14, height = 5, bg = "white")
cat("    ✓ Forest ATT\n")

# ---- 2: ----
cat("  ...\n")

violin_data <- matched_data %>%
  mutate(Group = ifelse(treatment == 1, "GGW", "Non-GGW")) %>%
  dplyr::select(Group, bin_change, prob_change, tss_change, weights) %>%
  pivot_longer(cols = c(bin_change, prob_change, tss_change),
               names_to = "Metric", values_to = "Change") %>%
  mutate(Group = factor(Group, levels = c("GGW", "Non-GGW")),
         Metric = factor(Metric,
                         levels = c("bin_change", "prob_change", "tss_change"),
                         labels = c("Binary-based", "Probability-based", "TSS-weighted")))

mean_data <- violin_data %>%
  group_by(Group, Metric) %>%
  summarise(mean_change = weighted.mean(Change, weights, na.rm = TRUE), .groups = "drop")

fig_violin <- ggplot(violin_data, aes(x = Group, y = Change, fill = Group)) +
  geom_violin(alpha = 0.6, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = 0.8) +
  geom_point(data = mean_data, aes(y = mean_change),
             shape = 23, fill = "black", color = "white", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("GGW" = color_ggw, "Non-GGW" = color_nonggw)) +
  labs(title = sprintf("%s: Richness Change Distribution (50-threshold, 1:1 Matched)", COUNTRY),
       subtitle = sprintf("N = %d pairs | Diamonds = weighted means", nrow(matched_data)/2),
       x = NULL, y = "Change in Species Richness") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "none")

ggsave(file.path(output_dir, paste0(COUNTRY, "_50threshold_Violin_Change.png")), fig_violin,
       width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, paste0(COUNTRY, "_50threshold_Violin_Change.pdf")), fig_violin,
       width = 12, height = 5, bg = "white")
cat("    ✓ Violin\n")

# ---- 3: ----
cat("  ...\n")

boxplot_data <- matched_data %>%
  mutate(Group = ifelse(treatment == 1, "GGW", "Non-GGW")) %>%
  dplyr::select(Group, weights, bin_before, bin_after, prob_before, prob_after, tss_before, tss_after) %>%
  pivot_longer(cols = -c(Group, weights),
               names_to = c("Metric", "Period"),
               names_pattern = "(bin|prob|tss)_(before|after)",
               values_to = "Richness") %>%
  mutate(Group = factor(Group, levels = c("GGW", "Non-GGW")),
         Period = factor(Period, levels = c("before", "after"), labels = c("2007-2015", "2016-2024")),
         Metric = factor(Metric, levels = c("bin", "prob", "tss"),
                         labels = c("Binary-based", "Probability-based", "TSS-weighted")))

summary_box <- boxplot_data %>%
  group_by(Group, Metric, Period) %>%
  summarise(mean = weighted.mean(Richness, weights, na.rm = TRUE), .groups = "drop")

fig_boxplot <- ggplot(boxplot_data, aes(x = Period, y = Richness, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, position = position_dodge(width = 0.75), width = 0.6) +
  geom_line(data = summary_box, aes(x = Period, y = mean, group = Group, color = Group),
            position = position_dodge(width = 0.75), linewidth = 1) +
  geom_point(data = summary_box, aes(x = Period, y = mean, color = Group),
             position = position_dodge(width = 0.75), size = 2.5, shape = 21, fill = "white") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("GGW" = color_ggw, "Non-GGW" = color_nonggw)) +
  scale_color_manual(values = c("GGW" = color_ggw, "Non-GGW" = color_nonggw)) +
  labs(title = sprintf("%s: Temporal Change (50-threshold, 1:1 Matched)", COUNTRY),
       subtitle = "Boxes show distribution; lines connect weighted means",
       x = "Period", y = "Species Richness") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(output_dir, paste0(COUNTRY, "_50threshold_Boxplot_Temporal.png")), fig_boxplot,
       width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, paste0(COUNTRY, "_50threshold_Boxplot_Temporal.pdf")), fig_boxplot,
       width = 12, height = 5, bg = "white")
cat("    ✓ Temporal\n")

# =============================================================================

# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Step 3 ！                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  :\n"))
cat(sprintf("    - %s\n", att3_file))
cat(sprintf("    - %s\n", dr_file))
cat(sprintf("\n:  Step4_scenario_DR_ATT_50threshold.R\n"))
