# =============================================================================
# Step 4: - ATT (Δbio)

# : 
#   - Step2{COUNTRY}_matched_with_scenarios_50threshold.csv
#   - 01_calculate_climate_trends.R


#   1. PSM+ (Δbio)
#   2. ATT
#   3. ATT ()
#   4. 

# :
#   :   Δbio2, Δbio4, Δbio13, Δbio18
#   : Δbio3, Δbio4, Δbio7, Δbio9, Δbio13, Δbio19


#   - {COUNTRY}_DR_scenario_ATT_50threshold.csv       (ATT)
#   - {COUNTRY}_scenario_boxplot_50threshold.png       ()
#   - {COUNTRY}_effects_barplot_50threshold.png        ()
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Step 4:  - ATT (Δbio)       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================

# =============================================================================

COUNTRY <- "Ethiopia"  # ← 

if(COUNTRY == "Senegal") {
  
  scenario_matched_file <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold/Senegal_matched_with_scenarios_50threshold.csv"
  climate_file <- "E:/2026.1.8_biomod2/senegal_results/climate_trends/Senegal_PSM_with_climate_trends.csv"
  climate_vars <- c("bio2_trend", "bio4_trend", "bio13_trend", "bio18_trend")
  id_col       <- "GID_4"
  output_dir   <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_50threshold"
  
} else if(COUNTRY == "Ethiopia") {
  
  scenario_matched_file <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold/Ethiopia_matched_with_scenarios_50threshold.csv"
  climate_file <- "E:/2026.1.8_biomod2/ethiopia_results/climate_trends/Ethiopia_PSM_with_climate_trends.csv"
  climate_vars <- c("bio3_trend", "bio4_trend", "bio7_trend", "bio9_trend", "bio13_trend", "bio19_trend")
  id_col       <- "GID_3"
  output_dir   <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_50threshold"

} else if(COUNTRY == "Nigeria") {
  
  matched_file <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold/Nigeria_matched_with_scenarios_50threshold.csv"
  climate_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"
  climate_vars <- c("bio1_trend", "bio12_trend")
  id_col       <- "GID_2"
  output_dir   <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_50threshold"
}

cat(sprintf(": %s\n", COUNTRY))
cat(sprintf(": %s\n", paste(climate_vars, collapse = ", ")))

# =============================================================================
# 4.1 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  4.1 () +                         \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

if(!file.exists(scenario_matched_file)) {
  stop(sprintf(": : %s\n Step2", scenario_matched_file))
}

matched_data <- read.csv(scenario_matched_file)
cat(sprintf("  (): %d \n", nrow(matched_data)))
cat(sprintf("    Treatment: %d, Control: %d\n",
            sum(matched_data$treatment == 1), sum(matched_data$treatment == 0)))


available_climate <- c()

if(file.exists(climate_file)) {
  climate_data <- read.csv(climate_file)
  climate_cols <- c(id_col, climate_vars)
  climate_subset <- climate_data[, intersect(climate_cols, names(climate_data))]
  climate_subset <- climate_subset[!duplicated(climate_subset[[id_col]]), ]
  
  # - matched_data
  existing_climate <- intersect(climate_vars, names(matched_data))
  if(length(existing_climate) > 0) {
    matched_data <- matched_data[, !(names(matched_data) %in% existing_climate)]
  }
  
  matched_data <- merge(matched_data, climate_subset, by = id_col, all.x = TRUE)
  available_climate <- intersect(climate_vars, names(matched_data))
  
  cat(sprintf("  ✓ \n"))
  cat(sprintf("  : %s\n", paste(available_climate, collapse = ", ")))
  
  for(cv in available_climate) {
    na_count <- sum(is.na(matched_data[[cv]]))
    cat(sprintf("    %s: %d  (%d NA)\n", cv, nrow(matched_data) - na_count, na_count))
  }
} else {
  cat("  ⚠ : ATT\n")
}

# =============================================================================
# 4.2 
# =============================================================================

add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

# ATT
calc_att <- function(data, var_name, var_label) {
  treat <- data[data$treatment == 1, ]
  ctrl  <- data[data$treatment == 0, ]
  
  treat_mean <- mean(treat[[var_name]], na.rm = TRUE)
  ctrl_mean  <- mean(ctrl[[var_name]], na.rm = TRUE)
  t_result   <- t.test(treat[[var_name]], ctrl[[var_name]])
  
  data.frame(
    Variable = var_label,
    Treatment_Mean = treat_mean,
    Control_Mean = ctrl_mean,
    ATT = treat_mean - ctrl_mean,
    SE = sqrt(var(treat[[var_name]], na.rm = TRUE)/nrow(treat) + 
              var(ctrl[[var_name]], na.rm = TRUE)/nrow(ctrl)),
    p_value = t_result$p.value
  )
}

# ATT
calc_dr_att <- function(data, var_name, var_label, climate_vars) {
  treat <- data[data$treatment == 1, ]
  ctrl  <- data[data$treatment == 0, ]
  simple_att <- mean(treat[[var_name]], na.rm = TRUE) - mean(ctrl[[var_name]], na.rm = TRUE)
  simple_t   <- t.test(treat[[var_name]], ctrl[[var_name]])
  simple_p   <- simple_t$p.value
  
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
  dr_se  <- coef_adj["treatment", "Std. Error"]
  dr_p   <- coef_adj["treatment", "Pr(>|t|)"]
  change_pct <- ifelse(abs(simple_att) > 0.001, (dr_att - simple_att) / abs(simple_att) * 100, NA)
  
  data.frame(Variable = var_label, Simple_ATT = simple_att, Simple_p = simple_p,
             DR_ATT = dr_att, DR_SE = dr_se, DR_p = dr_p, Change_pct = change_pct)
}

# =============================================================================
# 4.3 - ATT
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  4.3  - ATT                                          \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# --- Binary---
cat("【Binary】\n")
binary_effects_att <- rbind(
  calc_att(matched_data, "binary_total", "Total Change"),
  calc_att(matched_data, "binary_climate", "Climate Effect"),
  calc_att(matched_data, "binary_vegetation", "Vegetation Effect"),
  calc_att(matched_data, "binary_interaction", "Interaction Effect")
)
binary_effects_att$Significance <- add_sig(binary_effects_att$p_value)
binary_effects_att$Category <- "Binary Effects"

cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %10s %10s %10s %8s\n", "", "GGW", "Non-GGW", "ATT", "Sig"))
cat("─────────────────────────────────────────────────────────────────────────\n")
for(i in 1:nrow(binary_effects_att)) {
  cat(sprintf("%-25s %10.2f %10.2f %10.2f %8s\n",
              binary_effects_att$Variable[i], binary_effects_att$Treatment_Mean[i],
              binary_effects_att$Control_Mean[i], binary_effects_att$ATT[i],
              binary_effects_att$Significance[i]))
}

# --- Probability---
cat("\n【Probability】\n")
prob_effects_att <- rbind(
  calc_att(matched_data, "prob_total", "Total Change"),
  calc_att(matched_data, "prob_climate", "Climate Effect"),
  calc_att(matched_data, "prob_vegetation", "Vegetation Effect"),
  calc_att(matched_data, "prob_interaction", "Interaction Effect")
)
prob_effects_att$Significance <- add_sig(prob_effects_att$p_value)
prob_effects_att$Category <- "Probability Effects"

cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %10s %10s %10s %8s\n", "", "GGW", "Non-GGW", "ATT", "Sig"))
cat("─────────────────────────────────────────────────────────────────────────\n")
for(i in 1:nrow(prob_effects_att)) {
  cat(sprintf("%-25s %10.2f %10.2f %10.2f %8s\n",
              prob_effects_att$Variable[i], prob_effects_att$Treatment_Mean[i],
              prob_effects_att$Control_Mean[i], prob_effects_att$ATT[i],
              prob_effects_att$Significance[i]))
}

# --- TSS---
cat("\n【TSS-weighted】\n")
tss_effects_att <- rbind(
  calc_att(matched_data, "tss_total", "Total Change"),
  calc_att(matched_data, "tss_climate", "Climate Effect"),
  calc_att(matched_data, "tss_vegetation", "Vegetation Effect"),
  calc_att(matched_data, "tss_interaction", "Interaction Effect")
)
tss_effects_att$Significance <- add_sig(tss_effects_att$p_value)
tss_effects_att$Category <- "TSS Effects"

cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %10s %10s %10s %8s\n", "", "GGW", "Non-GGW", "ATT", "Sig"))
cat("─────────────────────────────────────────────────────────────────────────\n")
for(i in 1:nrow(tss_effects_att)) {
  cat(sprintf("%-25s %10.2f %10.2f %10.2f %8s\n",
              tss_effects_att$Variable[i], tss_effects_att$Treatment_Mean[i],
              tss_effects_att$Control_Mean[i], tss_effects_att$ATT[i],
              tss_effects_att$Significance[i]))
}

# ATT
all_simple_att <- rbind(binary_effects_att, prob_effects_att, tss_effects_att)

# =============================================================================
# 4.4 ATT - 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  4.4 ATT -  (Δbio)                  \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")


dr_vars <- list(
  # Binary
  list(var = "binary_climate",     label = "Climate Effect (Binary)",     cat = "Binary"),
  list(var = "binary_vegetation",  label = "Vegetation Effect (Binary)",  cat = "Binary"),
  list(var = "binary_interaction", label = "Interaction Effect (Binary)", cat = "Binary"),
  list(var = "binary_total",       label = "Total Change (Binary)",       cat = "Binary"),
  # Probability
  list(var = "prob_climate",     label = "Climate Effect (Prob)",     cat = "Probability"),
  list(var = "prob_vegetation",  label = "Vegetation Effect (Prob)",  cat = "Probability"),
  list(var = "prob_interaction", label = "Interaction Effect (Prob)", cat = "Probability"),
  list(var = "prob_total",       label = "Total Change (Prob)",       cat = "Probability"),
  # TSS
  list(var = "tss_climate",     label = "Climate Effect (TSS)",     cat = "TSS"),
  list(var = "tss_vegetation",  label = "Vegetation Effect (TSS)",  cat = "TSS"),
  list(var = "tss_interaction", label = "Interaction Effect (TSS)", cat = "TSS"),
  list(var = "tss_total",       label = "Total Change (TSS)",       cat = "TSS")
)

dr_results <- data.frame()
for(v in dr_vars) {
  if(v$var %in% names(matched_data)) {
    result <- calc_dr_att(matched_data, v$var, v$label, climate_vars)
    result$Category <- v$cat
    dr_results <- rbind(dr_results, result)
  }
}

dr_results$Simple_Sig <- add_sig(dr_results$Simple_p)
dr_results$DR_Sig <- ifelse(is.na(dr_results$DR_p), "", add_sig(dr_results$DR_p))

# =============================================================================
# 4.5 
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  %s  - ATT (50)                       ║\n", COUNTRY))
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

for(cat_name in c("Binary", "Probability", "TSS")) {
  cat_dr <- dr_results[dr_results$Category == cat_name, ]
  
  if(nrow(cat_dr) > 0) {
    cat(sprintf("【%s -  ()】\n", cat_name))
    cat("─────────────────────────────────────────────────────────────────────────\n")
    cat(sprintf("%-30s %12s %15s %10s %15s\n", "", "ATT", "ATT", "%", ""))
    cat("─────────────────────────────────────────────────────────────────────────\n")
    
    for(i in 1:nrow(cat_dr)) {
      r <- cat_dr[i, ]
      dr_str <- ifelse(is.na(r$DR_ATT), "NA", sprintf("%.2f%s", r$DR_ATT, r$DR_Sig))
      pct_str <- ifelse(is.na(r$Change_pct), "NA", sprintf("%.1f%%", r$Change_pct))
      

      if(grepl("Climate", r$Variable)) {
        if(is.na(r$DR_p)) { interp <- ""
        } else if(r$Simple_p < 0.05 && r$DR_p > 0.05) { interp <- "✓ "
        } else if(r$Simple_p > 0.05 && r$DR_p > 0.05) { interp <- "○ "
        } else if(!is.na(r$Change_pct) && abs(r$Change_pct) < 15) { interp <- ""
        } else { interp <- "" }
      } else if(grepl("Vegetation", r$Variable)) {
        if(is.na(r$DR_p)) { interp <- ""
        } else {
          dir_str <- ifelse(r$DR_ATT > 0, "+", "-")
          interp <- ifelse(r$DR_p < 0.05, sprintf("GGW(%s)", dir_str), "")
        }
      } else if(grepl("Interaction", r$Variable)) {
        if(is.na(r$DR_p)) { interp <- ""
        } else { interp <- ifelse(r$DR_p < 0.05, "", "") }
      } else {
        # Total Change
        if(is.na(r$DR_p)) { interp <- ""
        } else { interp <- ifelse(r$DR_p < 0.05, "", "") }
      }
      
      cat(sprintf("%-30s %10.2f%s %15s %10s %15s\n",
                  r$Variable, r$Simple_ATT, r$Simple_Sig, dr_str, pct_str, interp))
    }
    cat("─────────────────────────────────────────────────────────────────────────\n\n")
  }
}
cat("* p<0.05, ** p<0.01, *** p<0.001\n")

# =============================================================================
# 4.6 
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  %s                                                     ║\n", COUNTRY))
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# Binary
bin_climate <- dr_results[dr_results$Variable == "Climate Effect (Binary)", ]
if(nrow(bin_climate) > 0 && !is.na(bin_climate$DR_p)) {
  cat("  【 (Binary)】\n")
  cat(sprintf("    : %.2f%s\n", bin_climate$Simple_ATT, bin_climate$Simple_Sig))
  cat(sprintf("    : %.2f%s\n", bin_climate$DR_ATT, bin_climate$DR_Sig))
  cat(sprintf("    : %.1f%%\n", bin_climate$Change_pct))
  
  if(bin_climate$Simple_p < 0.05 && bin_climate$DR_p > 0.05) {
    cat("    → ✓ ！\n")
  } else if(bin_climate$Simple_p > 0.05) {
    cat("    → ○ \n")
  } else if(abs(bin_climate$Change_pct) < 15) {
    cat("    → \n")
  } else {
    cat("    → \n")
  }
}

# Binary
bin_veg <- dr_results[dr_results$Variable == "Vegetation Effect (Binary)", ]
if(nrow(bin_veg) > 0 && !is.na(bin_veg$DR_p)) {
  cat("\n  【 (GGW) (Binary)】\n")
  cat(sprintf("    : %.2f%s\n", bin_veg$DR_ATT, bin_veg$DR_Sig))
  
  if(bin_veg$DR_p < 0.05) {
    if(bin_veg$DR_ATT > 0) { cat("    → ✓ GGW\n")
    } else { cat("    → ⚠ GGW\n") }
  } else { cat("    → GGW\n") }
}

# =============================================================================
# 4.7 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  4.7                                                         \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ATT
all_simple_att$Country <- COUNTRY
simple_file <- file.path(output_dir, paste0(COUNTRY, "_ATT_scenario_simple_50threshold.csv"))
write.csv(all_simple_att, simple_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", simple_file))

# ATT
dr_results$Country <- COUNTRY
dr_file <- file.path(output_dir, paste0(COUNTRY, "_DR_scenario_ATT_50threshold.csv"))
write.csv(dr_results, dr_file, row.names = FALSE)
cat(sprintf("  ✓ %s\n", dr_file))

# =============================================================================
# 4.8 
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  4.8                                                           \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

color_ggw <- "#3B9AB2"
color_nonggw <- "#E78AC3"

# ---- 1: ----
cat("  ...\n")

scenario_plot_data <- matched_data %>%
  mutate(Group = ifelse(treatment == 1, "GGW", "Non-GGW")) %>%
  dplyr::select(Group, binary_climate, binary_vegetation, binary_interaction, binary_total) %>%
  pivot_longer(cols = -Group, names_to = "Effect", values_to = "Value") %>%
  mutate(Effect = case_when(
    Effect == "binary_total" ~ "Total Change",
    Effect == "binary_climate" ~ "Climate Effect",
    Effect == "binary_vegetation" ~ "Vegetation Effect",
    Effect == "binary_interaction" ~ "Interaction Effect"),
    Effect = factor(Effect, levels = c("Climate Effect", "Vegetation Effect", "Interaction Effect", "Total Change")))

fig_scenario <- ggplot(scenario_plot_data, aes(x = Effect, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("GGW" = color_ggw, "Non-GGW" = color_nonggw)) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x.npc = 0.5) +
  labs(title = sprintf("%s: Scenario Decomposition (50-threshold)", COUNTRY),
       subtitle = "Binary Richness - PSM 1:1 Matched",
       x = NULL, y = "Effect on Species Richness", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(output_dir, paste0(COUNTRY, "_scenario_boxplot_50threshold.png")),
       fig_scenario, width = 10, height = 6, dpi = 300, bg = "white")
cat("    ✓ \n")

# ---- 2: ----
cat("  ...\n")

summary_bar <- binary_effects_att %>%
  dplyr::select(Variable, Treatment_Mean, Control_Mean) %>%
  pivot_longer(cols = c(Treatment_Mean, Control_Mean),
               names_to = "Group", values_to = "Mean") %>%
  mutate(Group = ifelse(Group == "Treatment_Mean", "GGW", "Non-GGW"),
         Variable = factor(Variable, levels = c("Climate Effect", "Vegetation Effect",
                                                "Interaction Effect", "Total Change")))

fig_bar <- ggplot(summary_bar, aes(x = Variable, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(values = c("GGW" = color_ggw, "Non-GGW" = color_nonggw)) +
  labs(title = sprintf("%s: Mean Effects by Group (50-threshold)", COUNTRY),
       subtitle = "Binary Richness - PSM 1:1 Matched",
       x = NULL, y = "Mean Effect (Species/pixel)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, paste0(COUNTRY, "_effects_barplot_50threshold.png")),
       fig_bar, width = 10, height = 6, dpi = 300, bg = "white")
cat("    ✓ \n")

# ---- 3: (ATT vs DR ATT) ----
cat("  ...\n")

# Binary
binary_dr <- dr_results[dr_results$Category == "Binary", ]
if(nrow(binary_dr) > 0 && any(!is.na(binary_dr$DR_ATT))) {
  
  comparison_data <- binary_dr %>%
    dplyr::select(Variable, Simple_ATT, Simple_p, DR_ATT, DR_p) %>%
    pivot_longer(cols = c(Simple_ATT, DR_ATT),
                 names_to = "Method", values_to = "ATT") %>%
    mutate(
      p_value = ifelse(Method == "Simple_ATT", Simple_p, DR_p),
      SE = NA,
      Method = ifelse(Method == "Simple_ATT", "Simple PSM", "Doubly Robust"),
      Method = factor(Method, levels = c("Simple PSM", "Doubly Robust"))
    ) %>%
    filter(!is.na(ATT))
  
  fig_dr_comparison <- ggplot(comparison_data, aes(x = ATT, y = Variable, color = Method)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3.5, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = c("Simple PSM" = "#E74C3C", "Doubly Robust" = "#2ECC71")) +
    labs(title = sprintf("%s: Simple PSM vs Doubly Robust (50-threshold)", COUNTRY),
         subtitle = sprintf("Climate adjustment: %s", paste(climate_vars, collapse = ", ")),
         x = "ATT (Average Treatment Effect on Treated)",
         y = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(output_dir, paste0(COUNTRY, "_DR_comparison_50threshold.png")),
         fig_dr_comparison, width = 10, height = 5, dpi = 300, bg = "white")
  cat("    ✓ \n")
}

# =============================================================================

# =============================================================================

cat("\n")
cat(rep("═", 70), "\n", sep = "")
cat(sprintf("  %s 50！\n", COUNTRY))
cat(rep("═", 70), "\n", sep = "")

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  : %d\n", nrow(matched_data)))
cat(sprintf("    - GGW (Treatment=1): %d\n", sum(matched_data$treatment == 1)))
cat(sprintf("    - Non-GGW (Treatment=0): %d\n", sum(matched_data$treatment == 0)))
if(length(available_climate) > 0) {
  cat(sprintf("  : %s\n", paste(available_climate, collapse = ", ")))
}

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  - %s_ATT_scenario_simple_50threshold.csv  (ATT)\n", COUNTRY))
cat(sprintf("  - %s_DR_scenario_ATT_50threshold.csv      (ATT)\n", COUNTRY))
cat(sprintf("  - %s_scenario_boxplot_50threshold.png      ()\n", COUNTRY))
cat(sprintf("  - %s_effects_barplot_50threshold.png       ()\n", COUNTRY))
cat(sprintf("  - %s_DR_comparison_50threshold.png         (DR)\n", COUNTRY))
cat(sprintf("\n: %s\n", output_dir))

cat(sprintf("\n=== %s ！===\n", COUNTRY))
cat(sprintf(":  COUNTRY <- \"%s\" 4\n",
            ifelse(COUNTRY == "Senegal", "Ethiopia", "Senegal")))
