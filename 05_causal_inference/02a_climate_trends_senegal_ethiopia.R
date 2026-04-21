# =============================================================================
# (Δbio)

# 2000-2024BIO

# : Δbio2, Δbio4, Δbio13, Δbio18
# : Δbio3, Δbio4, Δbio7, Δbio9, Δbio13, Δbio19
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(exactextractr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║         (Δbio) -  &                 ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. 
# =============================================================================

# ----- -----
sen_bio_dir <- "E:/GGW_bird_analysis/ENV/Senegal"
sen_psm_file <- "E:/2026.1.8_biomod2/senegal_results/PSM_1to1/Senegal_1to1_matched_data.csv"
sen_boundary_file <- "E:/11.17progress/study_area/SEN_adm4.shp"
sen_output_dir <- "E:/2026.1.8_biomod2/senegal_results/climate_trends"

# ----- -----
eth_bio_dir <- "E:/GGW_bird_analysis/ENV/Ethiopia"
eth_psm_file <- "E:/2026.1.8_biomod2/ethiopia_results/PSM_1to1/Ethiopia_1to1_matched_data.csv"
eth_boundary_file <- "E:/11.17progress/study_area/ETH_adm3.shp"
eth_output_dir <- "E:/2026.1.8_biomod2/ethiopia_results/climate_trends"


dir.create(sen_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(eth_output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. 
# =============================================================================

# : Biomod2bio2, bio4, bio13, bio18
sen_bio_vars <- c(2, 4, 13, 18)

# : Biomod2bio3, bio4, bio7, bio9, bio13, bio19
eth_bio_vars <- c(3, 4, 7, 9, 13, 19)


years <- 2000:2024

# =============================================================================
# 3. 
# =============================================================================

#' bio()
#' 
#' @param bio_dir BIO
#' @param boundary_sf sf
#' @param bio_num bio(1-19)
#' @param years 
#' @param id_col ID
#' @return data.frame 
extract_bio_timeseries <- function(bio_dir, boundary_sf, bio_num, years, id_col) {
  
  cat(sprintf("     bio%d ...\n", bio_num))
  

  result <- data.frame(ID = boundary_sf[[id_col]])
  names(result)[1] <- id_col
  
  for(year in years) {
    # - TIFF
    bio_file <- file.path(bio_dir, paste0("BIO_", year, ".tif"))
    
    if(!file.exists(bio_file)) {
      cat(sprintf("      :  %s\n", bio_file))
      result[[paste0("y", year)]] <- NA
      next
    }
    

    r <- tryCatch({
      rast(bio_file)[[bio_num]]
    }, error = function(e) {
      cat(sprintf("       %s: %s\n", bio_file, e$message))
      return(NULL)
    })
    
    if(is.null(r)) {
      result[[paste0("y", year)]] <- NA
      next
    }
    

    values <- exact_extract(r, boundary_sf, fun = "mean")
    result[[paste0("y", year)]] <- values
  }
  
  return(result)
}

#' 
#' 
#' @param values ()
#' @param years 
#' @return ()
calc_trend_slope <- function(values, years) {
  # NA
  valid <- !is.na(values)
  if(sum(valid) < 3) return(NA)  # 3
  

  fit <- lm(values[valid] ~ years[valid])
  return(coef(fit)[2])
}

#' bio
#' 
#' @param bio_dir BIO
#' @param boundary_sf 
#' @param bio_vars bio
#' @param years 
#' @param id_col ID
#' @return data.frame bio
calculate_all_trends <- function(bio_dir, boundary_sf, bio_vars, years, id_col) {
  
  cat("  ...\n")
  cat(sprintf("  : bio%s\n", paste(bio_vars, collapse = ", bio")))
  cat(sprintf("  : %d-%d (%d)\n", min(years), max(years), length(years)))
  cat(sprintf("  : %d\n\n", nrow(boundary_sf)))
  

  result <- data.frame(ID = boundary_sf[[id_col]])
  names(result)[1] <- id_col
  
  # bio
  for(bio_num in bio_vars) {
    

    ts_data <- extract_bio_timeseries(bio_dir, boundary_sf, bio_num, years, id_col)
    

    year_cols <- paste0("y", years)
    trends <- apply(ts_data[, year_cols], 1, function(row) {
      calc_trend_slope(as.numeric(row), years)
    })
    

    result[[paste0("bio", bio_num, "_trend")]] <- trends
    
    cat(sprintf("    bio%d :  = %.6f,  = [%.6f, %.6f]\n",
                bio_num, 
                mean(trends, na.rm = TRUE),
                min(trends, na.rm = TRUE),
                max(trends, na.rm = TRUE)))
  }
  
  return(result)
}

# =============================================================================
# 4. 
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("                                                            \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")


cat("  ...\n")
sen_boundary <- st_read(sen_boundary_file, quiet = TRUE)
cat(sprintf("  : %d\n", nrow(sen_boundary)))

# ID()
sen_id_col <- intersect(c("GID_4", "GID_3", "GID_2", "NAME_4", "NAME_3"), names(sen_boundary))[1]
if(is.na(sen_id_col)) {
  cat("  : ID\n")
  sen_id_col <- names(sen_boundary)[1]
}
cat(sprintf("  ID: %s\n\n", sen_id_col))


sen_trends <- calculate_all_trends(
  bio_dir = sen_bio_dir,
  boundary_sf = sen_boundary,
  bio_vars = sen_bio_vars,
  years = years,
  id_col = sen_id_col
)

# PSM
cat("\n  PSM...\n")
sen_psm <- read.csv(sen_psm_file)
cat(sprintf("  PSM: %d\n", nrow(sen_psm)))

# PSMID
sen_psm_id_col <- intersect(c("GID_4", "GID_3", "GID_2", "NAME_4", "NAME_3", sen_id_col), names(sen_psm))[1]
cat(sprintf("  PSMID: %s\n", sen_psm_id_col))


if(sen_id_col != sen_psm_id_col) {
  names(sen_trends)[1] <- sen_psm_id_col
}

sen_psm_with_trends <- sen_psm %>%
  left_join(sen_trends, by = sen_psm_id_col)


n_matched <- sum(!is.na(sen_psm_with_trends$bio2_trend))
cat(sprintf("  : %d / %d (%.1f%%)\n", 
            n_matched, nrow(sen_psm_with_trends), 
            100 * n_matched / nrow(sen_psm_with_trends)))


write.csv(sen_trends, file.path(sen_output_dir, "Senegal_climate_trends_all_LGA.csv"), row.names = FALSE)
write.csv(sen_psm_with_trends, file.path(sen_output_dir, "Senegal_PSM_with_climate_trends.csv"), row.names = FALSE)

cat("\n  :\n")
cat(sprintf("    - %s\n", file.path(sen_output_dir, "Senegal_climate_trends_all_LGA.csv")))
cat(sprintf("    - %s\n", file.path(sen_output_dir, "Senegal_PSM_with_climate_trends.csv")))


cat("\n   (PSM):\n")
cat("  ─────────────────────────────────────────────────────────────────\n")
for(bio_num in sen_bio_vars) {
  col_name <- paste0("bio", bio_num, "_trend")
  if(col_name %in% names(sen_psm_with_trends)) {
    treat_mean <- mean(sen_psm_with_trends[[col_name]][sen_psm_with_trends$treatment == 1], na.rm = TRUE)
    ctrl_mean <- mean(sen_psm_with_trends[[col_name]][sen_psm_with_trends$treatment == 0], na.rm = TRUE)
    cat(sprintf("  bio%d_trend: GGW = %.6f, Non-GGW = %.6f,  = %.6f\n",
                bio_num, treat_mean, ctrl_mean, treat_mean - ctrl_mean))
  }
}

# =============================================================================
# 5. 
# =============================================================================

cat("\n\n═══════════════════════════════════════════════════════════════════════\n")
cat("                                                           \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")


cat("  ...\n")
eth_boundary <- st_read(eth_boundary_file, quiet = TRUE)
cat(sprintf("  : %d\n", nrow(eth_boundary)))

# ID
eth_id_col <- intersect(c("GID_3", "GID_2", "NAME_3", "NAME_2"), names(eth_boundary))[1]
if(is.na(eth_id_col)) {
  cat("  : ID\n")
  eth_id_col <- names(eth_boundary)[1]
}
cat(sprintf("  ID: %s\n\n", eth_id_col))


eth_trends <- calculate_all_trends(
  bio_dir = eth_bio_dir,
  boundary_sf = eth_boundary,
  bio_vars = eth_bio_vars,
  years = years,
  id_col = eth_id_col
)

# PSM
cat("\n  PSM...\n")
eth_psm <- read.csv(eth_psm_file)
cat(sprintf("  PSM: %d\n", nrow(eth_psm)))

# PSMID
eth_psm_id_col <- intersect(c("GID_3", "GID_2", "NAME_3", "NAME_2", eth_id_col), names(eth_psm))[1]
cat(sprintf("  PSMID: %s\n", eth_psm_id_col))


if(eth_id_col != eth_psm_id_col) {
  names(eth_trends)[1] <- eth_psm_id_col
}

eth_psm_with_trends <- eth_psm %>%
  left_join(eth_trends, by = eth_psm_id_col)


n_matched <- sum(!is.na(eth_psm_with_trends$bio3_trend))
cat(sprintf("  : %d / %d (%.1f%%)\n", 
            n_matched, nrow(eth_psm_with_trends), 
            100 * n_matched / nrow(eth_psm_with_trends)))


write.csv(eth_trends, file.path(eth_output_dir, "Ethiopia_climate_trends_all_LGA.csv"), row.names = FALSE)
write.csv(eth_psm_with_trends, file.path(eth_output_dir, "Ethiopia_PSM_with_climate_trends.csv"), row.names = FALSE)

cat("\n  :\n")
cat(sprintf("    - %s\n", file.path(eth_output_dir, "Ethiopia_climate_trends_all_LGA.csv")))
cat(sprintf("    - %s\n", file.path(eth_output_dir, "Ethiopia_PSM_with_climate_trends.csv")))


cat("\n   (PSM):\n")
cat("  ─────────────────────────────────────────────────────────────────\n")
for(bio_num in eth_bio_vars) {
  col_name <- paste0("bio", bio_num, "_trend")
  if(col_name %in% names(eth_psm_with_trends)) {
    treat_mean <- mean(eth_psm_with_trends[[col_name]][eth_psm_with_trends$treatment == 1], na.rm = TRUE)
    ctrl_mean <- mean(eth_psm_with_trends[[col_name]][eth_psm_with_trends$treatment == 0], na.rm = TRUE)
    cat(sprintf("  bio%d_trend: GGW = %.6f, Non-GGW = %.6f,  = %.6f\n",
                bio_num, treat_mean, ctrl_mean, treat_mean - ctrl_mean))
  }
}

# =============================================================================
# 6. 
# =============================================================================

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                          ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("【】\n")
cat(sprintf("  : Δbio%s\n", paste(sen_bio_vars, collapse = ", Δbio")))
cat(sprintf("  : %s\n", sen_output_dir))
cat("  :\n")
cat("    1. Senegal_climate_trends_all_LGA.csv - \n")
cat("    2. Senegal_PSM_with_climate_trends.csv - PSM+ ()\n")

cat("\n【】\n")
cat(sprintf("  : Δbio%s\n", paste(eth_bio_vars, collapse = ", Δbio")))
cat(sprintf("  : %s\n", eth_output_dir))
cat("  :\n")
cat("    1. Ethiopia_climate_trends_all_LGA.csv - \n")
cat("    2. Ethiopia_PSM_with_climate_trends.csv - PSM+ ()\n")

cat("\n【】\n")
cat("   *_PSM_with_climate_trends.csv \n")

cat("\n=== ！===\n")



# =============================================================================
# - & 

# PSM

# : Δbio2, Δbio4, Δbio13, Δbio18
# : Δbio3, Δbio4, Δbio7, Δbio9, Δbio13, Δbio19
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║            -  &                         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. 
# =============================================================================

# ----- -----
sen_data_file <- "E:/2026.1.8_biomod2/senegal_results/climate_trends/Senegal_PSM_with_climate_trends.csv"
sen_scenario_file <- "E:/2026.1.8_biomod2/senegal_results/scenario_decomposition/matched_LGA_effects.csv"
sen_output_dir <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_results"

# ----- -----
eth_data_file <- "E:/2026.1.8_biomod2/ethiopia_results/climate_trends/Ethiopia_PSM_with_climate_trends.csv"
eth_scenario_file <- "E:/2026.1.8_biomod2/ethiopia_results/scenario_decomposition/matched_LGA_effects.csv"
eth_output_dir <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_results"


dir.create(sen_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(eth_output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. 
# =============================================================================


sen_climate_vars <- c("bio2_trend", "bio4_trend", "bio13_trend", "bio18_trend")


eth_climate_vars <- c("bio3_trend", "bio4_trend", "bio7_trend", "bio9_trend", "bio13_trend", "bio19_trend")

# ()
outcome_vars <- c("bin_change", "prob_change", "tss_change")
outcome_labels <- c("Binary ()", "Probability (%)", "TSS-weighted")

# =============================================================================
# 3. 
# =============================================================================

#' SMD ()
calc_smd <- function(data, var_name, treat_col = "treatment") {
  treat <- data %>% filter(!!sym(treat_col) == 1)
  ctrl <- data %>% filter(!!sym(treat_col) == 0)
  
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
    Balanced = ifelse(abs(smd) < 0.25, "Yes", "No")
  ))
}

#' ATT
#' 
#' @param data PSM()
#' @param outcome 
#' @param climate_vars 
#' @param weight_col ()
#' @return list simpledoubly robust
doubly_robust_att <- function(data, outcome, climate_vars, weight_col = NULL) {
  

  if(!(outcome %in% names(data))) {
    return(list(simple = NULL, dr = NULL))
  }
  
  # NA
  valid_vars <- c("treatment", outcome, climate_vars)
  if(!is.null(weight_col)) valid_vars <- c(valid_vars, weight_col)
  data_clean <- data %>% 
    dplyr::select(all_of(valid_vars)) %>%
    na.omit()
  
  if(nrow(data_clean) < 10) {
    return(list(simple = NULL, dr = NULL))
  }
  
  # ----- Method 1: Simple PSM (t-test) -----
  treat <- data_clean %>% filter(treatment == 1)
  ctrl <- data_clean %>% filter(treatment == 0)
  
  t_result <- t.test(treat[[outcome]], ctrl[[outcome]])
  simple_att <- t_result$estimate[1] - t_result$estimate[2]
  simple_se <- sqrt(var(treat[[outcome]])/nrow(treat) + var(ctrl[[outcome]])/nrow(ctrl))
  
  simple_result <- data.frame(
    Method = "Simple PSM",
    ATT = simple_att,
    SE = simple_se,
    CI_Lower = t_result$conf.int[1],
    CI_Upper = t_result$conf.int[2],
    p_value = t_result$p.value
  )
  
  # ----- Method 2: Doubly Robust () -----

  available_climate_vars <- intersect(climate_vars, names(data_clean))
  
  if(length(available_climate_vars) == 0) {
    cat("    : \n")
    return(list(simple = simple_result, dr = NULL))
  }
  

  formula_str <- paste(outcome, "~ treatment +", paste(available_climate_vars, collapse = " + "))
  formula_dr <- as.formula(formula_str)
  

  if(!is.null(weight_col) && weight_col %in% names(data_clean)) {
    model_dr <- lm(formula_dr, data = data_clean, weights = data_clean[[weight_col]])
  } else {
    model_dr <- lm(formula_dr, data = data_clean)
  }
  
  coef_dr <- summary(model_dr)$coefficients
  
  dr_att <- coef_dr["treatment", "Estimate"]
  dr_se <- coef_dr["treatment", "Std. Error"]
  dr_p <- coef_dr["treatment", "Pr(>|t|)"]
  
  dr_result <- data.frame(
    Method = "Doubly Robust",
    ATT = dr_att,
    SE = dr_se,
    CI_Lower = dr_att - 1.96 * dr_se,
    CI_Upper = dr_att + 1.96 * dr_se,
    p_value = dr_p
  )
  
  return(list(
    simple = simple_result, 
    dr = dr_result,
    model = model_dr,
    climate_vars_used = available_climate_vars
  ))
}

#' 
run_doubly_robust_analysis <- function(data, climate_vars, outcome_vars, outcome_labels, 
                                       country_name, output_dir, weight_col = NULL) {
  
  cat(sprintf("\n═══════════════════════════════════════════════════════════════════════\n"))
  cat(sprintf("                    %s                                     \n", country_name))
  cat(sprintf("═══════════════════════════════════════════════════════════════════════\n\n"))
  
  # ----- -----
  cat("【】\n")
  cat("─────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-15s %12s %12s %8s %8s\n", "Variable", "GGW", "Non-GGW", "SMD", "Balanced"))
  cat("─────────────────────────────────────────────────────────────────────────\n")
  
  balance_results <- data.frame()
  for(var in climate_vars) {
    if(var %in% names(data)) {
      bal <- calc_smd(data, var)
      balance_results <- rbind(balance_results, bal)
      cat(sprintf("%-15s %12.6f %12.6f %8.3f %8s\n",
                  bal$Variable, bal$Mean_GGW, bal$Mean_NonGGW, bal$SMD, bal$Balanced))
    } else {
      cat(sprintf("%-15s %12s\n", var, ""))
    }
  }
  
  # ----- ATT-----
  cat("\n\n【ATT】\n")
  cat("─────────────────────────────────────────────────────────────────────────\n")
  
  att_results <- data.frame()
  
  for(i in seq_along(outcome_vars)) {
    outcome <- outcome_vars[i]
    label <- outcome_labels[i]
    
    cat(sprintf("\n[%s]\n", label))
    
    result <- doubly_robust_att(data, outcome, climate_vars, weight_col)
    
    if(is.null(result$simple)) {
      cat("  \n")
      next
    }
    
    # Simple PSM
    simple <- result$simple
    simple$Metric <- label
    simple$Sig <- ifelse(simple$p_value < 0.001, "***", 
                         ifelse(simple$p_value < 0.01, "**", 
                                ifelse(simple$p_value < 0.05, "*", "")))
    
    cat(sprintf("  Simple PSM:    ATT = %8.3f (SE = %.3f), 95%% CI [%.3f, %.3f], p %s %s\n",
                simple$ATT, simple$SE, simple$CI_Lower, simple$CI_Upper,
                ifelse(simple$p_value < 0.001, "< 0.001", sprintf("= %.4f", simple$p_value)),
                simple$Sig))
    
    att_results <- rbind(att_results, simple)
    
    # Doubly Robust
    if(!is.null(result$dr)) {
      dr <- result$dr
      dr$Metric <- label
      dr$Sig <- ifelse(dr$p_value < 0.001, "***", 
                       ifelse(dr$p_value < 0.01, "**", 
                              ifelse(dr$p_value < 0.05, "*", "")))
      
      cat(sprintf("  Doubly Robust: ATT = %8.3f (SE = %.3f), 95%% CI [%.3f, %.3f], p %s %s\n",
                  dr$ATT, dr$SE, dr$CI_Lower, dr$CI_Upper,
                  ifelse(dr$p_value < 0.001, "< 0.001", sprintf("= %.4f", dr$p_value)),
                  dr$Sig))
      

      pct_change <- (dr$ATT - simple$ATT) / abs(simple$ATT) * 100
      cat(sprintf("  → ATT: %.1f%%\n", pct_change))
      
      att_results <- rbind(att_results, dr)
    }
  }
  
  # ----- -----
  cat("\n\n【】\n")
  cat("─────────────────────────────────────────────────────────────────────────\n")
  
  # ATT
  att_file <- file.path(output_dir, paste0(country_name, "_ATT_doubly_robust.csv"))
  write.csv(att_results, att_file, row.names = FALSE)
  cat(sprintf("  - %s\n", att_file))
  

  if(nrow(balance_results) > 0) {
    balance_file <- file.path(output_dir, paste0(country_name, "_climate_balance.csv"))
    write.csv(balance_results, balance_file, row.names = FALSE)
    cat(sprintf("  - %s\n", balance_file))
  }
  
  return(list(
    att_results = att_results,
    balance_results = balance_results
  ))
}

# =============================================================================
# 4. 
# =============================================================================

if(file.exists(sen_data_file)) {
  
  sen_data <- read.csv(sen_data_file)
  cat(sprintf(": %d \n", nrow(sen_data)))
  

  cat(":\n")
  print(names(sen_data))
  

  sen_results <- run_doubly_robust_analysis(
    data = sen_data,
    climate_vars = sen_climate_vars,
    outcome_vars = outcome_vars,
    outcome_labels = outcome_labels,
    country_name = "Senegal",
    output_dir = sen_output_dir,
    weight_col = "weights"
  )
  
} else {
  cat(sprintf("\n: : %s\n", sen_data_file))
  cat(" 01_calculate_climate_trends.R \n")
}

# =============================================================================
# 5. 
# =============================================================================

if(file.exists(eth_data_file)) {
  
  eth_data <- read.csv(eth_data_file)
  cat(sprintf("\n: %d \n", nrow(eth_data)))
  

  cat(":\n")
  print(names(eth_data))
  

  eth_results <- run_doubly_robust_analysis(
    data = eth_data,
    climate_vars = eth_climate_vars,
    outcome_vars = outcome_vars,
    outcome_labels = outcome_labels,
    country_name = "Ethiopia",
    output_dir = eth_output_dir,
    weight_col = "weights"
  )
  
} else {
  cat(sprintf("\n: : %s\n", eth_data_file))
  cat(" 01_calculate_climate_trends.R \n")
}

# =============================================================================
# 6. 
# =============================================================================

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")


if(exists("sen_results") && exists("eth_results")) {
  

  all_att <- rbind(
    sen_results$att_results %>% mutate(Country = "Senegal"),
    eth_results$att_results %>% mutate(Country = "Ethiopia")
  )
  

  if(nrow(all_att) > 0) {
    
    fig_comparison <- ggplot(all_att, aes(x = ATT, y = Metric, color = Method, shape = Country)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), 
                     height = 0.2, linewidth = 0.8,
                     position = position_dodge(width = 0.6)) +
      geom_point(size = 3.5, position = position_dodge(width = 0.6)) +
      scale_color_manual(values = c("Simple PSM" = "#E74C3C", "Doubly Robust" = "#2ECC71")) +
      facet_wrap(~Country, ncol = 1) +
      labs(title = "GGW Effect on Species Richness: Simple PSM vs Doubly Robust",
           subtitle = "Doubly Robust = PSM + Climate Trend Adjustment",
           x = "ATT (Average Treatment Effect on Treated)",
           y = NULL) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = "#3498DB"),
            strip.text = element_text(color = "white", face = "bold"))
    

    comparison_file <- "E:/2026.1.8_biomod2/ATT_comparison_all_countries.png"
    ggsave(comparison_file, fig_comparison, width = 12, height = 8, dpi = 300, bg = "white")
    cat(sprintf(": %s\n", comparison_file))
  }
}

# =============================================================================
# 7. 
# =============================================================================

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                          ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  :   Δbio1, Δbio12\n")
cat("  :   Δbio2, Δbio4, Δbio13, Δbio18\n")
cat("  : Δbio3, Δbio4, Δbio7, Δbio9, Δbio13, Δbio19\n")

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  (Biomod2)\n")
cat("  \n")

cat("\n=== ！===\n")


# =============================================================================
# & ()


# - : ID=GID_4, treatment
# - : ID=GID_3, treatment
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║   -  &  ()     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. 
# =============================================================================

# ----- -----
sen_psm_file <- "E:/2026.1.8_biomod2/senegal_results/PSM_1to1/Senegal_1to1_matched_data.csv"
sen_effects_file <- "E:/2026.1.8_biomod2/senegal_results/scenario_decomposition/matched_LGA_effects.csv"
sen_climate_file <- "E:/2026.1.8_biomod2/senegal_results/climate_trends/Senegal_PSM_with_climate_trends.csv"
sen_output_dir <- "E:/2026.1.8_biomod2/senegal_results/doubly_robust_results"

# ----- -----
eth_psm_file <- "E:/2026.1.8_biomod2/ethiopia_results/PSM_1to1/Ethiopia_1to1_matched_data.csv"
eth_effects_file <- "E:/2026.1.8_biomod2/ethiopia_results/scenario_decomposition/matched_Woreda_effects.csv"
eth_climate_file <- "E:/2026.1.8_biomod2/ethiopia_results/climate_trends/Ethiopia_PSM_with_climate_trends.csv"
eth_output_dir <- "E:/2026.1.8_biomod2/ethiopia_results/doubly_robust_results"


dir.create(sen_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(eth_output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. 
# =============================================================================

# : bio2, bio4, bio13, bio18
sen_climate_vars <- c("bio2_trend", "bio4_trend", "bio13_trend", "bio18_trend")
sen_id_col <- "GID_4"

# : bio3, bio4, bio7, bio9, bio13, bio19
eth_climate_vars <- c("bio3_trend", "bio4_trend", "bio7_trend", "bio9_trend", "bio13_trend", "bio19_trend")
eth_id_col <- "GID_3"

# =============================================================================
# 3. 
# =============================================================================

add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

analyze_effect_DR <- function(data, effect_var, effect_name, climate_vars, treat_col = "treatment") {
  
  if(!(effect_var %in% names(data))) {
    cat(sprintf("    :  %s \n", effect_var))
    return(NULL)
  }
  
  if(!(treat_col %in% names(data))) {
    cat(sprintf("    : treatment %s \n", treat_col))
    return(NULL)
  }
  
  available_climate <- intersect(climate_vars, names(data))
  cols_needed <- c(treat_col, effect_var, available_climate)
  
  data_clean <- data[, cols_needed, drop = FALSE]
  data_clean <- na.omit(data_clean)
  
  if(nrow(data_clean) < 20) {
    cat(sprintf("    : %s  (%d)\n", effect_var, nrow(data_clean)))
    return(NULL)
  }
  
  treat <- data_clean[data_clean[[treat_col]] == 1, ]
  ctrl <- data_clean[data_clean[[treat_col]] == 0, ]
  

  t_result <- t.test(treat[[effect_var]], ctrl[[effect_var]])
  simple_diff <- mean(treat[[effect_var]], na.rm=TRUE) - mean(ctrl[[effect_var]], na.rm=TRUE)
  simple_p <- t_result$p.value
  

  if(length(available_climate) == 0) {
    return(data.frame(
      Effect = effect_name,
      Simple_Diff = simple_diff, Simple_p = simple_p,
      Adjusted_Diff = NA, Adjusted_SE = NA, Adjusted_p = NA
    ))
  }
  
  # treatment
  names(data_clean)[names(data_clean) == treat_col] <- "treatment_var"
  
  formula_adj <- as.formula(paste(effect_var, "~ treatment_var +", 
                                  paste(available_climate, collapse = " + ")))
  model_adj <- lm(formula_adj, data = data_clean)
  coef_adj <- summary(model_adj)$coefficients
  
  return(data.frame(
    Effect = effect_name,
    Simple_Diff = simple_diff,
    Simple_p = simple_p,
    Adjusted_Diff = coef_adj["treatment_var", "Estimate"],
    Adjusted_SE = coef_adj["treatment_var", "Std. Error"],
    Adjusted_p = coef_adj["treatment_var", "Pr(>|t|)"]
  ))
}

run_scenario_decomposition_DR <- function(psm_file, effects_file, climate_file,
                                          climate_vars, id_col, country_name, output_dir) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("                    %s                               \n", country_name))
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # ----- PSM-----
  if(!file.exists(psm_file)) {
    cat(sprintf(": PSM: %s\n", psm_file))
    return(NULL)
  }
  psm_data <- read.csv(psm_file)
  cat(sprintf("  PSM: %d \n", nrow(psm_data)))
  
  # ----- -----
  if(!file.exists(effects_file)) {
    cat(sprintf(": : %s\n", effects_file))
    return(NULL)
  }
  effects_data <- read.csv(effects_file)
  cat(sprintf("  : %d \n", nrow(effects_data)))
  
  # ----- -----
  climate_data <- NULL
  if(file.exists(climate_file)) {
    climate_data <- read.csv(climate_file)
    cat(sprintf("  : %d \n", nrow(climate_data)))
    available_climate <- intersect(climate_vars, names(climate_data))
    cat(sprintf("  : %s\n", paste(available_climate, collapse = ", ")))
    
    # ID
    climate_cols <- c(id_col, available_climate)
    climate_data <- climate_data[, intersect(climate_cols, names(climate_data)), drop = FALSE]
    climate_data <- climate_data[!duplicated(climate_data[[id_col]]), ]
  } else {
    cat(sprintf(": : %s\n", climate_file))
    available_climate <- c()
  }
  
  # ----- -----
  cat("\n【】\n")
  
  # treatment
  effects_cols <- c(id_col, "binary_climate", "binary_vegetation", "binary_interaction",
                    "prob_climate", "prob_vegetation", "prob_interaction",
                    "tss_climate", "tss_vegetation", "tss_interaction")
  effects_subset <- effects_data[, intersect(effects_cols, names(effects_data)), drop = FALSE]
  
  # PSM + 
  merged_data <- merge(psm_data, effects_subset, by = id_col, all.x = TRUE)
  cat(sprintf("  PSM + : %d \n", nrow(merged_data)))
  
  # treatment
  if("treatment" %in% names(merged_data)) {
    cat(sprintf("  treatment:  (0: %d, 1: %d)\n", 
                sum(merged_data$treatment == 0), sum(merged_data$treatment == 1)))
  }
  

  if(!is.null(climate_data) && nrow(climate_data) > 0) {
    merged_data <- merge(merged_data, climate_data, by = id_col, all.x = TRUE)
    cat(sprintf("  + : %d \n", nrow(merged_data)))
    
    # NA
    for(cv in available_climate) {
      if(cv %in% names(merged_data)) {
        na_count <- sum(is.na(merged_data[[cv]]))
        cat(sprintf("    %s: %d  (%d NA)\n", cv, 
                    nrow(merged_data) - na_count, na_count))
      }
    }
  }
  
  # ----- -----
  cat("\n【】\n")
  
  effect_definitions <- list(
    list(var = "binary_climate", name = "Climate Effect", metric = "Binary"),
    list(var = "binary_vegetation", name = "Vegetation Effect", metric = "Binary"),
    list(var = "binary_interaction", name = "Interaction Effect", metric = "Binary"),
    list(var = "prob_climate", name = "Climate Effect", metric = "Probability"),
    list(var = "prob_vegetation", name = "Vegetation Effect", metric = "Probability"),
    list(var = "prob_interaction", name = "Interaction Effect", metric = "Probability"),
    list(var = "tss_climate", name = "Climate Effect", metric = "TSS-weighted"),
    list(var = "tss_vegetation", name = "Vegetation Effect", metric = "TSS-weighted"),
    list(var = "tss_interaction", name = "Interaction Effect", metric = "TSS-weighted")
  )
  
  all_results <- data.frame()
  
  for(ef in effect_definitions) {
    result <- analyze_effect_DR(merged_data, ef$var, ef$name, climate_vars, "treatment")
    if(!is.null(result)) {
      result$Metric <- ef$metric
      all_results <- rbind(all_results, result)
    }
  }
  
  # ----- -----
  if(nrow(all_results) > 0) {
    
    for(metric in c("Binary", "Probability", "TSS-weighted")) {
      metric_results <- all_results[all_results$Metric == metric, ]
      
      if(nrow(metric_results) > 0) {
        cat(sprintf("\n【%s - 】\n", metric))
        cat("─────────────────────────────────────────────────────────────────────────\n")
        cat(sprintf("%-20s %12s %12s %20s\n", "", "", "", ""))
        cat("─────────────────────────────────────────────────────────────────────────\n")
        
        for(i in 1:nrow(metric_results)) {
          row <- metric_results[i, ]
          sig_simple <- add_sig(row$Simple_p)
          sig_adj <- ifelse(is.na(row$Adjusted_p), "", add_sig(row$Adjusted_p))
          

          if(!is.na(row$Adjusted_Diff) && abs(row$Simple_Diff) > 0.001) {
            pct_change <- (row$Adjusted_Diff - row$Simple_Diff) / abs(row$Simple_Diff) * 100
          } else {
            pct_change <- NA
          }
          

          if(grepl("Climate", row$Effect)) {
            if(is.na(row$Adjusted_p)) {
              interpretation <- ""
            } else if(row$Simple_p < 0.05 && row$Adjusted_p > 0.05) {
              interpretation <- "✓ "
            } else if(row$Simple_p > 0.05 && row$Adjusted_p > 0.05) {
              interpretation <- "○ "
            } else if(!is.na(pct_change) && abs(pct_change) < 15) {
              interpretation <- ""
            } else {
              interpretation <- ""
            }
          } else if(grepl("Vegetation", row$Effect)) {
            if(is.na(row$Adjusted_p)) {
              interpretation <- ""
            } else {
              dir_str <- ifelse(row$Adjusted_Diff > 0, "+", "-")
              interpretation <- ifelse(row$Adjusted_p < 0.05, 
                                       sprintf("GGW(%s)", dir_str), 
                                       "")
            }
          } else {
            if(is.na(row$Adjusted_p)) {
              interpretation <- ""
            } else {
              interpretation <- ifelse(row$Adjusted_p < 0.05, "", "")
            }
          }
          
          adj_str <- ifelse(is.na(row$Adjusted_Diff), "NA", sprintf("%.2f", row$Adjusted_Diff))
          
          cat(sprintf("%-20s %10.2f%s %10s%s %s\n",
                      row$Effect,
                      row$Simple_Diff, sig_simple,
                      adj_str, sig_adj,
                      interpretation))
        }
        cat("─────────────────────────────────────────────────────────────────────────\n")
      }
    }
    cat("* p<0.05, ** p<0.01, *** p<0.001\n")
    
    # ----- -----
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════════════╗\n")
    cat(sprintf("║  %s                                                     ║\n", country_name))
    cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")
    
    binary_climate <- all_results[all_results$Metric == "Binary" & all_results$Effect == "Climate Effect", ]
    binary_veg <- all_results[all_results$Metric == "Binary" & all_results$Effect == "Vegetation Effect", ]
    
    if(nrow(binary_climate) > 0 && !is.na(binary_climate$Adjusted_p)) {
      pct_change <- (binary_climate$Adjusted_Diff - binary_climate$Simple_Diff) / 
        abs(binary_climate$Simple_Diff) * 100
      
      cat("  【】\n")
      cat(sprintf("    : %.2f%s\n", binary_climate$Simple_Diff, add_sig(binary_climate$Simple_p)))
      cat(sprintf("    : %.2f%s\n", binary_climate$Adjusted_Diff, add_sig(binary_climate$Adjusted_p)))
      cat(sprintf("    : %.1f%%\n", pct_change))
      
      if(binary_climate$Simple_p < 0.05 && binary_climate$Adjusted_p > 0.05) {
        cat("    → ✓ ！\n")
      } else if(binary_climate$Simple_p > 0.05) {
        cat("    → ○ \n")
      } else if(abs(pct_change) < 15) {
        cat("    → \n")
      } else {
        cat("    → \n")
      }
    }
    
    if(nrow(binary_veg) > 0 && !is.na(binary_veg$Adjusted_p)) {
      cat("\n  【 (GGW)】\n")
      cat(sprintf("    : %.2f%s\n", binary_veg$Adjusted_Diff, add_sig(binary_veg$Adjusted_p)))
      
      if(binary_veg$Adjusted_p < 0.05) {
        if(binary_veg$Adjusted_Diff > 0) {
          cat("    → ✓ GGW\n")
        } else {
          cat("    → ⚠ GGW\n")
        }
      } else {
        cat("    → GGW\n")
      }
    }
    
    # ----- -----
    all_results$Country <- country_name
    output_file <- file.path(output_dir, paste0(country_name, "_scenario_decomposition_DR.csv"))
    write.csv(all_results, output_file, row.names = FALSE)
    cat(sprintf("\n  : %s\n", output_file))
    
    return(all_results)
  } else {
    cat("  \n")
    return(NULL)
  }
}

# =============================================================================
# 4. 
# =============================================================================

sen_results <- run_scenario_decomposition_DR(
  psm_file = sen_psm_file,
  effects_file = sen_effects_file,
  climate_file = sen_climate_file,
  climate_vars = sen_climate_vars,
  id_col = sen_id_col,
  country_name = "Senegal",
  output_dir = sen_output_dir
)

# =============================================================================
# 5. 
# =============================================================================

eth_results <- run_scenario_decomposition_DR(
  psm_file = eth_psm_file,
  effects_file = eth_effects_file,
  climate_file = eth_climate_file,
  climate_vars = eth_climate_vars,
  id_col = eth_id_col,
  country_name = "Ethiopia",
  output_dir = eth_output_dir
)

# =============================================================================
# 6. 
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("【 (Binary) - 】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-12s %12s %12s %12s %15s\n", 
            "", "", "", "%", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")

# ()
cat(sprintf("%-12s %10.2f*** %10.2f %10.1f%% %15s\n",
            "Nigeria", 10.70, 1.90, -82.2, "✓ "))


if(!is.null(sen_results)) {
  sen_climate <- sen_results[sen_results$Metric == "Binary" & sen_results$Effect == "Climate Effect", ]
  if(nrow(sen_climate) > 0 && !is.na(sen_climate$Adjusted_p)) {
    pct <- (sen_climate$Adjusted_Diff - sen_climate$Simple_Diff) / abs(sen_climate$Simple_Diff) * 100
    sig_b <- add_sig(sen_climate$Simple_p)
    sig_a <- add_sig(sen_climate$Adjusted_p)
    
    if(sen_climate$Simple_p < 0.05 && sen_climate$Adjusted_p > 0.05) {
      conclusion <- "✓ "
    } else if(sen_climate$Simple_p > 0.05) {
      conclusion <- "○ "
    } else if(abs(pct) < 15) {
      conclusion <- ""
    } else {
      conclusion <- ""
    }
    
    cat(sprintf("%-12s %10.2f%s %10.2f%s %10.1f%% %15s\n",
                "Senegal", 
                sen_climate$Simple_Diff, sig_b,
                sen_climate$Adjusted_Diff, sig_a,
                pct, conclusion))
  }
}


if(!is.null(eth_results)) {
  eth_climate <- eth_results[eth_results$Metric == "Binary" & eth_results$Effect == "Climate Effect", ]
  if(nrow(eth_climate) > 0 && !is.na(eth_climate$Adjusted_p)) {
    pct <- (eth_climate$Adjusted_Diff - eth_climate$Simple_Diff) / abs(eth_climate$Simple_Diff) * 100
    sig_b <- add_sig(eth_climate$Simple_p)
    sig_a <- add_sig(eth_climate$Adjusted_p)
    
    if(eth_climate$Simple_p < 0.05 && eth_climate$Adjusted_p > 0.05) {
      conclusion <- "✓ "
    } else if(eth_climate$Simple_p > 0.05) {
      conclusion <- "○ "
    } else if(abs(pct) < 15) {
      conclusion <- ""
    } else {
      conclusion <- ""
    }
    
    cat(sprintf("%-12s %10.2f%s %10.2f%s %10.1f%% %15s\n",
                "Ethiopia", 
                eth_climate$Simple_Diff, sig_b,
                eth_climate$Adjusted_Diff, sig_a,
                pct, conclusion))
  }
}

cat("─────────────────────────────────────────────────────────────────────────\n")


cat("\n【 (Binary) - GGW】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-12s %15s %15s\n", "", "", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")


cat(sprintf("%-12s %13.2f*** %15s\n", "Nigeria", 7.02, "GGW"))


if(!is.null(sen_results)) {
  sen_veg <- sen_results[sen_results$Metric == "Binary" & sen_results$Effect == "Vegetation Effect", ]
  if(nrow(sen_veg) > 0 && !is.na(sen_veg$Adjusted_p)) {
    sig <- add_sig(sen_veg$Adjusted_p)
    if(sen_veg$Adjusted_p < 0.05) {
      conclusion <- ifelse(sen_veg$Adjusted_Diff > 0, "GGW", "GGW")
    } else {
      conclusion <- ""
    }
    cat(sprintf("%-12s %13.2f%s %15s\n", "Senegal", sen_veg$Adjusted_Diff, sig, conclusion))
  }
}


if(!is.null(eth_results)) {
  eth_veg <- eth_results[eth_results$Metric == "Binary" & eth_results$Effect == "Vegetation Effect", ]
  if(nrow(eth_veg) > 0 && !is.na(eth_veg$Adjusted_p)) {
    sig <- add_sig(eth_veg$Adjusted_p)
    if(eth_veg$Adjusted_p < 0.05) {
      conclusion <- ifelse(eth_veg$Adjusted_Diff > 0, "GGW", "GGW")
    } else {
      conclusion <- ""
    }
    cat(sprintf("%-12s %13.2f%s %15s\n", "Ethiopia", eth_veg$Adjusted_Diff, sig, conclusion))
  }
}

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n")

# =============================================================================
# 7. 
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  1. ✓ :  → \n")
cat("     → \n\n")

cat("  2. ○ : \n")
cat("     → PSM\n\n")

cat("  3. : <15%\n")
cat("     → \n")
cat("     → ATT\n\n")

cat("  4. : >15%\n")
cat("     → \n")
cat("─────────────────────────────────────────────────────────────────────────\n")

cat("\n【SMD ≠ 】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  : SMD=0.42 (), ATT-43% → \n")
cat("  : SMD=1.16 (), ATT-5%   → \n")
cat("  : SMD=0.31 (), ATT+6% → \n\n")
cat("  :  =  × \n")
cat("        \n")
cat("─────────────────────────────────────────────────────────────────────────\n")

cat("\n=== ！===\n")