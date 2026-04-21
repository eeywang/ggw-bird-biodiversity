# ============================================================================
# & (- )
# 1. BIO1-19(2007-2015, 2016-2024)
# 2. (DEM, Slope, Aspect_sin, Aspect_cos)
# 3. VIF+ 
# ============================================================================


library(terra)
library(usdm)      # VIF
library(corrplot)
library(dplyr)

# ============================================================================
# PART 0: 
# ============================================================================

# ============================================
# !
# ============================================

# ()
SEN_BIO_dir <- "E:/GGW_bird_analysis/ENV/Senegal"                    # BIO_2000.tif - BIO_2024.tif
SEN_ENV_dir <- "E:/2026.1.8_biomod2/Senegal_ENV_TwoPeriod"           # DHI, SHDI, VCF, gHM, 
SEN_output_dir <- "E:/2026.1.8_biomod2/Senegal_ENV_TwoPeriod"

# ()
ETH_BIO_dir <- "E:/GGW_bird_analysis/ENV/Ethiopia"                   # BIO_2000.tif - BIO_2024.tif
ETH_ENV_dir <- "E:/2026.1.8_biomod2/Ethiopia_ENV_TwoPeriod"          # DHI, SHDI, VCF, gHM, 
ETH_output_dir <- "E:/2026.1.8_biomod2/Ethiopia_ENV_TwoPeriod"


cat("...\n")
stopifnot("BIO" = dir.exists(SEN_BIO_dir))
stopifnot("ENV" = dir.exists(SEN_ENV_dir))
stopifnot("BIO" = dir.exists(ETH_BIO_dir))
stopifnot("ENV" = dir.exists(ETH_ENV_dir))
cat("!\n")

# ============================================================================
# PART 1: BIO
# ============================================================================

calculate_bio_average <- function(bio_dir, years, output_dir, country_prefix, period_name) {
  
  cat("\n========================================\n")
  cat(sprintf(" %s %s BIO\n", country_prefix, period_name))
  cat(sprintf(": %d - %d\n", min(years), max(years)))
  cat("========================================\n")
  
  # BIO
  bio_files <- list()
  for (year in years) {
    file_path <- file.path(bio_dir, sprintf("BIO_%d.tif", year))
    if (file.exists(file_path)) {
      bio_files[[as.character(year)]] <- rast(file_path)
      cat(sprintf("  : BIO_%d.tif (%d )\n", year, nlyr(rast(file_path))))
    } else {
      warning(sprintf(": %s", file_path))
    }
  }
  
  if (length(bio_files) == 0) {
    stop("BIO!")
  }
  
  # BIO19
  n_bio <- nlyr(bio_files[[1]])
  cat(sprintf("\nBIO: %d\n", n_bio))
  
  # BIO
  bio_avg_list <- list()
  
  for (i in 1:n_bio) {
    # iBIO
    layers <- lapply(bio_files, function(x) x[[i]])
    stacked <- rast(layers)
    

    avg <- mean(stacked, na.rm = TRUE)
    names(avg) <- sprintf("bio%d", i)
    bio_avg_list[[i]] <- avg
    
    cat(sprintf("  bio%d: \n", i))
  }
  

  bio_avg <- rast(bio_avg_list)
  

  output_file <- file.path(output_dir, sprintf("%s_BIO_avg_%s.tif", country_prefix, period_name))
  writeRaster(bio_avg, output_file, overwrite = TRUE)
  cat(sprintf("\n: %s\n", output_file))
  
  return(bio_avg)
}

# ============================================================================
# PART 2: (- )
# ============================================================================

load_and_stack_env <- function(env_dir, country_prefix, period_suffix, bio_avg) {
  
  cat("\n========================================\n")
  cat(sprintf(" %s %s \n", country_prefix, period_suffix))
  cat("========================================\n")
  
  # afterbefore
  is_after <- grepl("after", period_suffix)
  
  # (- )
  env_files_dynamic <- list(
    DHI_mean = sprintf("%s_DHI_mean_%s.tif", country_prefix, period_suffix),
    DHI_min = sprintf("%s_DHI_min_%s.tif", country_prefix, period_suffix),
    DHI_cv = sprintf("%s_DHI_cv_%s.tif", country_prefix, period_suffix),
    SHDI = sprintf("%s_SHDI_%s.tif", country_prefix, period_suffix),
    Woody_pct = sprintf("%s_Woody_pct_%s.tif", country_prefix, 
                        ifelse(is_after, "after_2016_2023", "before_2007_2015")),
    Herbaceous_pct = sprintf("%s_Herbaceous_pct_%s.tif", country_prefix,
                             ifelse(is_after, "after_2016_2023", "before_2007_2015"))
  )
  
  # ()
  env_files_static <- list(
    gHM = sprintf("%s_gHM_static.tif", country_prefix),
    DEM = sprintf("%s_DEM_static.tif", country_prefix),
    Slope = sprintf("%s_Slope_static.tif", country_prefix),
    Aspect_sin = sprintf("%s_Aspect_sin_static.tif", country_prefix),
    Aspect_cos = sprintf("%s_Aspect_cos_static.tif", country_prefix)
  )
  

  env_files <- c(env_files_dynamic, env_files_static)
  

  env_layers <- list()
  
  for (var_name in names(env_files)) {
    file_path <- file.path(env_dir, env_files[[var_name]])
    if (file.exists(file_path)) {
      r <- rast(file_path)
      names(r) <- var_name
      env_layers[[var_name]] <- r
      cat(sprintf("  : %s\n", var_name))
    } else {
      cat(sprintf("  : %s  (%s)\n", var_name, env_files[[var_name]]))
    }
  }
  
  # bio_avg
  cat("\n...\n")
  
  env_resampled <- list()
  for (var_name in names(env_layers)) {
    env_resampled[[var_name]] <- resample(env_layers[[var_name]], bio_avg, method = "bilinear")
    cat(sprintf("  : %s\n", var_name))
  }
  
  # BIO
  all_env <- c(bio_avg, rast(env_resampled))
  
  cat(sprintf("\n: %d\n", nlyr(all_env)))
  cat(sprintf(": %s\n", paste(names(all_env), collapse = ", ")))
  
  return(all_env)
}

# ============================================================================
# PART 3: VIF
# ============================================================================

calculate_vif_and_select <- function(env_stack, country_name, output_dir) {
  
  cat("\n========================================\n")
  cat(sprintf("%s VIF\n", country_name))
  cat("========================================\n")
  
  # VIF
  set.seed(42)
  n_samples <- 5000
  
  cat(sprintf(" %d VIF...\n", n_samples))
  

  sample_data <- spatSample(env_stack, size = n_samples, method = "random", 
                            na.rm = TRUE, as.df = TRUE)
  
  cat(sprintf(": %d\n", nrow(sample_data)))
  
  # NA
  sample_data <- na.omit(sample_data)
  cat(sprintf("NA: %d\n", nrow(sample_data)))
  
  # ============================================
  # 3.1 
  # ============================================
  
  cat("\n---  ---\n")
  cor_matrix <- cor(sample_data, use = "complete.obs")
  

  cor_file <- file.path(output_dir, sprintf("%s_correlation_matrix.csv", country_name))
  write.csv(cor_matrix, cor_file)
  cat(sprintf(": %s\n", cor_file))
  

  png(file.path(output_dir, sprintf("%s_correlation_heatmap.png", country_name)),
      width = 1200, height = 1200, res = 100)
  corrplot(cor_matrix, method = "color", type = "upper", 
           tl.col = "black", tl.srt = 45, tl.cex = 0.7,
           title = sprintf("%s - Correlation Matrix", country_name),
           mar = c(0,0,2,0))
  dev.off()
  cat("\n")
  
  # (|r| > 0.7)
  high_cor <- which(abs(cor_matrix) > 0.7 & upper.tri(cor_matrix), arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("\n (|r| > 0.7):\n")
    high_cor_df <- data.frame(
      Var1 = rownames(cor_matrix)[high_cor[, 1]],
      Var2 = colnames(cor_matrix)[high_cor[, 2]],
      Correlation = sapply(1:nrow(high_cor), function(i) {
        cor_matrix[high_cor[i, 1], high_cor[i, 2]]
      })
    )
    high_cor_df <- high_cor_df[order(-abs(high_cor_df$Correlation)), ]
    print(high_cor_df)
    

    write.csv(high_cor_df, 
              file.path(output_dir, sprintf("%s_high_correlation_pairs.csv", country_name)),
              row.names = FALSE)
  }
  
  # ============================================
  # 3.2 VIF
  # ============================================
  
  cat("\n--- VIF ---\n")
  
  # usdmVIF
  vif_result <- vif(sample_data)
  print(vif_result)
  
  # VIF
  vif_file <- file.path(output_dir, sprintf("%s_VIF_all_variables.csv", country_name))
  write.csv(vif_result, vif_file, row.names = FALSE)
  cat(sprintf("\nVIF: %s\n", vif_file))
  
  # ============================================
  # 3.3 VIF(=10)
  # ============================================
  
  cat("\n--- VIF (=10) ---\n")
  
  vif_selected <- vifstep(sample_data, th = 10)
  print(vif_selected)
  

  selected_vars <- vif_selected@results$Variables
  cat(sprintf("\n (%d): %s\n", length(selected_vars), paste(selected_vars, collapse = ", ")))
  

  selected_file <- file.path(output_dir, sprintf("%s_VIF_selected_variables.csv", country_name))
  write.csv(vif_selected@results, selected_file, row.names = FALSE)
  cat(sprintf(": %s\n", selected_file))
  
  # ============================================
  # 3.4 
  # ============================================
  
  report <- list(
    country = country_name,
    total_variables = ncol(sample_data),
    selected_variables = selected_vars,
    n_selected = length(selected_vars),
    vif_all = vif_result,
    vif_selected = vif_selected@results,
    correlation = cor_matrix
  )
  
  return(report)
}

# ============================================================================
# PART 4: - 
# ============================================================================

cat("\n")
cat("############################################################\n")
cat("#                     #\n")
cat("############################################################\n")

# 4.1 BIO
SEN_bio_before <- calculate_bio_average(
  bio_dir = SEN_BIO_dir,
  years = 2007:2015,
  output_dir = SEN_output_dir,
  country_prefix = "SEN",
  period_name = "before_2007_2015"
)

SEN_bio_after <- calculate_bio_average(
  bio_dir = SEN_BIO_dir,
  years = 2016:2024,
  output_dir = SEN_output_dir,
  country_prefix = "SEN",
  period_name = "after_2016_2024"
)

# 4.2 (AfterVIF)
SEN_env_after <- load_and_stack_env(
  env_dir = SEN_ENV_dir,
  country_prefix = "SEN",
  period_suffix = "after_2016_2024",
  bio_avg = SEN_bio_after
)

# 4.3 VIF
SEN_vif_report <- calculate_vif_and_select(
  env_stack = SEN_env_after,
  country_name = "Senegal",
  output_dir = SEN_output_dir
)

# ============================================================================
# PART 5: - 
# ============================================================================

cat("\n")
cat("############################################################\n")
cat("#                    #\n")
cat("############################################################\n")

# 5.1 BIO
ETH_bio_before <- calculate_bio_average(
  bio_dir = ETH_BIO_dir,
  years = 2007:2015,
  output_dir = ETH_output_dir,
  country_prefix = "ETH",
  period_name = "before_2007_2015"
)

ETH_bio_after <- calculate_bio_average(
  bio_dir = ETH_BIO_dir,
  years = 2016:2024,
  output_dir = ETH_output_dir,
  country_prefix = "ETH",
  period_name = "after_2016_2024"
)

# 5.2 (AfterVIF)
ETH_env_after <- load_and_stack_env(
  env_dir = ETH_ENV_dir,
  country_prefix = "ETH",
  period_suffix = "after_2016_2024",
  bio_avg = ETH_bio_after
)

# 5.3 VIF
ETH_vif_report <- calculate_vif_and_select(
  env_stack = ETH_env_after,
  country_name = "Ethiopia",
  output_dir = ETH_output_dir
)

# ============================================================================
# PART 6: 
# ============================================================================

cat("\n")
cat("############################################################\n")
cat("#                      #\n")
cat("############################################################\n")

cat("\n==========  ==========\n")
cat(sprintf(": %d\n", SEN_vif_report$total_variables))
cat(sprintf("VIF: %d\n", SEN_vif_report$n_selected))
cat(sprintf(": %s\n", paste(SEN_vif_report$selected_variables, collapse = ", ")))

cat("\n==========  ==========\n")
cat(sprintf(": %d\n", ETH_vif_report$total_variables))
cat(sprintf("VIF: %d\n", ETH_vif_report$n_selected))
cat(sprintf(": %s\n", paste(ETH_vif_report$selected_variables, collapse = ", ")))


cat("\n==========  ==========\n")
common_vars <- intersect(SEN_vif_report$selected_variables, ETH_vif_report$selected_variables)
sen_only <- setdiff(SEN_vif_report$selected_variables, ETH_vif_report$selected_variables)
eth_only <- setdiff(ETH_vif_report$selected_variables, SEN_vif_report$selected_variables)

cat(sprintf(" (%d): %s\n", length(common_vars), paste(common_vars, collapse = ", ")))
cat(sprintf(" (%d): %s\n", length(sen_only), paste(sen_only, collapse = ", ")))
cat(sprintf(" (%d): %s\n", length(eth_only), paste(eth_only, collapse = ", ")))

# ============================================================================
# PART 7: (VIF)
# ============================================================================

cat("\n==========  ==========\n")

# - After
SEN_env_selected <- SEN_env_after[[SEN_vif_report$selected_variables]]
writeRaster(SEN_env_selected, 
            file.path(SEN_output_dir, "SEN_ENV_selected_after_2016_2024.tif"),
            overwrite = TRUE)
cat(" After \n")

# - Before()
SEN_env_before <- load_and_stack_env(
  env_dir = SEN_ENV_dir,
  country_prefix = "SEN",
  period_suffix = "before_2007_2015",
  bio_avg = SEN_bio_before
)
SEN_env_before_selected <- SEN_env_before[[SEN_vif_report$selected_variables]]
writeRaster(SEN_env_before_selected,
            file.path(SEN_output_dir, "SEN_ENV_selected_before_2007_2015.tif"),
            overwrite = TRUE)
cat(" Before \n")

# - After
ETH_env_selected <- ETH_env_after[[ETH_vif_report$selected_variables]]
writeRaster(ETH_env_selected,
            file.path(ETH_output_dir, "ETH_ENV_selected_after_2016_2024.tif"),
            overwrite = TRUE)
cat(" After \n")

# - Before()
ETH_env_before <- load_and_stack_env(
  env_dir = ETH_ENV_dir,
  country_prefix = "ETH",
  period_suffix = "before_2007_2015",
  bio_avg = ETH_bio_before
)
ETH_env_before_selected <- ETH_env_before[[ETH_vif_report$selected_variables]]
writeRaster(ETH_env_before_selected,
            file.path(ETH_output_dir, "ETH_ENV_selected_before_2007_2015.tif"),
            overwrite = TRUE)
cat(" Before \n")

# ============================================================================
# PART 8: 
# ============================================================================

cat("\n==========  ==========\n")


summary_report <- data.frame(
  Country = c(rep("Senegal", SEN_vif_report$n_selected),
              rep("Ethiopia", ETH_vif_report$n_selected)),
  Variable = c(SEN_vif_report$selected_variables,
               ETH_vif_report$selected_variables),
  VIF = c(SEN_vif_report$vif_selected$VIF,
          ETH_vif_report$vif_selected$VIF)
)

write.csv(summary_report, 
          file.path(SEN_output_dir, "Variable_Selection_Summary.csv"),
          row.names = FALSE)
write.csv(summary_report, 
          file.path(ETH_output_dir, "Variable_Selection_Summary.csv"),
          row.names = FALSE)

cat("\n############################################################\n")
cat("#                      !                            #\n")
cat("############################################################\n")

cat("\n:\n")
cat("
:
  - SEN_BIO_avg_before_2007_2015.tif    (BIO)
  - SEN_BIO_avg_after_2016_2024.tif     (BIO)
  - Senegal_correlation_matrix.csv       ()
  - Senegal_correlation_heatmap.png      ()
  - Senegal_high_correlation_pairs.csv   ()
  - Senegal_VIF_all_variables.csv        (VIF)
  - Senegal_VIF_selected_variables.csv   (VIF)
  - SEN_ENV_selected_after_2016_2024.tif ()
  - SEN_ENV_selected_before_2007_2015.tif(hindcast)

:
  - ETH_BIO_avg_before_2007_2015.tif    (BIO)
  - ETH_BIO_avg_after_2016_2024.tif     (BIO)
  - Ethiopia_correlation_matrix.csv      ()
  - Ethiopia_correlation_heatmap.png     ()
  - Ethiopia_high_correlation_pairs.csv  ()
  - Ethiopia_VIF_all_variables.csv       (VIF)
  - Ethiopia_VIF_selected_variables.csv  (VIF)
  - ETH_ENV_selected_after_2016_2024.tif ()
  - ETH_ENV_selected_before_2007_2015.tif(hindcast)
")




2026.1.9
# ============================================================================
# - ()
# 2007-20152016-2024
# ============================================================================

library(terra)
library(dplyr)
library(tidyr)

# ============================================================================

# ============================================================================

SEN_ENV_dir <- "E:/2026.1.8_biomod2/Senegal_ENV_TwoPeriod"
ETH_ENV_dir <- "E:/2026.1.8_biomod2/Ethiopia_ENV_TwoPeriod"

# ============================================================================

# ============================================================================

compare_periods <- function(env_dir, country_prefix, country_name) {
  
  cat("\n")
  cat("############################################################\n")
  cat(sprintf("#  %s \n", country_name))
  cat("############################################################\n")
  
  # BIO
  bio_before <- rast(file.path(env_dir, sprintf("%s_BIO_avg_before_2007_2015.tif", country_prefix)))
  bio_after <- rast(file.path(env_dir, sprintf("%s_BIO_avg_after_2016_2024.tif", country_prefix)))
  
  # (before)
  dhi_mean_before <- rast(file.path(env_dir, sprintf("%s_DHI_mean_before_2007_2015.tif", country_prefix)))
  dhi_min_before <- rast(file.path(env_dir, sprintf("%s_DHI_min_before_2007_2015.tif", country_prefix)))
  dhi_cv_before <- rast(file.path(env_dir, sprintf("%s_DHI_cv_before_2007_2015.tif", country_prefix)))
  shdi_before <- rast(file.path(env_dir, sprintf("%s_SHDI_before_2007_2015.tif", country_prefix)))
  woody_before <- rast(file.path(env_dir, sprintf("%s_Woody_pct_before_2007_2015.tif", country_prefix)))
  herb_before <- rast(file.path(env_dir, sprintf("%s_Herbaceous_pct_before_2007_2015.tif", country_prefix)))
  
  # (after)
  dhi_mean_after <- rast(file.path(env_dir, sprintf("%s_DHI_mean_after_2016_2024.tif", country_prefix)))
  dhi_min_after <- rast(file.path(env_dir, sprintf("%s_DHI_min_after_2016_2024.tif", country_prefix)))
  dhi_cv_after <- rast(file.path(env_dir, sprintf("%s_DHI_cv_after_2016_2024.tif", country_prefix)))
  shdi_after <- rast(file.path(env_dir, sprintf("%s_SHDI_after_2016_2024.tif", country_prefix)))
  woody_after <- rast(file.path(env_dir, sprintf("%s_Woody_pct_after_2016_2023.tif", country_prefix)))
  herb_after <- rast(file.path(env_dir, sprintf("%s_Herbaceous_pct_after_2016_2023.tif", country_prefix)))
  

  names(dhi_mean_before) <- "DHI_mean"; names(dhi_mean_after) <- "DHI_mean"
  names(dhi_min_before) <- "DHI_min"; names(dhi_min_after) <- "DHI_min"
  names(dhi_cv_before) <- "DHI_cv"; names(dhi_cv_after) <- "DHI_cv"
  names(shdi_before) <- "SHDI"; names(shdi_after) <- "SHDI"
  names(woody_before) <- "Woody_pct"; names(woody_after) <- "Woody_pct"
  names(herb_before) <- "Herbaceous_pct"; names(herb_after) <- "Herbaceous_pct"
  
  # ============================================

  # ============================================
  
  set.seed(42)
  n_samples <- 5000
  
  # (projectCRS)
  cat("\n...\n")
  
  template <- bio_before[[1]]
  
  dhi_mean_before <- project(dhi_mean_before, template, method = "bilinear")
  dhi_min_before <- project(dhi_min_before, template, method = "bilinear")
  dhi_cv_before <- project(dhi_cv_before, template, method = "bilinear")
  shdi_before <- project(shdi_before, template, method = "bilinear")
  woody_before <- project(woody_before, template, method = "bilinear")
  herb_before <- project(herb_before, template, method = "bilinear")
  
  template_after <- bio_after[[1]]
  
  dhi_mean_after <- project(dhi_mean_after, template_after, method = "bilinear")
  dhi_min_after <- project(dhi_min_after, template_after, method = "bilinear")
  dhi_cv_after <- project(dhi_cv_after, template_after, method = "bilinear")
  shdi_after <- project(shdi_after, template_after, method = "bilinear")
  woody_after <- project(woody_after, template_after, method = "bilinear")
  herb_after <- project(herb_after, template_after, method = "bilinear")
  

  env_before <- c(bio_before, dhi_mean_before, dhi_min_before, dhi_cv_before, 
                  shdi_before, woody_before, herb_before)
  env_after <- c(bio_after, dhi_mean_after, dhi_min_after, dhi_cv_after,
                 shdi_after, woody_after, herb_after)
  
  cat(sprintf("Before: %d\n", nlyr(env_before)))
  cat(sprintf("After: %d\n", nlyr(env_after)))
  

  cat("...\n")
  sample_before <- spatSample(env_before, size = n_samples, method = "random", na.rm = TRUE, as.df = TRUE)
  sample_after <- spatSample(env_after, size = n_samples, method = "random", na.rm = TRUE, as.df = TRUE)
  
  cat(sprintf("Before: %d\n", nrow(sample_before)))
  cat(sprintf("After: %d\n", nrow(sample_after)))
  
  # ============================================

  # ============================================
  
  var_names <- names(env_before)
  
  results <- data.frame(
    Variable = var_names,
    Before_Mean = NA_real_,
    Before_SD = NA_real_,
    Before_Min = NA_real_,
    Before_Max = NA_real_,
    Before_P5 = NA_real_,
    Before_P95 = NA_real_,
    After_Mean = NA_real_,
    After_SD = NA_real_,
    After_Min = NA_real_,
    After_Max = NA_real_,
    After_P5 = NA_real_,
    After_P95 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(var_names)) {
    var_name <- var_names[i]
    
    before_vals <- sample_before[[var_name]]
    after_vals <- sample_after[[var_name]]
    
    results$Before_Mean[i] <- mean(before_vals, na.rm = TRUE)
    results$Before_SD[i] <- sd(before_vals, na.rm = TRUE)
    results$Before_Min[i] <- min(before_vals, na.rm = TRUE)
    results$Before_Max[i] <- max(before_vals, na.rm = TRUE)
    results$Before_P5[i] <- quantile(before_vals, 0.05, na.rm = TRUE)
    results$Before_P95[i] <- quantile(before_vals, 0.95, na.rm = TRUE)
    
    results$After_Mean[i] <- mean(after_vals, na.rm = TRUE)
    results$After_SD[i] <- sd(after_vals, na.rm = TRUE)
    results$After_Min[i] <- min(after_vals, na.rm = TRUE)
    results$After_Max[i] <- max(after_vals, na.rm = TRUE)
    results$After_P5[i] <- quantile(after_vals, 0.05, na.rm = TRUE)
    results$After_P95[i] <- quantile(after_vals, 0.95, na.rm = TRUE)
  }
  
  # ============================================

  # ============================================
  
  results$Mean_Change_Pct <- (results$After_Mean - results$Before_Mean) / abs(results$Before_Mean) * 100
  
  results$Range_Overlap <- pmax(0, pmin(results$After_Max, results$Before_Max) - 
                                  pmax(results$After_Min, results$Before_Min)) / 
    (pmax(results$After_Max, results$Before_Max) - pmin(results$After_Min, results$Before_Min))
  
  results$Standardized_Change <- abs(results$After_Mean - results$Before_Mean) / results$Before_SD
  
  # SD0NA
  results$Standardized_Change[is.na(results$Standardized_Change) | is.infinite(results$Standardized_Change)] <- 0
  

  results$Extrapolation_Risk <- ifelse(
    results$After_P5 < results$Before_P5 | results$After_P95 > results$Before_P95,
    "HIGH",
    ifelse(abs(results$Mean_Change_Pct) > 20, "MEDIUM", "LOW")
  )
  
  # ============================================

  # ============================================
  
  cat("\n==========  ==========\n\n")
  

  results_sorted <- results[order(-results$Standardized_Change), ]
  

  print_df <- data.frame(
    Variable = results_sorted$Variable,
    Before_Mean = round(results_sorted$Before_Mean, 2),
    After_Mean = round(results_sorted$After_Mean, 2),
    Change_Pct = round(results_sorted$Mean_Change_Pct, 1),
    Std_Change = round(results_sorted$Standardized_Change, 3),
    Overlap = round(results_sorted$Range_Overlap, 3),
    Risk = results_sorted$Extrapolation_Risk
  )
  
  print(print_df)
  
  # ============================================

  # ============================================
  
  cat("\n==========  ==========\n")
  
  high_risk <- results[results$Extrapolation_Risk == "HIGH" | results$Standardized_Change > 0.5, ]
  high_risk <- high_risk[order(-high_risk$Standardized_Change), ]
  
  if (nrow(high_risk) > 0) {
    cat("\n[!]  ():\n")
    for (i in 1:nrow(high_risk)) {
      cat(sprintf("  - %s:  %.1f%%,  %.3f,  %.3f\n",
                  high_risk$Variable[i], 
                  high_risk$Mean_Change_Pct[i],
                  high_risk$Standardized_Change[i],
                  high_risk$Range_Overlap[i]))
    }
  } else {
    cat("\n[OK] \n")
  }
  
  medium_risk <- results[results$Extrapolation_Risk == "MEDIUM" & results$Standardized_Change <= 0.5, ]
  
  if (nrow(medium_risk) > 0) {
    cat("\n[*]  ():\n")
    for (i in 1:nrow(medium_risk)) {
      cat(sprintf("  - %s:  %.1f%%\n",
                  medium_risk$Variable[i],
                  medium_risk$Mean_Change_Pct[i]))
    }
  }
  
  low_risk <- results[results$Extrapolation_Risk == "LOW", ]
  
  cat(sprintf("\n[OK]  (%d): %s\n", 
              nrow(low_risk), 
              paste(low_risk$Variable, collapse = ", ")))
  

  output_file <- file.path(env_dir, sprintf("%s_Period_Comparison.csv", country_name))
  write.csv(results, output_file, row.names = FALSE)
  cat(sprintf("\n: %s\n", output_file))
  
  return(results)
}

# ============================================================================

# ============================================================================


cat("\n...\n")
SEN_comparison <- compare_periods(SEN_ENV_dir, "SEN", "Senegal")


cat("\n...\n")
ETH_comparison <- compare_periods(ETH_ENV_dir, "ETH", "Ethiopia")

# ============================================================================

# ============================================================================

cat("\n")
cat("############################################################\n")
cat("#                    #\n")
cat("############################################################\n")

cat("\n\n")
cat("\n1.  (Standardized_Change < 0.3)\n")
cat("2. Range_Overlap < 0.7\n")
cat("3. \n")
cat("   - \n")
cat("   - \n")
cat("   - \n")

cat("\n############################################################\n")
cat("#                      #\n")
cat("############################################################\n")