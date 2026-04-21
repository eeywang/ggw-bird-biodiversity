#!/usr/bin/env Rscript
# =============================================================================
# 02_biomod2_modeling.R — Ensemble SDM training + current & hindcast projections
#
# Trains 5 algorithms (GLM, GAM, RF, MAXNET, GBM) per species, builds
# ensemble models, and projects to both After (2016-2024) and Before (2007-2015).
#
# Modelling threshold: ≥50 occurrence records (consistent across all 3 countries).
#
# Usage (SLURM array):  Rscript 02_biomod2_modeling.R <species_index> <country>
# Usage (sequential):   Rscript 02_biomod2_modeling.R <country>
# =============================================================================

source("00_hpc_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
  library(dplyr)
  library(mgcv)
  library(randomForest)
  library(maxnet)
  library(gbm)
})

cat("biomod2 version:", as.character(packageVersion("biomod2")), "\n")

# =============================================================================
# Load data
# =============================================================================

env_stack  <- rast(hpc_paths$env_after)
env_before <- rast(hpc_paths$env_before)
cat("After env:", nlyr(env_stack), "| Before env:", nlyr(env_before), "\n")

if (nlyr(env_stack) != THRESHOLDS$n_env_vars)
  stop("Expected ", THRESHOLDS$n_env_vars, " env variables, got ", nlyr(env_stack))

occ_data_raw <- read.csv(hpc_paths$occurrence)

# Filter to training period (2016-2024)
occ_data <- occ_data_raw %>% filter(year >= 2016, year <= 2024)
cat("Training period records:", format(nrow(occ_data), big.mark = ","), "\n")

# Build species list: ≥30 records
sp_counts <- occ_data %>%
  group_by(species) %>% summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
species_list <- sp_counts %>%
  filter(n >= THRESHOLDS$min_occurrences) %>% pull(species)
cat("Species for modelling (>=", THRESHOLDS$min_occurrences, "):",
    length(species_list), "\n")

# =============================================================================
# Single-species function
# =============================================================================

run_single_species <- function(sp_name) {

  sp_name_clean <- gsub(" ", "_", sp_name)
  cat("\n==================================================\n")
  cat("Species:", sp_name, "| Country:", COUNTRY, "\n")
  cat("==================================================\n")

  sp_dir <- file.path(hpc_paths$models, sp_name_clean)
  dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    # 1. Occurrence points
    sp_occ <- occ_data %>%
      filter(species == sp_name) %>%
      dplyr::select(decimalLongitude, decimalLatitude) %>% distinct()

    occ_env <- extract(env_stack, as.matrix(sp_occ))
    sp_occ <- sp_occ[complete.cases(occ_env), ]
    cat("  Valid points:", nrow(sp_occ), "\n")

    if (nrow(sp_occ) < THRESHOLDS$min_valid_points) {
      cat("  SKIP: <", THRESHOLDS$min_valid_points, "valid points\n")
      return(NULL)
    }

    # 2. Format data
    myData <- BIOMOD_FormatingData(
      resp.var = rep(1, nrow(sp_occ)), resp.xy = as.matrix(sp_occ),
      resp.name = sp_name_clean, expl.var = env_stack,
      PA.nb.rep = THRESHOLDS$n_pa_sets,
      PA.nb.absences = min(nrow(sp_occ) * 5, 5000),
      PA.strategy = "random", filter.raster = TRUE
    )

    # 3. Train models
    start_time <- Sys.time()
    myModel <- BIOMOD_Modeling(
      bm.format = myData, modeling.id = "prod", models = SELECTED_MODELS,
      CV.strategy = "kfold", CV.nb.rep = THRESHOLDS$cv_runs,
      CV.k = THRESHOLDS$cv_folds, var.import = 3,
      metric.eval = c("TSS", "ROC"), seed.val = 42, do.progress = FALSE
    )
    model_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

    # 4. Evaluate — tiered model selection
    eval_df <- get_evaluations(myModel)
    mean_tss <- mean(eval_df$validation[eval_df$metric.eval == "TSS"], na.rm = TRUE)
    mean_auc <- mean(eval_df$validation[eval_df$metric.eval == "ROC"], na.rm = TRUE)
    cat("  TSS:", round(mean_tss, 3), "| AUC:", round(mean_auc, 3), "\n")

    good_models <- eval_df %>%
      filter(metric.eval == "TSS", validation >= THRESHOLDS$tss_min) %>%
      pull(full.name)
    if (length(good_models) == 0) {
      good_models <- eval_df %>%
        filter(metric.eval == "TSS") %>%
        arrange(desc(validation)) %>% head(5) %>% pull(full.name)
    }
    cat("  Selected models:", length(good_models), "\n")

    # 5. Project: After (current) + Before (hindcast)
    cat("  Projecting After_2016_2024...\n")
    proj_after <- BIOMOD_Projection(
      bm.mod = myModel, proj.name = "After_2016_2024",
      new.env = env_stack, models.chosen = good_models,
      metric.binary = "TSS", compress = TRUE, output.format = ".tif"
    )
    cat("  Projecting Before_2007_2015 (hindcast)...\n")
    proj_before <- BIOMOD_Projection(
      bm.mod = myModel, proj.name = "Before_2007_2015",
      new.env = env_before, models.chosen = good_models,
      metric.binary = "TSS", compress = TRUE, output.format = ".tif"
    )

    # 6. Ensemble + ensemble projections
    if (length(good_models) >= 2) {
      cat("  Ensemble modelling...\n")
      myEnsemble <- BIOMOD_EnsembleModeling(
        bm.mod = myModel, models.chosen = good_models, em.by = "all",
        em.algo = c("EMwmean", "EMca"), metric.select = "TSS",
        metric.select.thresh = THRESHOLDS$ensemble_thresh,
        metric.eval = c("TSS", "ROC"), var.import = 3
      )
      BIOMOD_EnsembleForecasting(bm.em = myEnsemble, bm.proj = proj_after,
                                  metric.binary = "TSS", output.format = ".tif")
      BIOMOD_EnsembleForecasting(bm.em = myEnsemble, bm.proj = proj_before,
                                  metric.binary = "TSS", output.format = ".tif")
    }

    # 7. Save
    write.csv(eval_df, file.path(sp_dir, "evaluation.csv"), row.names = FALSE)
    result_row <- data.frame(
      species = sp_name, country = COUNTRY, n_occ = nrow(sp_occ),
      mean_TSS = mean_tss, mean_AUC = mean_auc,
      n_good_models = length(good_models), time_min = model_time,
      status = "success"
    )
    idx <- if (!is.null(SPECIES_INDEX)) SPECIES_INDEX else
           which(species_list == sp_name)
    write.csv(result_row,
              file.path(hpc_paths$results, paste0("result_", idx, ".csv")),
              row.names = FALSE)
    cat("  Done (", round(model_time, 1), "min)\n")
    return(result_row)

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    idx <- if (!is.null(SPECIES_INDEX)) SPECIES_INDEX else
           which(species_list == sp_name)
    write.csv(data.frame(species = sp_name, country = COUNTRY,
                         status = "error", error = conditionMessage(e)),
              file.path(hpc_paths$results, paste0("result_", idx, ".csv")),
              row.names = FALSE)
    return(NULL)
  })
}

# =============================================================================
# Execute
# =============================================================================

if (!is.null(SPECIES_INDEX)) {
  if (SPECIES_INDEX > length(species_list)) {
    cat("Index", SPECIES_INDEX, ">", length(species_list), "species. Exiting.\n")
    quit(save = "no", status = 0)
  }
  run_single_species(species_list[SPECIES_INDEX])
} else {
  cat("Sequential mode: all", length(species_list), "species\n")
  for (i in seq_along(species_list)) {
    SPECIES_INDEX <<- i
    cat("\nProgress:", i, "/", length(species_list), "\n")
    run_single_species(species_list[i])
  }
}
cat("\n=== Modelling complete for", COUNTRY, "===\n")
