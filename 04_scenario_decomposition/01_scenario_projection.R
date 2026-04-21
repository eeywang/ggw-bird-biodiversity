#!/usr/bin/env Rscript
# =============================================================================
# 01_scenario_projection.R — Project SDMs onto counterfactual scenarios S1 & S2
#
# S1 = climate-only change (P2 climate + P1 vegetation)
# S2 = vegetation-only change (P1 climate + P2 vegetation)
#
# Usage: Rscript 01_scenario_projection.R <species_index> <country>
# =============================================================================

source("../03_species_distribution_models/00_hpc_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
})

env_s1 <- rast(hpc_paths$scenario_s1)
env_s2 <- rast(hpc_paths$scenario_s2)
cat("S1 variables:", nlyr(env_s1), "| S2 variables:", nlyr(env_s2), "\n")

# Find species with trained models
species_dirs <- list.dirs(hpc_paths$work_dir, recursive = FALSE, full.names = FALSE)
species_dirs <- species_dirs[grepl("^[A-Z]", species_dirs)]
species_dirs <- species_dirs[sapply(species_dirs, function(d) {
  file.exists(file.path(hpc_paths$work_dir, d, paste0(d, ".v2.models.out")))
})]
species_dirs <- sort(species_dirs)
cat("Species with models:", length(species_dirs), "\n")

project_species <- function(idx) {
  if (idx > length(species_dirs)) {
    cat("Index exceeds species count, skipping\n"); return(NULL)
  }

  sp_name <- species_dirs[idx]
  cat(sprintf("\n[%d/%d] %s\n", idx, length(species_dirs), sp_name))

  model_file    <- file.path(hpc_paths$work_dir, sp_name, paste0(sp_name, ".v2.models.out"))
  ensemble_file <- file.path(hpc_paths$work_dir, sp_name, paste0(sp_name, ".v2.ensemble.models.out"))

  if (!file.exists(model_file) || !file.exists(ensemble_file)) {
    cat("  Skip: model files missing\n"); return(NULL)
  }

  tryCatch({
    myModel    <- get(load(model_file))
    myEnsemble <- get(load(ensemble_file))

    # S1 projection
    cat("  S1_climate_only...\n")
    proj_s1 <- BIOMOD_Projection(
      bm.mod = myModel, proj.name = "S1_climate_only", new.env = env_s1,
      metric.binary = "TSS", build.clamping.mask = TRUE, output.format = ".tif"
    )
    BIOMOD_EnsembleForecasting(
      bm.em = myEnsemble, bm.proj = proj_s1,
      metric.binary = "TSS", output.format = ".tif"
    )

    # S2 projection
    cat("  S2_vegetation_only...\n")
    proj_s2 <- BIOMOD_Projection(
      bm.mod = myModel, proj.name = "S2_vegetation_only", new.env = env_s2,
      metric.binary = "TSS", build.clamping.mask = TRUE, output.format = ".tif"
    )
    BIOMOD_EnsembleForecasting(
      bm.em = myEnsemble, bm.proj = proj_s2,
      metric.binary = "TSS", output.format = ".tif"
    )

    cat("  Done\n")
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
}

# Execute
if (!is.null(SPECIES_INDEX)) {
  project_species(SPECIES_INDEX)
} else {
  for (i in seq_along(species_dirs)) project_species(i)
}
cat("\n=== Scenario projection complete for", COUNTRY, "===\n")
