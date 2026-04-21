#!/usr/bin/env Rscript
# =============================================================================
# 00_hpc_config.R — Shared HPC configuration for 3-country biomod2 pipeline
#
# Sourced by all scripts in 03_species_distribution_models/ and
# 04_scenario_decomposition/. Sets country, paths, and thresholds.
#
# Usage:  Rscript 02_biomod2_modeling.R  <species_index>  <country>
#         e.g.   Rscript 02_biomod2_modeling.R  1  Nigeria
#
# The COUNTRY argument is always the LAST command-line argument.
# On HPC (Aire), pass via SLURM: --export=COUNTRY=Nigeria
# =============================================================================

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")

# ---- Parse country from command line or environment ----
args <- commandArgs(trailingOnly = TRUE)

# Country: last argument, or env var, or default
COUNTRY <- Sys.getenv("COUNTRY", unset = "")
if (COUNTRY == "" && length(args) >= 1) {
  # If last arg looks like a country name (starts with capital, not a number)
  last_arg <- args[length(args)]
  if (grepl("^[A-Z]", last_arg) && !grepl("^[0-9]", last_arg)) {
    COUNTRY <- last_arg
    args <- args[-length(args)]  # remove country from args
  }
}
if (COUNTRY == "") COUNTRY <- "Nigeria"  # default for testing

# Species index: first numeric argument (for SLURM array jobs)
SPECIES_INDEX <- NULL
if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) {
  SPECIES_INDEX <- as.integer(args[1])
}

cat("========================================\n")
cat("  Country:", COUNTRY, "\n")
if (!is.null(SPECIES_INDEX)) cat("  Species index:", SPECIES_INDEX, "\n")
cat("========================================\n\n")

# ---- Validate country ----
valid_countries <- c("Nigeria", "Senegal", "Ethiopia")
if (!COUNTRY %in% valid_countries) {
  stop("Invalid country: ", COUNTRY,
       ". Must be one of: ", paste(valid_countries, collapse = ", "))
}

# ---- Paths ----
SCRATCH <- Sys.getenv("SCRATCH", unset = "/mnt/scratch/eeywang")
work_dir <- file.path(SCRATCH, "ggw_birds", COUNTRY)

if (!dir.exists(work_dir)) {
  # Fallback: try flat structure (ggw_birds_nigeria etc.)
  work_dir_alt <- file.path(SCRATCH, paste0("ggw_birds_", tolower(COUNTRY)))
  if (dir.exists(work_dir_alt)) {
    work_dir <- work_dir_alt
  } else {
    stop("Working directory not found: ", work_dir, " or ", work_dir_alt)
  }
}

setwd(work_dir)
cat("Working directory:", work_dir, "\n")

hpc_paths <- list(
  work_dir     = work_dir,
  env_after    = file.path(work_dir, "data", "env_stack_after_FINAL.tif"),
  env_before   = file.path(work_dir, "data", "env_stack_before_FINAL.tif"),
  occurrence   = file.path(work_dir, "data", "occurrence_data.csv"),
  species_list = file.path(work_dir, "data", "species_for_modeling.csv"),
  guild_file   = file.path(work_dir, "data",
                           paste0(COUNTRY, "_validated_functional_groups.xlsx")),
  scenario_s1  = file.path(work_dir, "data", "scenario_stacks",
                           "S1_climate_only_climate_P2_veg_P1.tif"),
  scenario_s2  = file.path(work_dir, "data", "scenario_stacks",
                           "S2_vegetation_only_climate_P1_veg_P2.tif"),
  models       = file.path(work_dir, "models"),
  results      = file.path(work_dir, "results"),
  logs         = file.path(work_dir, "logs")
)

for (p in c(hpc_paths$models, hpc_paths$results, hpc_paths$logs)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# KEY THRESHOLDS — clearly documented for manuscript reproducibility
# =============================================================================
THRESHOLDS <- list(
  # ---- Occurrence threshold ----
  # Only species with ≥50 occurrence records are modelled.
  # This ensures adequate sample size for 5-algorithm ensemble SDMs
  # with cross-validation.
  min_occurrences    = 50,   # Minimum records to model a species
  min_valid_points   = 20,   # Minimum valid points after env-NA removal

  # ---- Model selection ----
  tss_min            = 0.4,  # Minimum TSS for model selection (else top-5)
  ensemble_thresh    = 0.4,  # TSS threshold for ensemble model inclusion

  # ---- Cross-validation ----
  n_pa_sets          = 2,    # Number of pseudo-absence sets
  cv_folds           = 3,    # Cross-validation folds
  cv_runs            = 2,    # Cross-validation repetitions

  # ---- Environment ----
  n_env_vars         = 12    # Expected number of environmental variables
)

cat("\nThresholds:\n")
cat("  Min occurrences: >=", THRESHOLDS$min_occurrences, "records\n")
cat("  TSS min:", THRESHOLDS$tss_min, "\n")

# ---- Selected algorithms ----
SELECTED_MODELS <- c("GLM", "GAM", "RF", "MAXNET", "GBM")

cat("  Algorithms:", paste(SELECTED_MODELS, collapse = ", "), "\n\n")
