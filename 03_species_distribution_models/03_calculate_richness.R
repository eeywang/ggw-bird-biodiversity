#!/usr/bin/env Rscript
# =============================================================================
# 03_calculate_richness.R — Species richness for S0 (Before) and S3 (After)
#
# Stacks ensemble predictions across all successfully modelled species.
# No additional occurrence filter is needed here because the ≥50 threshold
# was already applied at the modelling stage (02_biomod2_modeling.R).
#
# Outputs: binary, probability-weighted, and TSS-weighted richness rasters
#          for both periods, plus change maps.
#
# Usage: Rscript 03_calculate_richness.R <country>
# =============================================================================

source("00_hpc_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
})

richness_dir <- file.path(hpc_paths$results, "richness")
dir.create(richness_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Find all species with successful ensemble projections
# =============================================================================
cat("\n=== Collecting species with ensemble projections ===\n")

all_dirs <- list.dirs(hpc_paths$work_dir, recursive = FALSE, full.names = TRUE)
species_dirs <- all_dirs[grepl("/[A-Z][a-z]+\\.[a-z]", all_dirs)]

after_bin_f <- after_prob_f <- before_bin_f <- before_prob_f <- c()
species_tss <- c()
included_species <- c()

for (sp_dir in species_dirs) {
  sp <- basename(sp_dir)

  ab <- file.path(sp_dir, "proj_After_2016_2024",
                  paste0("proj_After_2016_2024_", sp, "_ensemble_TSSbin.tif"))
  ap <- file.path(sp_dir, "proj_After_2016_2024",
                  paste0("proj_After_2016_2024_", sp, "_ensemble.tif"))
  bb <- file.path(sp_dir, "proj_Before_2007_2015",
                  paste0("proj_Before_2007_2015_", sp, "_ensemble_TSSbin.tif"))
  bp <- file.path(sp_dir, "proj_Before_2007_2015",
                  paste0("proj_Before_2007_2015_", sp, "_ensemble.tif"))

  if (all(file.exists(c(ab, ap, bb, bp)))) {
    after_bin_f  <- c(after_bin_f, ab);   after_prob_f <- c(after_prob_f, ap)
    before_bin_f <- c(before_bin_f, bb);  before_prob_f <- c(before_prob_f, bp)
    included_species <- c(included_species, sp)

    # TSS weight
    ef <- list.files(sp_dir, pattern = "evaluation.csv",
                     recursive = TRUE, full.names = TRUE)
    tv <- if (length(ef) > 0) {
      ev <- read.csv(ef[1])
      mean(ev$validation[ev$metric.eval == "TSS"], na.rm = TRUE)
    } else 0.5
    species_tss <- c(species_tss, ifelse(is.na(tv), 0.5, tv))
  }
}

cat("Species included:", length(included_species),
    "(all modelled with >=", THRESHOLDS$min_occurrences, "records)\n")
write.csv(data.frame(species = included_species),
          file.path(richness_dir, "included_species.csv"), row.names = FALSE)

# =============================================================================
# Compute richness
# =============================================================================
cat("\n=== Computing richness ===\n")

template <- rast(after_bin_f[1])[[1]]
r_bin_after <- r_bin_before <- template * 0
r_prob_after <- r_prob_before <- template * 0
r_tss_after <- r_tss_before <- template * 0

for (i in seq_along(after_bin_f)) {
  tryCatch({
    a_bin <- ifel(rast(after_bin_f[i])[[1]] > 0, 1, 0)
    b_bin <- ifel(rast(before_bin_f[i])[[1]] > 0, 1, 0)
    a_prob <- rast(after_prob_f[i])[[1]] / 1000
    b_prob <- rast(before_prob_f[i])[[1]] / 1000

    a_bin[is.na(a_bin)] <- 0;  b_bin[is.na(b_bin)] <- 0
    a_prob[is.na(a_prob)] <- 0; b_prob[is.na(b_prob)] <- 0

    r_bin_after  <- r_bin_after + a_bin
    r_bin_before <- r_bin_before + b_bin
    r_prob_after  <- r_prob_after + a_prob
    r_prob_before <- r_prob_before + b_prob
    r_tss_after  <- r_tss_after + (a_prob * species_tss[i])
    r_tss_before <- r_tss_before + (b_prob * species_tss[i])

    if (i %% 50 == 0) cat("  ", i, "/", length(after_bin_f), "\n")
  }, error = function(e) cat("  Error at", i, ":", conditionMessage(e), "\n"))
}

# =============================================================================
# Save
# =============================================================================
cat("\n=== Saving ===\n")

for (nm in c("binary", "prob", "tss")) {
  for (pd in c("after", "before")) {
    writeRaster(get(paste0("r_", nm, "_", pd)),
      file.path(richness_dir, paste0("richness_", nm, "_", pd, ".tif")),
      overwrite = TRUE)
  }
  writeRaster(get(paste0("r_", nm, "_after")) - get(paste0("r_", nm, "_before")),
    file.path(richness_dir, paste0("richness_", nm, "_change.tif")),
    overwrite = TRUE)
}

cat("\n=== Done for", COUNTRY, "===\n")
cat("Species:", length(included_species), "| Output:", richness_dir, "\n")
