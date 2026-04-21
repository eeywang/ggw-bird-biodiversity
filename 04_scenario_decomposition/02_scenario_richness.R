#!/usr/bin/env Rscript
# =============================================================================
# 02_scenario_richness.R — S1/S2 richness + 2×2 effect decomposition
#
# All successfully modelled species are included (already filtered to ≥50
# at the modelling stage).
#
# Usage: Rscript 02_scenario_richness.R <country>
# =============================================================================

source("../03_species_distribution_models/00_hpc_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
})

out_dir <- file.path(hpc_paths$results, "scenario_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load S0/S3 from richness directory
# =============================================================================
rich_dir <- file.path(hpc_paths$results, "richness")
S0_bin  <- rast(file.path(rich_dir, "richness_binary_before.tif"))
S3_bin  <- rast(file.path(rich_dir, "richness_binary_after.tif"))
S0_prob <- rast(file.path(rich_dir, "richness_prob_before.tif"))
S3_prob <- rast(file.path(rich_dir, "richness_prob_after.tif"))
S0_tss  <- rast(file.path(rich_dir, "richness_tss_before.tif"))
S3_tss  <- rast(file.path(rich_dir, "richness_tss_after.tif"))

cat("S0 range:", round(minmax(S0_bin)[1]), "-", round(minmax(S0_bin)[2]), "\n")
cat("S3 range:", round(minmax(S3_bin)[1]), "-", round(minmax(S3_bin)[2]), "\n")

# =============================================================================
# 2. Collect S1/S2 ensemble projections
# =============================================================================
cat("\n=== Computing S1/S2 richness ===\n")

all_dirs <- list.dirs(hpc_paths$work_dir, recursive = FALSE, full.names = TRUE)
species_dirs <- all_dirs[grepl("/[A-Z][a-z]+\\.[a-z]", all_dirs)]

get_species_tss <- function(sp_dir) {
  ef <- list.files(sp_dir, pattern = "evaluation.csv",
                   recursive = TRUE, full.names = TRUE)
  if (length(ef) > 0) {
    ev <- read.csv(ef[1])
    tv <- mean(ev$validation[ev$metric.eval == "TSS"], na.rm = TRUE)
    if (!is.na(tv)) return(tv)
  }
  0.5
}

get_ensemble_bin <- function(sp_dir, proj_name) {
  pdir <- file.path(sp_dir, paste0("proj_", proj_name))
  if (!dir.exists(pdir)) return(NULL)
  f <- list.files(pdir, pattern = "ensemble_TSSbin\\.tif$", full.names = TRUE)
  if (length(f) == 0) return(NULL)
  tryCatch({ r <- rast(f[1]); if (nlyr(r) > 1) r <- r[[1]]; r },
           error = function(e) NULL)
}

get_ensemble_prob <- function(sp_dir, proj_name) {
  pdir <- file.path(sp_dir, paste0("proj_", proj_name))
  if (!dir.exists(pdir)) return(NULL)
  f <- list.files(pdir, pattern = "ensemble\\.tif$", full.names = TRUE)
  f <- f[!grepl("TSSbin|ROCbin", f)]
  if (length(f) == 0) return(NULL)
  tryCatch({ r <- rast(f[1]); if (nlyr(r) > 1) r <- r[[1]]; r / 1000 },
           error = function(e) NULL)
}

template <- S0_bin
S1 <- list(binary = template * 0, prob = template * 0, tss = template * 0)
S2 <- list(binary = template * 0, prob = template * 0, tss = template * 0)
cnt_s1 <- cnt_s2 <- 0

for (i in seq_along(species_dirs)) {
  sp_dir <- species_dirs[i]
  tv <- get_species_tss(sp_dir)

  for (sc in c("S1_climate_only", "S2_vegetation_only")) {
    b <- get_ensemble_bin(sp_dir, sc)
    p <- get_ensemble_prob(sp_dir, sc)
    if (!is.null(b) && !is.null(p)) {
      b[is.na(b)] <- 0; p[is.na(p)] <- 0
      tgt <- if (grepl("S1", sc)) "S1" else "S2"
      S <- get(tgt)
      S$binary <- S$binary + ifel(b > 0, 1, 0)
      S$prob   <- S$prob + p
      S$tss    <- S$tss + p * tv
      assign(tgt, S)
      if (tgt == "S1") cnt_s1 <- cnt_s1 + 1 else cnt_s2 <- cnt_s2 + 1
    }
  }
  if (i %% 50 == 0) cat("  ", i, "/", length(species_dirs), "\n")
}
cat("S1:", cnt_s1, "species | S2:", cnt_s2, "species\n")

for (nm in c("binary", "prob", "tss")) {
  writeRaster(S1[[nm]], file.path(out_dir, paste0("richness_S1_", nm, ".tif")), overwrite = TRUE)
  writeRaster(S2[[nm]], file.path(out_dir, paste0("richness_S2_", nm, ".tif")), overwrite = TRUE)
}

# =============================================================================
# 3. Effect decomposition
# =============================================================================
cat("\n=== Effect decomposition ===\n")

decompose <- function(s0, s1, s2, s3, type) {
  total <- s3 - s0; clim <- s1 - s0; veg <- s2 - s0
  inter <- total - clim - veg

  stats <- data.frame(
    Type = type, Country = COUNTRY,
    Effect = c("Total", "Climate", "Vegetation", "Interaction"),
    Mean = sapply(list(total, clim, veg, inter),
                  function(r) global(r, "mean", na.rm = TRUE)[[1]]),
    SD = sapply(list(total, clim, veg, inter),
                function(r) global(r, "sd", na.rm = TRUE)[[1]])
  )

  for (e in list(c("total", total), c("climate", clim),
                 c("vegetation", veg), c("interaction", inter)))
    writeRaster(e[[2]], file.path(out_dir,
      paste0("effect_", tolower(type), "_", e[[1]], ".tif")), overwrite = TRUE)

  cat("  ", type, ": Total=", round(stats$Mean[1], 2),
      "Clim=", round(stats$Mean[2], 2),
      "Veg=", round(stats$Mean[3], 2),
      "Inter=", round(stats$Mean[4], 2), "\n")
  stats
}

all_stats <- rbind(
  decompose(S0_bin,  S1$binary, S2$binary, S3_bin,  "Binary"),
  decompose(S0_prob, S1$prob,   S2$prob,   S3_prob, "Probability"),
  decompose(S0_tss,  S1$tss,    S2$tss,    S3_tss,  "TSS_weighted")
)
write.csv(all_stats, file.path(out_dir, "effect_summary.csv"), row.names = FALSE)

cat("\n=== Done for", COUNTRY, "===\n")
print(all_stats)
