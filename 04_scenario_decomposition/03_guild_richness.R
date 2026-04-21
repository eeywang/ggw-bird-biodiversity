#!/usr/bin/env Rscript
# =============================================================================
# 03_guild_richness.R — Calculate richness by habitat guild
#                        (Woodland, Open-habitat, Wetland) for S0 and S3
#
# Usage: Rscript 03_guild_richness.R Nigeria
# =============================================================================

source("../03_species_distribution_models/00_hpc_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(readxl)
})

output_dir <- file.path(hpc_paths$results, "guild_richness")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

periods <- c("S0" = "Before_2007_2015", "S3" = "Current")
guilds  <- c("Woodland", "Open-habitat", "Wetland")

# ---- Load species classification ----
cat("Loading guild classification for", COUNTRY, "...\n")
guild_data <- read_excel(hpc_paths$guild_file, sheet = "Species Classification")
colnames(guild_data) <- c("No", "Species", "Habitat_Guild",
                           "Trophic_Niche", "Migratory_Status")
guild_data$species_folder <- gsub(" ", ".", guild_data$Species)
cat("  Total species:", nrow(guild_data), "\n")
print(table(guild_data$Habitat_Guild))

# Check which species have models
guild_data$has_model <- sapply(guild_data$species_folder, function(sp) {
  file.exists(file.path(hpc_paths$work_dir, sp, "proj_Current",
                        paste0("proj_Current_", sp, "_ensemble.tif")))
})
cat("With models:", sum(guild_data$has_model), "/", nrow(guild_data), "\n")
guild_data <- guild_data[guild_data$has_model, ]

# ---- Calculate per guild ----
for (guild in guilds) {
  sp_list <- guild_data$species_folder[guild_data$Habitat_Guild == guild]
  cat("\n=====", guild, "(", length(sp_list), "species) =====\n")
  if (length(sp_list) == 0) next

  for (period_label in names(periods)) {
    period_folder <- periods[period_label]
    cat("  Period:", period_label, "\n")
    prob_stack <- list()

    for (sp in sp_list) {
      prob_path <- file.path(hpc_paths$work_dir, sp,
                             paste0("proj_", period_folder),
                             paste0("proj_", period_folder, "_", sp, "_ensemble.tif"))
      if (file.exists(prob_path)) {
        tryCatch({
          r <- rast(prob_path)
          if (global(r, "max", na.rm = TRUE)[1, 1] > 10) r <- r / 1000
          prob_stack[[length(prob_stack) + 1]] <- r
        }, error = function(e) cat("    WARNING:", sp, "\n"))
      }
    }
    cat("    Loaded:", length(prob_stack), "species\n")

    if (length(prob_stack) > 0) {
      richness <- prob_stack[[1]]
      if (length(prob_stack) > 1)
        for (i in 2:length(prob_stack)) richness <- richness + prob_stack[[i]]

      guild_name <- gsub("-", "", guild)
      out_path <- file.path(output_dir,
                            paste0("richness_prob_", guild_name, "_", period_label, ".tif"))
      writeRaster(richness, out_path, overwrite = TRUE)
      cat("    Saved:", basename(out_path), "\n")
    }
    rm(prob_stack); gc()
  }

  # Change map
  guild_name <- gsub("-", "", guild)
  s3_f <- file.path(output_dir, paste0("richness_prob_", guild_name, "_S3.tif"))
  s0_f <- file.path(output_dir, paste0("richness_prob_", guild_name, "_S0.tif"))
  if (file.exists(s3_f) && file.exists(s0_f)) {
    change <- rast(s3_f) - rast(s0_f)
    writeRaster(change,
      file.path(output_dir, paste0("richness_prob_change_", guild_name, ".tif")),
      overwrite = TRUE)
    cat("  Change saved\n")
  }
}

cat("\n=== Guild richness complete for", COUNTRY, "===\n")
print(list.files(output_dir, pattern = "\\.tif$"))
