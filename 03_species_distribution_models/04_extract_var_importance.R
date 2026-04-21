#!/usr/bin/env Rscript
# =============================================================================
# 05_extract_var_importance.R — Extract and summarise variable importance
#                                across all modelled species
#
# Usage: Rscript 05_extract_var_importance.R Nigeria
# =============================================================================

source("00_hpc_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(dplyr)
})

species_dirs <- list.dirs(hpc_paths$work_dir, recursive = FALSE)
species_dirs <- species_dirs[grepl("/[A-Z]", species_dirs)]

all_importance <- list()

for (sp_dir in species_dirs) {
  model_file <- list.files(sp_dir, pattern = "\\.models\\.out$", full.names = TRUE)
  model_file <- model_file[!grepl("ensemble", model_file)]

  if (length(model_file) > 0) {
    tryCatch({
      myModel <- get(load(model_file[1]))
      var_imp <- get_variables_importance(myModel)
      var_imp$species <- basename(sp_dir)
      var_imp$country <- COUNTRY
      all_importance[[basename(sp_dir)]] <- var_imp
    }, error = function(e) NULL)
  }
}

result <- bind_rows(all_importance)
write.csv(result, file.path(hpc_paths$results, "all_species_variable_importance.csv"),
          row.names = FALSE)

summary_imp <- result %>%
  group_by(expl.var) %>%
  summarise(mean_importance = mean(var.imp, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_importance))

cat("\nVariable importance summary for", COUNTRY, ":\n")
print(summary_imp)
write.csv(summary_imp, file.path(hpc_paths$results, "variable_importance_summary.csv"),
          row.names = FALSE)
