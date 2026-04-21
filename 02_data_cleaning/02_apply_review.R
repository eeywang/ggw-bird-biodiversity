# =============================================================================
# 02_apply_review.R — Apply manual review decisions, produce final
#                      resident bird datasets for SDM modelling
#
# Outputs per country in {step2_root}/{Country}/:
#   - {country}_final_residents.csv      ← FINAL data for biomod2
#   - {country}_final_species_list.csv
#   - {country}_migration_final.csv      — full audit trail
#   - {country}_removed_species.csv
# =============================================================================

# Source setup (works whether run from project root or this subfolder)
if (!file.exists("00_setup.R") && file.exists("../00_setup.R")) setwd("..")
source("00_setup.R")
source("functions.R")
source("manual_review_decisions.R")

cat("\n")
cat("============================================================\n")
cat("  Step 2: Apply Manual Review -> Final Resident Data\n")
cat("  Focal countries:", paste(focal_country_names, collapse = ", "), "\n")
cat("============================================================\n")

results_step2 <- list()

for (country_name in focal_country_names) {

  output_dir <- file.path(paths$step2_root, country_name)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  step1_dir <- file.path(paths$step1_root, country_name, "step1_review")

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("  Processing:", country_name, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  # Read step 1 outputs
  cleaned_file   <- file.path(step1_dir, paste0(tolower(country_name), "_cleaned_step1.csv"))
  migration_file <- file.path(step1_dir, paste0(tolower(country_name), "_migration_all_species.csv"))

  if (!file.exists(cleaned_file) || !file.exists(migration_file)) {
    cat("  Step 1 files not found. Skipping.\n"); next
  }

  data_cleaned    <- read_csv(cleaned_file, show_col_types = FALSE)
  migration_table <- read_csv(migration_file, show_col_types = FALSE)
  cat("  Loaded:", format(nrow(data_cleaned), big.mark = ","), "records\n")

  # Retrieve decisions
  decisions <- country_decisions[[country_name]]

  # Apply review
  result <- apply_manual_review(
    data_cleaned     = data_cleaned,
    migration_table  = migration_table,
    unknown_decisions = decisions$unknown_decisions,
    migrants_to_keep  = decisions$migrants_to_keep,
    error_records     = decisions$error_records,
    country_name      = country_name
  )

  # Save
  cat("\n  Saving final outputs...\n")

  write_csv(result$data_final,
    file.path(output_dir, paste0(tolower(country_name), "_final_residents.csv")))
  cat("    Final:", format(nrow(result$data_final), big.mark = ","),
      "records,", n_distinct(result$data_final$species), "species\n")

  write_csv(result$species_list,
    file.path(output_dir, paste0(tolower(country_name), "_final_species_list.csv")))
  write_csv(result$migration_table,
    file.path(output_dir, paste0(tolower(country_name), "_migration_final.csv")))

  removed <- result$removed_species %>%
    select(any_of(c("species", "order", "family",
                    "Migration_status", "Migration_status_updated", "n_records"))) %>%
    distinct() %>% arrange(desc(n_records))
  write_csv(removed,
    file.path(output_dir, paste0(tolower(country_name), "_removed_species.csv")))

  results_step2[[country_name]] <- result
}

# Final summary
cat("\n============================================================\n")
cat("  Pipeline complete!\n")
cat("============================================================\n\n")

for (cn in focal_country_names) {
  if (!is.null(results_step2[[cn]])) {
    cat(sprintf("  %s: %s records, %d species\n", cn,
                format(nrow(results_step2[[cn]]$data_final), big.mark = ","),
                n_distinct(results_step2[[cn]]$data_final$species)))
  }
}
cat("\n  Final data:", paths$step2_root, "\n")
cat("  Use *_final_residents.csv for SDM modelling with biomod2.\n")
cat("============================================================\n")
