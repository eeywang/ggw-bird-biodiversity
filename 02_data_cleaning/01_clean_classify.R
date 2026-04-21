# =============================================================================
# 01_clean_classify.R — Coordinate cleaning + migration status classification
#                        for 3 focal countries (Nigeria, Senegal, Ethiopia)
#
# 7-step cleaning: missing coords → equal lat/lon → zero coords →
#   bounding box → centroid → spatio-temporal dedup → 5 km thinning
#
# Outputs per country in {step1_root}/{Country}/step1_review/:
#   - {country}_cleaned_step1.csv
#   - {country}_cleaning_log.csv
#   - {country}_migration_all_species.csv
#   - {country}_REVIEW_unknown_species.csv   ← inspect these!
#   - {country}_REVIEW_full_migrants.csv     ← inspect these!
#   - {country}_TO_REMOVE_non_african.csv
#   - {country}_TO_REMOVE_extinct.csv
# =============================================================================

# Source setup (works whether run from project root or this subfolder)
if (!file.exists("00_setup.R") && file.exists("../00_setup.R")) setwd("..")
source("00_setup.R")
source("functions.R")

cat("\n")
cat("============================================================\n")
cat("  Step 1: Coordinate Cleaning + Migration Classification\n")
cat("  Focal countries:", paste(focal_country_names, collapse = ", "), "\n")
cat("============================================================\n")

results_step1 <- list()

for (country_name in focal_country_names) {

  params <- focal_countries[[country_name]]
  output_dir <- file.path(paths$step1_root, country_name, "step1_review")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n", paste(rep("*", 70), collapse = ""), "\n")
  cat("  Processing:", country_name, "\n")
  cat(paste(rep("*", 70), collapse = ""), "\n")

  # ---- Read combined data ----
  data_file <- file.path(paths$combined,
                         paste0(tolower(country_name), "_birds_combined.rds"))
  if (!file.exists(data_file))
    data_file <- sub("\\.rds$", ".csv", data_file)
  if (!file.exists(data_file)) {
    cat("  Data not found:", data_file, "\n  Skipping.\n"); next
  }

  data_raw <- if (grepl("\\.rds$", data_file)) readRDS(data_file)
              else read_csv(data_file, show_col_types = FALSE)
  cat("  Loaded:", format(nrow(data_raw), big.mark = ","), "records,",
      n_distinct(data_raw$species), "species\n")

  # ---- Coordinate cleaning ----
  cleaning_result <- clean_coordinates_pipeline(data_raw, params)
  data_cleaned <- cleaning_result$data

  # ---- Migration classification ----
  migration_result <- classify_migration_status(
    data_cleaned, paths$avonet, country_name, shared_lists
  )

  # ---- Save outputs ----
  cat("\n  Saving step 1 outputs...\n")

  write_csv(data_cleaned,
    file.path(output_dir, paste0(tolower(country_name), "_cleaned_step1.csv")))
  write_csv(cleaning_result$cleaning_log,
    file.path(output_dir, paste0(tolower(country_name), "_cleaning_log.csv")))
  write_csv(migration_result$migration_table,
    file.path(output_dir, paste0(tolower(country_name), "_migration_all_species.csv")))

  write_csv(migration_result$unknown_species,
    file.path(output_dir, paste0(tolower(country_name), "_REVIEW_unknown_species.csv")))
  cat("    REVIEW — unknown:", nrow(migration_result$unknown_species), "\n")

  write_csv(migration_result$full_migrants,
    file.path(output_dir, paste0(tolower(country_name), "_REVIEW_full_migrants.csv")))
  cat("    REVIEW — full migrants:", nrow(migration_result$full_migrants), "\n")

  write_csv(migration_result$non_african,
    file.path(output_dir, paste0(tolower(country_name), "_TO_REMOVE_non_african.csv")))
  write_csv(migration_result$extinct,
    file.path(output_dir, paste0(tolower(country_name), "_TO_REMOVE_extinct.csv")))

  results_step1[[country_name]] <- list(
    data_cleaned = data_cleaned,
    migration_result = migration_result,
    cleaning_log = cleaning_result$cleaning_log
  )
  cat("  Done.\n")
}

cat("\n============================================================\n")
cat("  Step 1 complete. Next:\n")
cat("  1. Inspect *_REVIEW_*.csv files\n")
cat("  2. Record decisions in manual_review_decisions.R\n")
cat("  3. Run 02_data_cleaning/02_apply_review.R\n")
cat("============================================================\n")
