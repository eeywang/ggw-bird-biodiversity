# =============================================================================
# main.R — Master orchestration for GGW bird data pipeline
#
# Run each step in order by uncommenting the relevant source() line.
# Or source individual scripts directly from their subfolders.
#
# Before running:
#   1. Edit paths in 00_setup.R
#   2. Set GBIF credentials in .Renviron (for step 01 only)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  GGW Bird Data Pipeline — Great Green Wall Avian Biodiversity\n")
cat("================================================================\n\n")

# Step 01: Download GBIF data for 11 GGW countries
# source("01_data_acquisition/01_download_gbif.R")

# Step 02: Standardise eBird + merge with GBIF (11 countries)
# source("01_data_acquisition/02_standardize_merge.R")

# Step 03: Coordinate cleaning + migration classification (3 focal countries)
# source("02_data_cleaning/01_clean_classify.R")

# Step 04: Apply manual review → final resident data (3 focal countries)
# source("02_data_cleaning/02_apply_review.R")

cat("  Steps are commented out by default.\n")
cat("  Uncomment the step(s) you want to run.\n")
cat("================================================================\n")
