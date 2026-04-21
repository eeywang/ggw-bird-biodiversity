# =============================================================================
# 01_download_gbif.R — Download bird occurrence data from GBIF
#                       for 11 GGW countries (Aves, 2007-2024)
#
# Outputs: {gbif_raw}/gbif_{country}_standardized.rds / .csv
#
# Note: Each GBIF download request takes ~5-15 min to prepare.
#       The script skips already-downloaded countries.
# =============================================================================

# Source setup (works whether run from project root or this subfolder)
if (!file.exists("00_setup.R") && file.exists("../00_setup.R")) setwd("..")
source("00_setup.R")

library(rgbif)

cat("\n")
cat("============================================================\n")
cat("  GBIF Download — 11 GGW Countries (Aves, 2007-2024)\n")
cat("============================================================\n\n")

gbif_country_codes <- c(
  Senegal = "SN", Mauritania = "MR", Mali = "ML",
  `Burkina Faso` = "BF", Niger = "NE", Nigeria = "NG",
  Chad = "TD", Sudan = "SD", Eritrea = "ER",
  Ethiopia = "ET", Djibouti = "DJ"
)

dir.create(paths$gbif_raw, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(gbif_country_codes)) {

  country_name <- names(gbif_country_codes)[i]
  country_code <- gbif_country_codes[i]

  cat(sprintf("\n[%d/%d] %s (%s)\n", i, length(gbif_country_codes),
              country_name, country_code))
  cat(paste(rep("-", 50), collapse = ""), "\n")

  output_rds <- file.path(paths$gbif_raw,
                          paste0("gbif_", tolower(country_name), "_standardized.rds"))
  if (file.exists(output_rds)) {
    cat("  Already exists; skipping.\n")
    next
  }

  # Request download
  cat("  Requesting GBIF download...\n")
  download_key <- occ_download(
    pred("country", country_code),
    pred("classKey", 212),        # Aves
    pred("hasCoordinate", TRUE),
    pred_gte("year", 2007),
    pred_lte("year", 2024),
    format = "SIMPLE_CSV"
  )

  cat("  Download key:", download_key, "\n")
  cat("  Waiting for GBIF to prepare data...\n")
  occ_download_wait(download_key, status_ping = 30)

  # Download & import
  cat("  Downloading...\n")
  data_file <- occ_download_get(download_key, path = paths$gbif_raw)
  gbif_raw  <- occ_download_import(data_file)

  cat("  Raw records:", format(nrow(gbif_raw), big.mark = ","), "\n")

  # Standardise columns
  gbif_final <- gbif_raw %>%
    transmute(
      species, decimalLatitude, decimalLongitude,
      eventDate = as.Date(eventDate),
      basisOfRecord, datasetName = "GBIF",
      occurrenceID, gbifID = as.character(gbifID),
      country = country_name, countryCode = country_code,
      stateProvince, taxonRank, order, family, year
    ) %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude))

  cat("  Standardised:", format(nrow(gbif_final), big.mark = ","),
      "records |", n_distinct(gbif_final$species), "species\n")

  saveRDS(gbif_final, output_rds)
  fwrite(gbif_final, sub("\\.rds$", ".csv", output_rds))
  cat("  Saved.\n")
}

cat("\n============================================================\n")
cat("  GBIF download complete.\n")
cat("  Output folder:", paths$gbif_raw, "\n")
cat("============================================================\n")
