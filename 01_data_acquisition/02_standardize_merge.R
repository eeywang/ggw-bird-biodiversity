# =============================================================================
# 02_standardize_merge.R — Standardise eBird + merge with GBIF
#                           for all 11 GGW countries
#
# Outputs per country:
#   - {ebird_raw}/ebird_{country}_standardized.rds / .csv
#   - {combined}/{country}_birds_combined.rds / .csv
#   - {combined}/GGW_all_countries_birds_combined.rds / .csv
#   - {combined}/country_summary_statistics.csv
# =============================================================================

# Source setup (works whether run from project root or this subfolder)
if (!file.exists("00_setup.R") && file.exists("../00_setup.R")) setwd("..")
source("00_setup.R")
source("functions.R")

cat("\n")
cat("============================================================\n")
cat("  eBird Standardisation + GBIF Merge — 11 GGW Countries\n")
cat("============================================================\n\n")

all_results  <- list()
summary_list <- list()

for (i in seq_len(nrow(countries_all))) {

  cc   <- countries_all$code[i]
  cn   <- countries_all$name[i]
  epat <- countries_all$ebird_pattern[i]

  cat(sprintf("\n[%d/%d] %s (%s)\n", i, nrow(countries_all), cn, cc))
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Standardise eBird
  ebird_data <- standardize_ebird(
    country_code    = cc,
    country_name    = cn,
    ebird_pattern   = epat,
    ebird_raw_path  = paths$ebird_raw,
    ebird_output_path = paths$ebird_raw
  )

  # Merge with GBIF
  combined_data <- merge_and_deduplicate(
    country_code         = cc,
    country_name         = cn,
    gbif_path            = paths$gbif_raw,
    ebird_data           = ebird_data,
    combined_output_path = paths$combined
  )

  if (!is.null(combined_data)) {
    all_results[[cn]] <- combined_data

    gbif_file <- file.path(paths$gbif_raw,
                           paste0("gbif_", tolower(cn), "_standardized.rds"))
    n_gbif  <- if (file.exists(gbif_file)) nrow(readRDS(gbif_file)) else 0
    n_ebird <- if (!is.null(ebird_data)) nrow(ebird_data) else 0

    summary_list[[i]] <- tibble(
      country = cn, code = cc,
      n_gbif = n_gbif, n_ebird = n_ebird,
      n_combined = nrow(combined_data),
      n_species  = n_distinct(combined_data$species)
    )
  }
}

# Summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

summary_df <- bind_rows(summary_list)
print(summary_df)
fwrite(summary_df, file.path(paths$combined, "country_summary_statistics.csv"))

# Merge all countries
cat("\nMerging all countries...\n")
all_combined <- bind_rows(all_results)
cat("  Total records:", format(nrow(all_combined), big.mark = ","), "\n")
cat("  Total species:", n_distinct(all_combined$species), "\n")

saveRDS(all_combined, file.path(paths$combined, "GGW_all_countries_birds_combined.rds"))
fwrite(all_combined,  file.path(paths$combined, "GGW_all_countries_birds_combined.csv"))

# Source breakdown
cat("\nRecords by country and source:\n")
source_stats <- all_combined %>%
  group_by(country, countryCode, datasetName) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = datasetName, values_from = n, values_fill = 0) %>%
  mutate(total = GBIF + eBird) %>%
  arrange(desc(total))
print(source_stats)
fwrite(source_stats, file.path(paths$combined, "source_statistics_by_country.csv"))

# Focal countries
cat("\nFocal countries:\n")
focal_detail <- all_combined %>%
  filter(country %in% focal_country_names) %>%
  group_by(country) %>%
  summarise(n_records = n(), n_species = n_distinct(species),
            pct_of_total = round(n() / nrow(all_combined) * 100, 1),
            .groups = "drop") %>%
  arrange(desc(n_records))
print(focal_detail)
fwrite(focal_detail, file.path(paths$combined, "focal_countries_statistics.csv"))

cat("\n  Done. Output folder:", paths$combined, "\n")
