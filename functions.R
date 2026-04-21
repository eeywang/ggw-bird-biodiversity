# =============================================================================
# functions.R — Shared helper functions for GGW bird data pipeline
# =============================================================================

# ---- 1. eBird standardisation ----

#' Standardise raw eBird data to GBIF-compatible format
#'
#' Reads a raw eBird .txt file, filters for approved species-level records,
#' and returns a data.frame with the same columns as the GBIF download.
#'
#' @param country_code ISO2 country code (e.g. "NG")
#' @param country_name Country name (e.g. "Nigeria")
#' @param ebird_pattern Folder/file name pattern under `ebird_raw_path`
#' @param ebird_raw_path Root folder containing eBird download folders
#' @param ebird_output_path Folder to write standardised .rds / .csv
#' @return data.frame or NULL on error
standardize_ebird <- function(country_code, country_name, ebird_pattern,
                              ebird_raw_path, ebird_output_path) {

  cat("\n--- eBird standardisation:", country_name, "---\n")

  output_rds <- file.path(ebird_output_path,
                          paste0("ebird_", tolower(country_name), "_standardized.rds"))
  if (file.exists(output_rds)) {
    cat("  Already processed; loading cached file.\n")
    return(readRDS(output_rds))
  }

  ebird_file <- file.path(ebird_raw_path, ebird_pattern, paste0(ebird_pattern, ".txt"))
  if (!file.exists(ebird_file)) {
    cat("  eBird file not found:", ebird_file, "\n  Skipping.\n")
    return(NULL)
  }

  tryCatch({
    cat("  Reading...\n")
    ebird_raw <- data.table::fread(
      ebird_file,
      select = ebird_cols,
      quote  = "",
      na.strings = c("", "NA")
    )
    cat("  Raw records:", format(nrow(ebird_raw), big.mark = ","), "\n")

    ebird_final <- ebird_raw %>%
      filter(APPROVED == 1, CATEGORY == "species") %>%
      transmute(
        species          = `SCIENTIFIC NAME`,
        decimalLatitude  = LATITUDE,
        decimalLongitude = LONGITUDE,
        eventDate        = as.Date(`OBSERVATION DATE`),
        basisOfRecord    = "HUMAN_OBSERVATION",
        datasetName      = "eBird",
        occurrenceID     = `GLOBAL UNIQUE IDENTIFIER`,
        gbifID           = `GLOBAL UNIQUE IDENTIFIER`,
        country          = country_name,
        countryCode      = country_code,
        stateProvince    = STATE,
        taxonRank        = "SPECIES",
        order            = NA_character_,
        family           = NA_character_,
        year             = year(eventDate)
      ) %>%
      filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
             year >= 2007, year <= 2024)

    saveRDS(ebird_final, output_rds)
    data.table::fwrite(ebird_final, sub("\\.rds$", ".csv", output_rds))

    cat("  Standardised records:", format(nrow(ebird_final), big.mark = ","),
        "| Species:", n_distinct(ebird_final$species), "\n")
    return(ebird_final)

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}


# ---- 2. GBIF + eBird merge & deduplication ----

#' Merge GBIF and eBird data for one country, deduplicating by
#' (species, lat, lon, date). GBIF records are always kept; only unique eBird
#' records are appended.
#'
#' @param country_code ISO2 code
#' @param country_name Country name
#' @param gbif_path Folder containing `gbif_{country}_standardized.rds`
#' @param ebird_data data.frame from `standardize_ebird()`
#' @param combined_output_path Folder to write combined .rds / .csv
#' @return data.frame or NULL
merge_and_deduplicate <- function(country_code, country_name,
                                  gbif_path, ebird_data,
                                  combined_output_path) {

  cat("  Merging GBIF + eBird:", country_name, "\n")

  output_rds <- file.path(combined_output_path,
                          paste0(tolower(country_name), "_birds_combined.rds"))
  if (file.exists(output_rds)) {
    cat("  Already merged; loading cached file.\n")
    return(readRDS(output_rds))
  }

  # Load GBIF
  gbif_file <- file.path(gbif_path,
                         paste0("gbif_", tolower(country_name), "_standardized.rds"))
  gbif_data <- if (file.exists(gbif_file)) readRDS(gbif_file) else NULL

  # Ensure countryCode column
 if (!is.null(gbif_data) && !"countryCode" %in% names(gbif_data)) {
    gbif_data$countryCode <- country_code
  }

  # Handle missing data
  if (is.null(gbif_data) && (is.null(ebird_data) || nrow(ebird_data) == 0)) {
    cat("  No data for", country_name, "\n")
    return(NULL)
  }
  if (is.null(gbif_data)) {
    saveRDS(ebird_data, output_rds)
    data.table::fwrite(ebird_data, sub("\\.rds$", ".csv", output_rds))
    return(ebird_data)
  }
  if (is.null(ebird_data) || nrow(ebird_data) == 0) {
    saveRDS(gbif_data, output_rds)
    data.table::fwrite(gbif_data, sub("\\.rds$", ".csv", output_rds))
    return(gbif_data)
  }

  n_gbif  <- nrow(gbif_data)
  n_ebird <- nrow(ebird_data)

  # Deduplicate: keep all GBIF, add only unique eBird
  gbif_data  <- gbif_data  %>% mutate(match_key = paste(species, decimalLatitude, decimalLongitude, eventDate, sep = "_"))
  ebird_data <- ebird_data %>% mutate(match_key = paste(species, decimalLatitude, decimalLongitude, eventDate, sep = "_"))

  ebird_unique <- ebird_data %>%
    filter(!match_key %in% gbif_data$match_key) %>%
    select(-match_key)
  gbif_data <- gbif_data %>% select(-match_key)

  combined <- bind_rows(gbif_data, ebird_unique)

  cat("    GBIF:", format(n_gbif, big.mark = ","),
      "| eBird unique:", format(nrow(ebird_unique), big.mark = ","),
      "| Total:", format(nrow(combined), big.mark = ","), "\n")

  saveRDS(combined, output_rds)
  data.table::fwrite(combined, sub("\\.rds$", ".csv", output_rds))
  return(combined)
}


# ---- 3. Coordinate cleaning ----

#' Seven-step coordinate cleaning for one country.
#' Steps: (1) missing coords, (2) equal lat/lon, (3) zero coords,
#' (4) geographic bounding box, (5) country centroid, (6) spatio-temporal
#' deduplication, (7) spatial thinning (~5 km grid).
#'
#' @param data Occurrence data.frame
#' @param params Country parameter list (from `focal_countries`)
#' @return list(data, cleaning_log, out_of_range)
clean_coordinates_pipeline <- function(data, params) {

  cat("\n=== Coordinate cleaning:", params$name, "===\n")
  cat("  Initial records:", format(nrow(data), big.mark = ","), "\n")

  log <- tibble(Step = character(), Description = character(),
                Removed = integer(), Remaining = integer())
  initial <- nrow(data)

  record_step <- function(step, desc, before, after) {
    removed <- before - after
    cat("  Step", step, "-", desc, ": removed", removed, "\n")
    log <<- bind_rows(log, tibble(
      Step = step, Description = desc,
      Removed = removed, Remaining = after
    ))
  }

  # Step 1: missing coordinates
  before <- nrow(data)
  data <- data %>% filter(!is.na(decimalLatitude), !is.na(decimalLongitude))
  record_step("1", "Missing coordinates", before, nrow(data))

  # Step 2: equal lat/lon
  before <- nrow(data)
  data <- data %>% filter(decimalLongitude != decimalLatitude)
  record_step("2", "Equal lat/lon", before, nrow(data))

  # Step 3: zero coordinates
  before <- nrow(data)
  data <- data %>%
    filter(!(decimalLongitude == 0 & decimalLatitude == 0),
           decimalLongitude != 0, decimalLatitude != 0)
  record_step("3", "Zero coordinates", before, nrow(data))

  # Step 4: geographic bounding box
  before <- nrow(data)
  out_of_range <- data %>%
    filter(decimalLongitude < params$lon_min | decimalLongitude > params$lon_max |
           decimalLatitude  < params$lat_min | decimalLatitude  > params$lat_max)
  data <- data %>%
    filter(decimalLongitude >= params$lon_min, decimalLongitude <= params$lon_max,
           decimalLatitude  >= params$lat_min, decimalLatitude  <= params$lat_max)
  record_step("4", "Outside geographic range", before, nrow(data))

  # Step 5: country centroid
  before <- nrow(data)
  buf <- params$centroid_buffer
  data <- data %>%
    filter(!(abs(decimalLongitude - params$centroid_lon) < buf &
             abs(decimalLatitude  - params$centroid_lat) < buf))
  record_step("5", "Country centroid", before, nrow(data))

  # Step 6: spatio-temporal duplicates
  before <- nrow(data)
  data <- data %>%
    distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)
  record_step("6", "Spatio-temporal duplicates", before, nrow(data))

  # Step 7: spatial thinning (~5 km)
  before <- nrow(data)
  res <- params$thinning_resolution
  data <- data %>%
    mutate(grid_lon = round(decimalLongitude / res) * res,
           grid_lat = round(decimalLatitude  / res) * res,
           grid_id  = paste(grid_lon, grid_lat, sep = "_")) %>%
    group_by(species, grid_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(-grid_lon, -grid_lat, -grid_id)
  record_step("7", "Spatial thinning (5 km)", before, nrow(data))

  cat("  Final:", format(nrow(data), big.mark = ","),
      "| Retention:", round(nrow(data) / initial * 100, 1), "%\n")

  list(data = data, cleaning_log = log, out_of_range = out_of_range)
}


# ---- 4. Migration status classification ----

#' Classify migration status for all species in a cleaned dataset using:
#'   (1) AVONET eBird sheet, (2) AVONET BirdLife sheet, (3) synonym mapping,
#'   (4) manually curated lists from `shared_lists`.
#'
#' Returns a list containing the full migration table and subsets for review.
#'
#' @param data Cleaned occurrence data.frame
#' @param avonet_file Path to AVONET Excel file
#' @param country_name Country name (for logging)
#' @param lists `shared_lists` object from 00_setup.R
#' @return list with migration_table, status_summary, unknown_species, etc.
classify_migration_status <- function(data, avonet_file, country_name, lists) {

  cat("\n=== Migration classification:", country_name, "===\n")

  # Requires readxl
  if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required")

  avonet_ebird    <- readxl::read_excel(avonet_file, sheet = "AVONET2_eBird")
  avonet_birdlife <- readxl::read_excel(avonet_file, sheet = "AVONET1_BirdLife")

  unique_species <- unique(data$species)
  cat("  Species to classify:", length(unique_species), "\n")

  # Step 1: eBird sheet match
  sp_mig <- tibble(species = unique_species) %>%
    left_join(avonet_ebird %>% select(Species2, Migration) %>%
              rename(species = Species2), by = "species")
  cat("  Step 1 (eBird):", sum(!is.na(sp_mig$Migration)), "matched\n")

  # Step 2: BirdLife sheet supplement
  not_found <- sp_mig %>% filter(is.na(Migration)) %>% pull(species)
  if (length(not_found) > 0) {
    bl_match <- tibble(species = not_found) %>%
      left_join(avonet_birdlife %>% select(Species1, Migration) %>%
                rename(species = Species1), by = "species")
    sp_mig <- sp_mig %>%
      left_join(bl_match %>% rename(Migration_BL = Migration), by = "species") %>%
      mutate(Migration = coalesce(Migration, Migration_BL)) %>%
      select(species, Migration)
    cat("  Step 2 (BirdLife):", sum(!is.na(bl_match$Migration)), "new matches\n")
  }

  # Step 3: synonym matching
  not_found2 <- sp_mig %>% filter(is.na(Migration)) %>% pull(species)
  synonym_sp <- not_found2[not_found2 %in% names(lists$synonym_map)]
  matched_syn <- 0
  for (sp in synonym_sp) {
    avonet_name <- lists$synonym_map[sp]
    mig_val <- avonet_ebird$Migration[avonet_ebird$Species2 == avonet_name]
    if (length(mig_val) == 0)
      mig_val <- avonet_birdlife$Migration[avonet_birdlife$Species1 == avonet_name]
    if (length(mig_val) > 0 && !is.na(mig_val[1])) {
      sp_mig$Migration[sp_mig$species == sp] <- mig_val[1]
      matched_syn <- matched_syn + 1
    }
  }
  cat("  Step 3 (synonyms):", matched_syn, "matched\n")

  # Add taxonomy
  sp_tax <- data %>% distinct(species, order, family) %>% filter(!is.na(species))
  sp_mig <- sp_mig %>% left_join(sp_tax, by = "species")

  # Step 4: final classification
  sp_mig <- sp_mig %>%
    mutate(
      Migration = as.numeric(Migration),
      Migration_final = case_when(
        species %in% lists$non_african_species    ~ 99L,
        species %in% lists$extinct_species         ~ 98L,
        species %in% lists$palearctic_winterers    ~ 2L,
        species %in% lists$species_keep_regardless ~ 2L,
        species %in% lists$manual_african_residents ~ 1L,
        !is.na(Migration)                          ~ as.integer(Migration),
        family %in% lists$african_resident_families ~ 1L,
        TRUE                                        ~ 0L
      ),
      Migration_status = case_when(
        Migration_final == 99L ~ "NON-AFRICAN (remove)",
        Migration_final == 98L ~ "EXTINCT (remove)",
        Migration_final == 1L  ~ "Sedentary",
        Migration_final == 2L  ~ "Partial migrant",
        Migration_final == 3L  ~ "Full migrant",
        Migration_final == 0L  ~ "UNKNOWN (needs review)",
        TRUE ~ "ERROR"
      ),
      Migration_source = case_when(
        species %in% lists$non_african_species     ~ "Non-African species list",
        species %in% lists$extinct_species          ~ "Extinct species list",
        species %in% lists$palearctic_winterers     ~ "Palearctic winterer list",
        species %in% lists$species_keep_regardless  ~ "Keep regardless (complex taxonomy)",
        species %in% lists$manual_african_residents  ~ "Manual African resident list",
        !is.na(Migration) & Migration_final %in% 1:3 ~ "AVONET",
        family %in% lists$african_resident_families  ~ "Family-based assumption",
        TRUE ~ "UNKNOWN - needs manual review"
      )
    )

  # Add per-species record counts
  sp_counts <- data %>% count(species, name = "n_records")
  sp_mig <- sp_mig %>% left_join(sp_counts, by = "species")

  # Summary
  status_summary <- sp_mig %>%
    group_by(Migration_status) %>%
    summarise(n_species = n(), n_records = sum(n_records, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(pct = round(n_species / sum(n_species) * 100, 1)) %>%
    arrange(desc(n_species))

  cat("\n  Migration status summary:\n")
  print(status_summary)

  unknown <- sp_mig %>% filter(Migration_status == "UNKNOWN (needs review)") %>%
    arrange(desc(n_records))
  full_mig <- sp_mig %>% filter(Migration_status == "Full migrant") %>%
    arrange(desc(n_records))
  cat("  UNKNOWN:", nrow(unknown), "| Full migrant:", nrow(full_mig), "\n")

  list(
    migration_table = sp_mig,
    status_summary  = status_summary,
    unknown_species = unknown,
    full_migrants   = full_mig,
    non_african     = sp_mig %>% filter(Migration_status == "NON-AFRICAN (remove)"),
    extinct         = sp_mig %>% filter(Migration_status == "EXTINCT (remove)")
  )
}


# ---- 5. Apply manual review decisions ----

#' Apply manual review decisions and produce final resident-bird dataset.
#'
#' @param data_cleaned Cleaned occurrence data (from step 1)
#' @param migration_table Full migration table (from step 1)
#' @param unknown_decisions data.frame(species, decision) for UNKNOWN species
#' @param migrants_to_keep Character vector of full-migrant species to retain
#' @param error_records Character vector of erroneous non-African full migrants
#' @param country_name Country name for logging
#' @return list(data_final, species_list, migration_table_updated, ...)
apply_manual_review <- function(data_cleaned, migration_table,
                                unknown_decisions, migrants_to_keep,
                                error_records, country_name) {

  cat("\n=== Applying manual review:", country_name, "===\n")

  # Parse unknown decisions
  unknown_to_remove <- unknown_decisions %>%
    filter(decision == "remove") %>% pull(species)
  unknown_keep_resident <- unknown_decisions %>%
    filter(decision == "keep_resident") %>% pull(species)
  unknown_keep_partial <- unknown_decisions %>%
    filter(decision == "keep_partial") %>% pull(species)

  # Full migrants: keep listed, remove the rest
  fm_species <- migration_table %>%
    filter(Migration_status == "Full migrant") %>%
    pull(species) %>% unique()
  fm_keep   <- intersect(fm_species, migrants_to_keep)
  fm_remove <- setdiff(fm_species, migrants_to_keep)
  fm_remove <- unique(c(fm_remove, error_records))

  # Already-flagged removals from step 1
  non_african <- migration_table %>%
    filter(Migration_status == "NON-AFRICAN (remove)") %>% pull(species)
  extinct <- migration_table %>%
    filter(Migration_status == "EXTINCT (remove)") %>% pull(species)

  all_remove <- unique(c(unknown_to_remove, fm_remove, non_african, extinct))
  all_remove <- all_remove[!is.na(all_remove)]

  cat("  Species to remove:", length(all_remove), "\n")

  # Filter data
  data_final <- data_cleaned %>%
    filter(!is.na(species), !species %in% all_remove)

  cat("  Records:", format(nrow(data_cleaned), big.mark = ","),
      "->", format(nrow(data_final), big.mark = ","), "\n")
  cat("  Species:", n_distinct(data_cleaned$species),
      "->", n_distinct(data_final$species), "\n")

  # Update migration table
  mig_updated <- migration_table %>%
    mutate(
      Migration_final_updated = case_when(
        species %in% unknown_keep_resident ~ 1L,
        species %in% unknown_keep_partial  ~ 2L,
        species %in% fm_keep               ~ 2L,
        TRUE ~ as.integer(Migration_final)
      ),
      Migration_status_updated = case_when(
        species %in% unknown_keep_resident ~ "Sedentary",
        species %in% unknown_keep_partial  ~ "Partial migrant",
        species %in% fm_keep               ~ "Partial migrant (African population)",
        species %in% all_remove            ~ "REMOVED",
        TRUE ~ Migration_status
      ),
      keep = !species %in% all_remove & !is.na(species)
    )

  final_species <- mig_updated %>%
    filter(keep) %>%
    select(species, order, family,
           Migration_final = Migration_final_updated,
           Migration_status = Migration_status_updated,
           Migration_source, n_records) %>%
    distinct() %>%
    arrange(order, family, species)

  final_status <- final_species %>%
    group_by(Migration_status) %>%
    summarise(n_species = n(), n_records = sum(n_records, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(pct = round(n_species / sum(n_species) * 100, 1)) %>%
    arrange(desc(n_species))
  cat("\n  Final migration status:\n")
  print(final_status)

  list(
    data_final       = data_final,
    species_list     = final_species,
    migration_table  = mig_updated,
    removed_species  = mig_updated %>% filter(!keep | is.na(species)),
    summary          = final_status
  )
}

cat("=== functions.R loaded ===\n")
