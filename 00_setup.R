# =============================================================================
# 00_setup.R — Central configuration for GGW bird data pipeline
# All paths, country parameters, and shared reference lists
#
# Usage: source("00_setup.R") at the top of every script
#        (always run from Github_coding/ as working directory)
# =============================================================================

# ---- Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Resolve common namespace conflicts
select  <- dplyr::select
filter  <- dplyr::filter

# ---- Root directory (auto-detect or set manually) ----
if (requireNamespace("here", quietly = TRUE)) {
  project_root <- here::here()
} else {
  # *** MODIFY THIS to match your local setup ***
  project_root <- "D:/OneDrive - University of Leeds/biodiversity_manscript/Github_coding"
}

# ---- File paths ----
# All downstream scripts reference `paths$...` — only edit paths here.
paths <- list(
  # Raw data inputs
  ebird_raw      = "E:/eBird_data",
  gbif_raw       = "E:/GBIF_data",
  avonet         = "E:/12.1progress_biomod2/AVONET_Supplementary_dataset_1.xlsx",

  # Intermediate outputs
  combined       = "E:/bird_dataset_combined",
  step1_root     = "E:/GGW_bird_analysis",
  step2_root     = "E:/GGW_bird_analysis/step2_final",

  # Within-repo processed data
  data_processed = file.path(project_root, "data_processed"),

  # Boundary shapefiles (for QC plots; optional)
  boundaries = list(
    Nigeria  = list(
      country = "E:/2025.11.2progress/studyarea/Nigeria_boundary.shp",
      ggw     = "E:/2025.11.2progress/studyarea/Nigeria_GGW_boundary.shp"
    ),
    Senegal  = list(
      country = "E:/GGW_boundaries/Senegal_boundary.shp",
      ggw     = "E:/GGW_boundaries/Senegal_GGW_boundary.shp"
    ),
    Ethiopia = list(
      country = "E:/GGW_boundaries/Ethiopia_boundary.shp",
      ggw     = "E:/GGW_boundaries/Ethiopia_GGW_boundary.shp"
    )
  )
)

# Create output directories if they don't exist
for (d in c(paths$combined, paths$step1_root, paths$step2_root, paths$data_processed)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- GBIF credentials ----
# Set your own credentials in .Renviron, NOT here (avoid committing to Git)
# GBIF_USER  = "your_username"
# GBIF_PWD   = "your_password"
# GBIF_EMAIL = "your_email"

# ---- 11 GGW countries ----
countries_all <- tibble::tibble(
  code = c("SN", "MR", "ML", "BF", "NE", "NG", "TD", "SD", "ER", "ET", "DJ"),
  name = c("Senegal", "Mauritania", "Mali", "Burkina Faso", "Niger",
           "Nigeria", "Chad", "Sudan", "Eritrea", "Ethiopia", "Djibouti"),
  ebird_pattern = c(
    "ebd_SN_200701_202412_smp_relOct-2025",
    "ebd_MR_200701_202412_smp_relOct-2025",
    "ebd_ML_200701_202412_smp_relOct-2025",
    "ebd_BF_200701_202412_smp_relOct-2025",
    "ebd_NE_200701_202412_smp_relOct-2025",
    "ebd_NG_200701_202412_smp_relOct-2025",
    "ebd_TD_200701_202412_smp_relOct-2025",
    "ebd_SD_200701_202412_smp_relOct-2025",
    "ebd_ER_200701_202412_smp_relOct-2025",
    "ebd_ET_200001_202412_smp_relOct-2025",
    "ebd_DJ_200701_202412_smp_relOct-2025"
  )
)

# ---- 3 focal countries — cleaning parameters ----
focal_country_names <- c("Nigeria", "Senegal", "Ethiopia")

focal_countries <- list(
  Nigeria = list(
    name = "Nigeria", code = "NG",
    lon_min =  2.5,  lon_max = 15.0,
    lat_min =  4.0,  lat_max = 14.5,
    centroid_lon = 8.68, centroid_lat = 9.08,
    centroid_buffer = 0.05,
    thinning_resolution = 0.05
  ),
  Senegal = list(
    name = "Senegal", code = "SN",
    lon_min = -17.6, lon_max = -11.3,
    lat_min =  12.3, lat_max =  16.7,
    centroid_lon = -14.45, centroid_lat = 14.50,
    centroid_buffer = 0.05,
    thinning_resolution = 0.05
  ),
  Ethiopia = list(
    name = "Ethiopia", code = "ET",
    lon_min = 33.0,  lon_max = 48.0,
    lat_min =  3.4,  lat_max = 15.0,
    centroid_lon = 38.75, centroid_lat = 9.02,
    centroid_buffer = 0.05,
    thinning_resolution = 0.05
  )
)

# ---- Shared reference lists ----
shared_lists <- list(

  # -- Taxonomic synonym map (old name -> AVONET name) --
  synonym_map = c(
    # Ardeidae
    "Ardea ibis"               = "Bubulcus ibis",
    "Egretta intermedia"       = "Ardea intermedia",
    "Dupetor flavicollis"      = "Ixobrychus flavicollis",
    "Butorides atricapilla"    = "Butorides striata",
    "Botaurus sturmii"         = "Ixobrychus sturmii",
    "Botaurus minutus"         = "Ixobrychus minutus",
    # Raptors
    "Aquila clanga"            = "Clanga clanga",
    "Aquila pomarina"          = "Clanga pomarina",
    "Milvus aegyptius"         = "Milvus migrans",
    "Gyps rueppellii"          = "Gyps rueppelli",
    "Astur melanoleucus"       = "Accipiter melanoleucus",
    # Rallidae
    "Amaurornis flavirostra"   = "Zapornia flavirostra",
    "Porzana pusilla"          = "Zapornia pusilla",
    "Porzana parva"            = "Zapornia parva",
    "Gallinula angulata"       = "Paragallinula angulata",
    # Paridae
    "Parus albiventris"        = "Melaniparus albiventris",
    "Parus funereus"           = "Melaniparus funereus",
    "Parus guineensis"         = "Melaniparus guineensis",
    "Parus leucomelas"         = "Melaniparus leucomelas",
    # Sylviidae
    "Sylvia communis"          = "Curruca communis",
    "Sylvia curruca"           = "Curruca curruca",
    # Others
    "Upupa africana"           = "Upupa epops",
    "Spilopelia senegalensis"  = "Streptopelia senegalensis",
    "Bradornis pallidus"       = "Melaenornis pallidus",
    "Icthyophaga vocifer"      = "Haliaeetus vocifer",
    "Ketupa lactea"            = "Bubo lacteus",
    "Tachyspiza badia"         = "Accipiter badius",
    "Tachyspiza erythropus"    = "Accipiter erythropus",
    "Anarhynchus marginatus"   = "Charadrius marginatus",
    "Thinornis tricollaris"    = "Charadrius tricollaris",
    "Thinornis forbesi"        = "Charadrius forbesi",
    "Cecropis melanocrissus"   = "Cecropis senegalensis",
    # Additional for Senegal/Ethiopia
    "Cercotrichas galactotes"  = "Erythropygia galactotes",
    "Iduna pallida"            = "Hippolais pallida",
    "Curruca nana"             = "Sylvia nana",
    "Curruca cantillans"       = "Sylvia cantillans"
  ),

  # -- Non-African species (erroneous records, always remove) --
  non_african_species = c(
    "Amaurornis akool", "Amaurornis olivieri", "Canirallus kioloides",
    "Fulica caribaea", "Gallinula melanops", "Gallinula mortierii",
    "Gallinula pacifica", "Gallinula ventralis", "Gallirallus striatus",
    "Laterallus spilonotus", "Nesoclopeus poecilopterus", "Porphyrio mantelli",
    "Porzana atra", "Porzana bicolor", "Porzana cinerea", "Porzana flaviventer",
    "Porzana fusca", "Porzana paykullii", "Porzana tabuensis",
    "Grus antigone", "Grus canadensis", "Grus leucogeranus",
    "Grus rubicunda", "Grus vipio",
    "Phalacrocorax aristotelis", "Phalacrocorax auritus",
    "Phalacrocorax brasilianus", "Phalacrocorax gaimardi",
    "Phalacrocorax magellanicus", "Phalacrocorax pelagicus",
    "Phalacrocorax penicillatus", "Phalacrocorax urile",
    "Rallus limicola", "Botaurus lentiginosus", "Ixobrychus sinensis",
    "Chroicocephalus saundersi", "Vanellus macropterus",
    "Phylloscopus sibillatrix"
  ),

  # -- Extinct species --
  extinct_species = c(
    "Rhodonessa caryophyllacea", "Haematopus meadewaldoi",
    "Tadorna cristata", "Prosobonia cancellata", "Prosobonia leucoptera"
  ),

  # -- Palearctic winterers (keep as Partial migrant) --
  palearctic_winterers = c(
    "Hirundo rustica", "Delichon urbicum", "Riparia riparia",
    "Ficedula hypoleuca", "Muscicapa striata", "Phoenicurus phoenicurus",
    "Oenanthe oenanthe", "Saxicola rubetra",
    "Phylloscopus trochilus", "Phylloscopus sibilatrix", "Phylloscopus collybita",
    "Sylvia borin", "Sylvia atricapilla", "Sylvia communis",
    "Curruca communis", "Curruca curruca", "Curruca cantillans",
    "Acrocephalus schoenobaenus", "Acrocephalus scirpaceus",
    "Hippolais icterina", "Hippolais polyglotta",
    "Motacilla flava", "Motacilla cinerea", "Motacilla alba",
    "Anthus trivialis", "Anthus pratensis", "Anthus cervinus",
    "Lanius collurio", "Lanius senator", "Lanius minor",
    "Crex crex", "Ciconia ciconia",
    "Circus aeruginosus", "Circus pygargus", "Circus macrourus",
    "Pernis apivorus", "Falco subbuteo", "Falco vespertinus",
    "Calidris minuta", "Calidris ferruginea", "Calidris temminckii",
    "Philomachus pugnax", "Tringa glareola", "Tringa ochropus",
    "Actitis hypoleucos", "Numenius phaeopus",
    "Merops apiaster", "Upupa epops", "Cuculus canorus",
    "Oriolus oriolus", "Cercotrichas galactotes"
  ),

  # -- African resident families (infer unknown species) --
  african_resident_families = c(
    "Estrildidae", "Cisticolidae", "Pycnonotidae", "Picidae",
    "Buphagidae", "Campephagidae", "Platysteiridae", "Malaconotidae",
    "Laniariidae", "Prionopidae", "Oriolidae", "Dicruridae",
    "Monarchidae", "Nectariniidae", "Ploceidae", "Viduidae",
    "Bucerotidae", "Coliidae", "Musophagidae", "Indicatoridae",
    "Lybiidae", "Sturnidae", "Ardeidae", "Podicipedidae"
  ),

  # -- Manually confirmed African residents --
  manual_african_residents = c(
    "Crex egregia", "Crecopsis egregia", "Porphyrio alleni",
    "Aenigmatolimnas marginalis", "Limnocorax flavirostris",
    "Sarothrura pulchra", "Sarothrura elegans",
    "Ardea purpurea", "Ardeola ralloides", "Nycticorax nycticorax",
    "Ixobrychus minutus", "Ixobrychus sturmii",
    "Milvus aegyptius", "Aquila rapax", "Circaetus spectabilis",
    "Phalacrocorax lucidus", "Phalacrocorax africanus",
    "Leptoptilos crumeniferus", "Mycteria ibis",
    "Rhinoptilus africanus", "Rhinoptilus chalcopterus",
    "Charadrius marginatus", "Charadrius tricollaris", "Charadrius forbesi",
    "Riparia cincta", "Petrochelidon fuliginosa",
    "Cecropis senegalensis", "Cecropis abyssinica",
    "Muscicapa caerulescens", "Muscicapa comitata",
    "Muscicapa infuscata", "Muscicapa olivascens", "Muscicapa ussheri",
    "Myioparus griseigularis", "Myioparus plumbeus",
    "Bradornis pallidus", "Melaenornis pallidus",
    "Pseudalethe poliocephala", "Pseudoalcippe abyssinica",
    "Cossyphicula isabellae", "Dioptrornis chocolatinus",
    "Erythropygia hartlaubi", "Pentholaea albifrons", "Thamnolaea coronata",
    "Salpornis spilonotus", "Anthus longicaudatus",
    "Icthyophaga vocifer", "Ketupa lactea", "Tachyspiza badia",
    "Tachyspiza erythropus", "Tauraco violaceus",
    "Anarhynchus marginatus", "Thinornis tricollaris", "Thinornis forbesi",
    "Crecopsis egregia", "Artomyias fuliginosa", "Aerospiza castanilius",
    "Amirafra rufocinnamomea", "Corypha kurrae", "Thinornis dubius",
    "Cecropis melanocrissus",
    "Poicephalus flavifrons", "Tauraco ruspolii", "Zavattariornis stresemanni"
  ),

  # -- Species to keep regardless --
  species_keep_regardless = c(
    "Erythropygia galactotes"
  )
)

# ---- Standardised column order ----
fields_needed <- c(
  "species", "decimalLatitude", "decimalLongitude", "eventDate",
  "basisOfRecord", "datasetName", "occurrenceID", "gbifID",
  "country", "countryCode", "stateProvince", "taxonRank",
  "order", "family", "year"
)

# ---- eBird columns to read ----
ebird_cols <- c(
  "GLOBAL UNIQUE IDENTIFIER", "SCIENTIFIC NAME", "CATEGORY",
  "LATITUDE", "LONGITUDE", "OBSERVATION DATE",
  "COUNTRY", "STATE", "APPROVED"
)

cat("=== 00_setup.R loaded ===\n")
cat("  Project root:", project_root, "\n")
cat("  Focal countries:", paste(focal_country_names, collapse = ", "), "\n")
