# =============================================================================
# + PSM+ 

# 1: 2000-2004bio1bio12
# 2: 2000-2024(Δbio1, Δbio12)
# 3: Full Matching1:1 Matching
# 4: 
# =============================================================================

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(exactextractr)
library(car)
library(MatchIt)
library(cobalt)
library(weights)
library(survey)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ncdf4)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  :  + PSM +                          ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# =============================================================================

# =============================================================================

cat("【】\n")
cat(rep("=", 60), "\n", sep = "")


nigeria_shp <- "E:/2025.11.2progress/studyarea/Nigeria_boundary.shp"
adm2_file <- "E:/11.17progress/study_area/Nigeria_adm2.shp"


bio_dir <- "E:/12.1progress_biomod2/ENV"
temp_dir <- file.path(bio_dir, "TerraClimate_raw")
trend_dir <- file.path(bio_dir, "climate_trend")
cov_dir <- "E:/2026.1.8_biomod2/NGA_covariates"
richness_dir <- "E:/12.1progress_biomod2/result/richness"


output_base <- "E:/2026.1.8_biomod2/nigeria_results"
output_full <- file.path(output_base, "PSM_full_with_trend")
output_1to1 <- file.path(output_base, "PSM_1to1_with_trend")
output_scenario <- file.path(output_base, "scenario_decomposition_with_trend")

dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(trend_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_full, showWarnings = FALSE, recursive = TRUE)
dir.create(output_1to1, showWarnings = FALSE, recursive = TRUE)
dir.create(output_scenario, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2000-2004TerraClimatebio1bio12
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  1: 2000-2004bio1bio12                     ║\n")
cat("╚", rep("═", 70), "╝\n\n")


nigeria_sf <- st_read(nigeria_shp, quiet = TRUE)
nigeria_vect <- vect(nigeria_sf)

bbox <- st_bbox(nigeria_sf)
nigeria_extent <- c(
  floor(bbox["xmin"]) - 1,
  ceiling(bbox["xmax"]) + 1,
  floor(bbox["ymin"]) - 1,
  ceiling(bbox["ymax"]) + 1
)

cat("  : ", round(nigeria_extent[1],1), "-", round(nigeria_extent[2],1),
    ", ", round(nigeria_extent[3],1), "-", round(nigeria_extent[4],1), "\n")

# TerraClimate
download_terraclimate <- function(variable, year, extent, output_folder) {
  base_url <- "https://climate.northwestknowledge.net/TERRACLIMATE-DATA"
  filename <- sprintf("TerraClimate_%s_%d.nc", variable, year)
  url <- paste0(base_url, "/", filename)
  
  local_file <- file.path(output_folder, filename)
  
  if (file.exists(local_file) && file.size(local_file) > 1000000) {
    r <- tryCatch({
      rast(local_file)
    }, error = function(e) {
      file.remove(local_file)
      return(NULL)
    })
    
    if (!is.null(r)) {
      ext_obj <- ext(extent[1], extent[2], extent[3], extent[4])
      return(crop(r, ext_obj))
    }
  }
  
  cat("    :", filename, "...")
  tryCatch({
    download.file(url, local_file, mode = "wb", quiet = TRUE)
    cat(" ✓\n")
  }, error = function(e) {
    cat(" ✗\n")
    stop(sprintf(": %s", e$message))
  })
  
  r <- rast(local_file)
  ext_obj <- ext(extent[1], extent[2], extent[3], extent[4])
  return(crop(r, ext_obj))
}

# bio1 () bio12 () 
calculate_bio1_bio12 <- function(tmin, tmax, prec) {
  # Bio1: = mean((tmax + tmin) / 2)
  tmean <- (tmin + tmax) / 2
  bio1 <- mean(tmean)
  
  # Bio12: = sum(prec)
  bio12 <- sum(prec)
  
  return(list(bio1 = bio1, bio12 = bio12))
}

# 2000-2004
years_to_download <- 2000:2004

cat("\n  2000-2004bio1bio12...\n")

for (year in years_to_download) {
  # BIO
  bio_file <- file.path(bio_dir, sprintf("BIO_%d.tif", year))
  
  if (file.exists(bio_file)) {
    cat("  ", year, ": BIO\n")
    next
  }
  
  cat("  ", year, "...\n")
  
  tryCatch({

    tmin <- download_terraclimate("tmin", year, nigeria_extent, temp_dir)
    tmax <- download_terraclimate("tmax", year, nigeria_extent, temp_dir)
    prec <- download_terraclimate("ppt", year, nigeria_extent, temp_dir)
    
    # TerraClimate°C*10
    tmin <- tmin / 10
    tmax <- tmax / 10
    
    # bio1bio12
    bio <- calculate_bio1_bio12(tmin, tmax, prec)
    

    bio1 <- crop(bio$bio1, nigeria_vect)
    bio1 <- mask(bio1, nigeria_vect)
    bio12 <- crop(bio$bio12, nigeria_vect)
    bio12 <- mask(bio12, nigeria_vect)
    
    # bio1bio12
    writeRaster(bio1, file.path(trend_dir, sprintf("bio1_%d.tif", year)), overwrite = TRUE)
    writeRaster(bio12, file.path(trend_dir, sprintf("bio12_%d.tif", year)), overwrite = TRUE)
    
    cat("    ✓ : bio1_", year, ".tif, bio12_", year, ".tif\n", sep = "")
    
    gc()
    
  }, error = function(e) {
    cat("    ✗ :", e$message, "\n")
  })
}

# =============================================================================
# 2000-2024bio1bio12
# =============================================================================

cat("\n【2000-2024bio1bio12】\n")
cat(rep("=", 60), "\n", sep = "")

all_years <- 2000:2024
bio1_list <- list()
bio12_list <- list()

for (year in all_years) {
  
  # 1: trend_dir
  bio1_file <- file.path(trend_dir, sprintf("bio1_%d.tif", year))
  bio12_file <- file.path(trend_dir, sprintf("bio12_%d.tif", year))
  
  if (file.exists(bio1_file) && file.exists(bio12_file)) {
    bio1_list[[as.character(year)]] <- rast(bio1_file)
    bio12_list[[as.character(year)]] <- rast(bio12_file)
    next
  }
  
  # 2: BIO
  bio_file <- file.path(bio_dir, sprintf("BIO_%d.tif", year))
  
  if (file.exists(bio_file)) {
    tryCatch({
      bio <- rast(bio_file)
      bio1_list[[as.character(year)]] <- bio[[1]]   # BIO11
      bio12_list[[as.character(year)]] <- bio[[12]] # BIO1212
    }, error = function(e) {
      cat("  ⚠️", year, "\n")
    })
  }
}

cat("  Bio1:", length(bio1_list), "\n")
cat("  Bio12:", length(bio12_list), "\n")

if (length(bio1_list) > 0) {
  cat("  :", paste(names(bio1_list), collapse = ", "), "\n")
}

# =============================================================================
# (2000-2024)
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  2: 2000-2024 (Δbio1, Δbio12)              ║\n")
cat("╚", rep("═", 70), "╝\n\n")


calc_trend_raster <- function(rast_list, var_name) {
  
  years <- as.numeric(names(rast_list))
  n <- length(years)
  
  if (n < 10) {
    cat("  ⚠️", var_name, ":  (", n, ")10\n")
    return(NULL)
  }
  
  cat("  ", var_name, " (", min(years), "-", max(years), ", n=", n, ")...\n")
  

  ref <- rast_list[[1]]
  aligned <- lapply(rast_list, function(r) {
    if (!compareGeom(r, ref, stopOnError = FALSE)) {
      r <- resample(r, ref)
    }
    return(r)
  })
  

  stack <- rast(aligned)
  
  # ()
  slope_func <- function(vals) {
    if (sum(!is.na(vals)) < 10) return(NA)
    valid <- !is.na(vals)
    x <- years[valid]
    y <- vals[valid]
    n <- length(x)
    (n * sum(x * y) - sum(x) * sum(y)) / (n * sum(x^2) - sum(x)^2)
  }
  
  trend <- app(stack, slope_func)
  names(trend) <- paste0(var_name, "_trend")
  

  vals <- values(trend, na.rm = TRUE)
  cat("    Mean:", round(mean(vals), 5), 
      ", Range:", round(min(vals), 5), "to", round(max(vals), 5), "\n")
  
  return(trend)
}


bio1_trend <- calc_trend_raster(bio1_list, "bio1")
bio12_trend <- calc_trend_raster(bio12_list, "bio12")


if (!is.null(bio1_trend)) {
  writeRaster(bio1_trend, file.path(trend_dir, "bio1_trend_2000_2024.tif"), overwrite = TRUE)
  cat("  ✓ : bio1_trend_2000_2024.tif (, °C/year)\n")
}

if (!is.null(bio12_trend)) {
  writeRaster(bio12_trend, file.path(trend_dir, "bio12_trend_2000_2024.tif"), overwrite = TRUE)
  cat("  ✓ : bio12_trend_2000_2024.tif (, mm/year)\n")
}

# =============================================================================
# LGA
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  3: PSM                                  ║\n")
cat("╚", rep("═", 70), "╝\n\n")

cat("【LGA】\n")

nga_adm2 <- st_read(adm2_file, quiet = TRUE)
cat("  LGA:", nrow(nga_adm2), "\n")

# GGW262LGA
ggw_gid2 <- c(
  "NGA.2.1_1", "NGA.2.2_1", "NGA.2.3_1", "NGA.2.4_1", "NGA.2.5_1", "NGA.2.6_1",
  "NGA.2.7_1", "NGA.2.8_1", "NGA.2.9_1", "NGA.2.10_1", "NGA.2.11_1", "NGA.2.12_1",
  "NGA.2.13_1", "NGA.2.14_1", "NGA.2.15_1", "NGA.2.16_1", "NGA.2.17_1", "NGA.2.18_1",
  "NGA.2.19_1", "NGA.2.20_1", "NGA.2.21_1", "NGA.5.1_1", "NGA.5.2_1", "NGA.5.3_1",
  "NGA.5.4_1", "NGA.5.5_1", "NGA.5.6_1", "NGA.5.7_1", "NGA.5.8_1", "NGA.5.9_1",
  "NGA.5.10_1", "NGA.5.11_1", "NGA.5.12_1", "NGA.5.13_1", "NGA.5.14_1", "NGA.5.15_1",
  "NGA.5.16_1", "NGA.5.17_1", "NGA.5.18_1", "NGA.5.19_1", "NGA.5.20_1", "NGA.8.1_1",
  "NGA.8.2_1", "NGA.8.3_1", "NGA.8.4_1", "NGA.8.5_1", "NGA.8.6_1", "NGA.8.7_1",
  "NGA.8.8_1", "NGA.8.9_1", "NGA.8.10_1", "NGA.8.11_1", "NGA.8.12_1", "NGA.8.13_1",
  "NGA.8.14_1", "NGA.8.15_1", "NGA.8.16_1", "NGA.8.17_1", "NGA.8.18_1", "NGA.8.19_1",
  "NGA.8.20_1", "NGA.8.21_1", "NGA.8.22_1", "NGA.8.23_1", "NGA.8.24_1", "NGA.8.25_1",
  "NGA.8.26_1", "NGA.8.27_1", "NGA.8.28_1", "NGA.16.1_1", "NGA.16.2_1", "NGA.16.3_1",
  "NGA.16.4_1", "NGA.16.5_1", "NGA.16.6_1", "NGA.16.7_1", "NGA.16.8_1", "NGA.16.9_1",
  "NGA.16.10_1", "NGA.16.11_1", "NGA.18.1_1", "NGA.18.2_1", "NGA.18.3_1", "NGA.18.4_1",
  "NGA.18.5_1", "NGA.18.6_1", "NGA.18.7_1", "NGA.18.8_1", "NGA.18.9_1", "NGA.18.10_1",
  "NGA.18.11_1", "NGA.18.12_1", "NGA.18.13_1", "NGA.18.14_1", "NGA.18.15_1", "NGA.18.16_1",
  "NGA.18.17_1", "NGA.18.18_1", "NGA.18.19_1", "NGA.18.20_1", "NGA.18.21_1", "NGA.18.22_1",
  "NGA.18.23_1", "NGA.18.24_1", "NGA.18.25_1", "NGA.18.26_1", "NGA.18.27_1", "NGA.19.5_1",
  "NGA.19.17_1", "NGA.20.1_1", "NGA.20.2_1", "NGA.20.3_1", "NGA.20.4_1", "NGA.20.5_1",
  "NGA.20.6_1", "NGA.20.7_1", "NGA.20.8_1", "NGA.20.9_1", "NGA.20.10_1", "NGA.20.11_1",
  "NGA.20.12_1", "NGA.20.13_1", "NGA.20.14_1", "NGA.20.15_1", "NGA.20.16_1", "NGA.20.17_1",
  "NGA.20.18_1", "NGA.20.19_1", "NGA.20.20_1", "NGA.20.21_1", "NGA.20.22_1", "NGA.20.23_1",
  "NGA.20.24_1", "NGA.20.25_1", "NGA.20.26_1", "NGA.20.27_1", "NGA.20.28_1", "NGA.20.29_1",
  "NGA.20.30_1", "NGA.20.31_1", "NGA.20.32_1", "NGA.20.33_1", "NGA.20.34_1", "NGA.20.35_1",
  "NGA.20.36_1", "NGA.20.37_1", "NGA.20.38_1", "NGA.20.39_1", "NGA.20.40_1", "NGA.20.41_1",
  "NGA.20.42_1", "NGA.20.43_1", "NGA.20.44_1", "NGA.21.1_1", "NGA.21.2_1", "NGA.21.3_1",
  "NGA.21.4_1", "NGA.21.5_1", "NGA.21.6_1", "NGA.21.7_1", "NGA.21.8_1", "NGA.21.9_1",
  "NGA.21.10_1", "NGA.21.11_1", "NGA.21.12_1", "NGA.21.13_1", "NGA.21.14_1", "NGA.21.15_1",
  "NGA.21.16_1", "NGA.21.17_1", "NGA.21.18_1", "NGA.21.19_1", "NGA.21.20_1", "NGA.21.21_1",
  "NGA.21.22_1", "NGA.21.23_1", "NGA.21.24_1", "NGA.21.25_1", "NGA.21.26_1", "NGA.21.27_1",
  "NGA.21.28_1", "NGA.21.29_1", "NGA.21.30_1", "NGA.21.31_1", "NGA.21.32_1", "NGA.21.33_1",
  "NGA.21.34_1", "NGA.22.1_1", "NGA.22.2_1", "NGA.22.3_1", "NGA.22.4_1", "NGA.22.5_1",
  "NGA.22.6_1", "NGA.22.7_1", "NGA.22.8_1", "NGA.22.9_1", "NGA.22.10_1", "NGA.22.11_1",
  "NGA.22.12_1", "NGA.22.13_1", "NGA.22.14_1", "NGA.22.15_1", "NGA.22.16_1", "NGA.22.17_1",
  "NGA.22.18_1", "NGA.22.19_1", "NGA.22.20_1", "NGA.22.21_1", "NGA.27.2_1", "NGA.27.15_1",
  "NGA.27.21_1", "NGA.32.4_1", "NGA.32.7_1", "NGA.34.1_1", "NGA.34.2_1", "NGA.34.3_1",
  "NGA.34.4_1", "NGA.34.5_1", "NGA.34.6_1", "NGA.34.7_1", "NGA.34.8_1", "NGA.34.9_1",
  "NGA.34.10_1", "NGA.34.11_1", "NGA.34.12_1", "NGA.34.13_1", "NGA.34.14_1", "NGA.34.15_1",
  "NGA.34.16_1", "NGA.34.17_1", "NGA.34.18_1", "NGA.34.19_1", "NGA.34.20_1", "NGA.34.21_1",
  "NGA.34.22_1", "NGA.34.23_1", "NGA.35.2_1", "NGA.35.4_1", "NGA.35.8_1", "NGA.36.1_1",
  "NGA.36.2_1", "NGA.36.3_1", "NGA.36.4_1", "NGA.36.5_1", "NGA.36.6_1", "NGA.36.7_1",
  "NGA.36.8_1", "NGA.36.9_1", "NGA.36.10_1", "NGA.36.11_1", "NGA.36.12_1", "NGA.36.13_1",
  "NGA.36.14_1", "NGA.36.15_1", "NGA.36.16_1", "NGA.36.17_1", "NGA.37.1_1", "NGA.37.2_1",
  "NGA.37.3_1", "NGA.37.4_1", "NGA.37.5_1", "NGA.37.6_1", "NGA.37.7_1", "NGA.37.8_1",
  "NGA.37.9_1", "NGA.37.10_1", "NGA.37.11_1", "NGA.37.12_1", "NGA.37.13_1", "NGA.37.14_1"
)

nga_adm2$treatment <- ifelse(nga_adm2$GID_2 %in% ggw_gid2, 1, 0)
cat("  GGW (Treatment=1):", sum(nga_adm2$treatment == 1), "\n")
cat("  GGW (Treatment=0):", sum(nga_adm2$treatment == 0), "\n")

# =============================================================================

# =============================================================================

cat("\n【】\n")

raster_files <- list(
  temp_mean = file.path(cov_dir, "NGA_temp_mean.tif"),
  prec_mean = file.path(cov_dir, "NGA_prec_mean.tif"),
  prec_cv = file.path(cov_dir, "NGA_prec_cv.tif"),
  aridity_mean = file.path(cov_dir, "NGA_aridity_mean.tif"),
  elev_mean = file.path(cov_dir, "NGA_elev_mean.tif"),
  slope_mean = file.path(cov_dir, "NGA_slope_mean.tif"),
  roughness_mean = file.path(cov_dir, "NGA_roughness_mean.tif"),
  ndvi_mean = file.path(cov_dir, "NGA_ndvi_mean_corrected.tif"),
  pop_mean = file.path(cov_dir, "NGA_pop_mean.tif")
)

for(var_name in names(raster_files)) {
  if(file.exists(raster_files[[var_name]])) {
    r <- rast(raster_files[[var_name]])
    nga_adm2[[var_name]] <- exact_extract(r, nga_adm2, 'mean')
    cat("  ✓", var_name, "(mean:", round(mean(nga_adm2[[var_name]], na.rm=TRUE), 3), ")\n")
  }
}

nga_adm2$log_pop <- log(nga_adm2$pop_mean + 1)

# =============================================================================
# ()
# =============================================================================

cat("\n【 ()】\n")

if (!is.null(bio1_trend)) {
  nga_adm2$bio1_trend <- exact_extract(bio1_trend, nga_adm2, 'mean')
  cat("  ✓ bio1_trend (mean:", round(mean(nga_adm2$bio1_trend, na.rm=TRUE), 5), "°C/year)\n")
}

if (!is.null(bio12_trend)) {
  nga_adm2$bio12_trend <- exact_extract(bio12_trend, nga_adm2, 'mean')
  cat("  ✓ bio12_trend (mean:", round(mean(nga_adm2$bio12_trend, na.rm=TRUE), 3), "mm/year)\n")
}

# =============================================================================

# =============================================================================

cat("\n【】\n")

lulc_file <- file.path(cov_dir, "NGA_lulc_2005.tif")

if(file.exists(lulc_file)) {
  lulc <- rast(lulc_file)
  
  lulc_classes <- list(
    forest = c(1, 2, 3, 4, 5), shrubland = c(6, 7), savanna = c(8, 9),
    grassland = c(10), wetland = c(11), cropland = c(12, 14),
    urban = c(13), barren = c(16), water = c(17)
  )
  
  calc_shannon <- function(props) {
    props <- props[props > 0]
    if(length(props) == 0) return(0)
    -sum(props * log(props))
  }
  
  n_adm <- nrow(nga_adm2)
  nga_adm2$lulc_cropland <- NA
  nga_adm2$lulc_diversity <- NA
  nga_vect_lulc <- vect(nga_adm2)
  
  pb <- txtProgressBar(min = 0, max = n_adm, style = 3)
  
  for(i in 1:n_adm) {
    tryCatch({
      adm_i <- nga_vect_lulc[i, ]
      lulc_crop <- crop(lulc, adm_i)
      lulc_mask <- mask(lulc_crop, adm_i)
      vals <- values(lulc_mask, na.rm = TRUE)
      
      if(length(vals) > 0) {
        total_pixels <- length(vals)
        props <- sapply(lulc_classes, function(x) sum(vals %in% x) / total_pixels)
        nga_adm2$lulc_cropland[i] <- props["cropland"]
        nga_adm2$lulc_diversity[i] <- calc_shannon(props)
      }
    }, error = function(e) {})
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat("\n  ✓ lulc_cropland, lulc_diversity\n")
}

# =============================================================================

# =============================================================================

cat("\n【】\n")

richness_files <- list(
  bin_before = "richness_binary_before_2007_2015.tif",
  bin_after = "richness_binary_after_2016_2024.tif",
  bin_change = "change_binary.tif",
  prob_before = "richness_probability_before_2007_2015.tif",
  prob_after = "richness_probability_after_2016_2024.tif",
  prob_change = "change_probability.tif",
  tss_before = "richness_tss_weighted_before_2007_2015.tif",
  tss_after = "richness_tss_weighted_after_2016_2024.tif",
  tss_change = "change_tss_weighted.tif"
)

ref_rast <- NULL
for(name in names(richness_files)) {
  fpath <- file.path(richness_dir, richness_files[[name]])
  if(file.exists(fpath)) {
    if(is.null(ref_rast)) {
      ref_rast <- rast(fpath)
      nga_adm2_rich <- st_transform(nga_adm2, crs(ref_rast))
    }
    r <- rast(fpath)
    nga_adm2[[name]] <- exact_extract(r, nga_adm2_rich, 'mean')
    cat("  ✓", name, "\n")
  }
}


nga_cov <- st_drop_geometry(nga_adm2)

# =============================================================================
# VIF()
# =============================================================================

cat("\n【VIF】\n")

# + 
psm_covars_initial <- c("temp_mean", "prec_mean", "prec_cv", "aridity_mean",
                        "elev_mean", "slope_mean", "roughness_mean",
                        "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity",
                        "bio1_trend", "bio12_trend")

available_covars <- psm_covars_initial[psm_covars_initial %in% names(nga_cov)]
cat("   (", length(available_covars), "):\n")
cat("  ", paste(available_covars, collapse = ", "), "\n")

vif_data <- nga_cov[, c("treatment", available_covars)]
vif_data <- na.omit(vif_data)
cat("  :", nrow(vif_data), "\n")

vif_formula <- as.formula(paste("treatment ~", paste(available_covars, collapse = " + ")))
vif_model <- glm(vif_formula, data = vif_data, family = binomial)
vif_values <- car::vif(vif_model)

cat("\n  VIF:\n")
print(round(sort(vif_values, decreasing = TRUE), 2))

# VIF > 10
final_covars <- available_covars

while(max(vif_values) > 10) {
  high_vif_vars <- names(sort(vif_values, decreasing = TRUE))
  

  remove_var <- NULL
  for(v in high_vif_vars) {
    if(vif_values[v] > 10 && !(v %in% c("bio1_trend", "bio12_trend"))) {
      remove_var <- v
      break
    }
  }
  
  if(is.null(remove_var)) {
    remove_var <- high_vif_vars[1]
  }
  
  cat("  :", remove_var, "(VIF =", round(vif_values[remove_var], 2), ")\n")
  final_covars <- setdiff(final_covars, remove_var)
  
  vif_formula <- as.formula(paste("treatment ~", paste(final_covars, collapse = " + ")))
  vif_model <- glm(vif_formula, data = vif_data, family = binomial)
  vif_values <- car::vif(vif_model)
}

cat("\n   (", length(final_covars), "):\n")
cat("  ", paste(final_covars, collapse = ", "), "\n")

# =============================================================================

# =============================================================================

match_vars <- c("treatment", "GID_2", "NAME_2",
                "bin_change", "prob_change", "tss_change",
                "bin_before", "bin_after", "prob_before", "prob_after",
                "tss_before", "tss_after", final_covars)

nga_match <- nga_cov[, match_vars[match_vars %in% names(nga_cov)]]
nga_match <- na.omit(nga_match)

cat("\n  :", nrow(nga_match), "LGAs\n")
cat("    Treatment=1:", sum(nga_match$treatment == 1), "\n")
cat("    Treatment=0:", sum(nga_match$treatment == 0), "\n")

psm_formula <- as.formula(paste("treatment ~", paste(final_covars, collapse = " + ")))

# =============================================================================
# ATT
# =============================================================================

calc_att <- function(matched_data, final_covars, method_name) {
  treated <- matched_data[matched_data$treatment == 1, ]
  control <- matched_data[matched_data$treatment == 0, ]
  
  results <- data.frame()
  metrics <- c("bin_change", "prob_change", "tss_change")
  labels <- c("Binary-based", "Probability-based", "TSS-weighted")
  
  for(i in 1:3) {
    if(!(metrics[i] %in% names(matched_data))) next
    
    # Weighted t-test
    wt <- wtd.t.test(treated[[metrics[i]]], control[[metrics[i]]],
                     weight = treated$weights, weighty = control$weights, samedata = FALSE)
    att <- wt$additional["Difference"]
    se <- wt$additional["Std. Err"]
    p <- wt$coefficients["p.value"]
    results <- rbind(results, data.frame(
      Richness_Metric = labels[i], Method = "Weighted t-test",
      ATT = att, SE = se, CI_Lower = att-1.96*se, CI_Upper = att+1.96*se, p_value = p))
    
    # WLS + Covariates
    design <- svydesign(ids = ~1, weights = ~weights, data = matched_data)
    adj_covars <- final_covars[1:min(4, length(final_covars))]
    wls <- svyglm(as.formula(paste(metrics[i], "~ treatment +", paste(adj_covars, collapse=" + "))), design = design)
    coef_sum <- summary(wls)$coefficients
    att <- coef_sum["treatment","Estimate"]
    se <- coef_sum["treatment","Std. Error"]
    p <- coef_sum["treatment","Pr(>|t|)"]
    results <- rbind(results, data.frame(
      Richness_Metric = labels[i], Method = "WLS + Covariates",
      ATT = att, SE = se, CI_Lower = att-1.96*se, CI_Upper = att+1.96*se, p_value = p))
    
    # Mixed Effects
    lmer_mod <- lmer(as.formula(paste(metrics[i], "~ treatment + (1|subclass)")),
                     data = matched_data, weights = weights)
    coef_sum <- summary(lmer_mod)$coefficients
    att <- coef_sum["treatment","Estimate"]
    se <- coef_sum["treatment","Std. Error"]
    p <- coef_sum["treatment","Pr(>|t|)"]
    results <- rbind(results, data.frame(
      Richness_Metric = labels[i], Method = "Mixed Effects",
      ATT = att, SE = se, CI_Lower = att-1.96*se, CI_Upper = att+1.96*se, p_value = p))
  }
  
  results$Significance <- ifelse(results$p_value < 0.001, "***",
                                 ifelse(results$p_value < 0.01, "**",
                                        ifelse(results$p_value < 0.05, "*", "")))
  results$Matching <- method_name
  rownames(results) <- NULL
  return(results)
}

# =============================================================================
# Full Matching
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  Full Matching ()                                ║\n")
cat("╚", rep("═", 70), "╝\n\n")

set.seed(123)
psm_full <- matchit(psm_formula, data = nga_match, method = "full",
                    distance = "glm", link = "logit", estimand = "ATT")

print(summary(psm_full))
matched_full <- match.data(psm_full)

# Love plot
png(file.path(output_full, "balance_loveplot.png"), width = 2400, height = 1600, res = 300)
love.plot(psm_full, binary = "std", thresholds = c(m = 0.1), var.order = "unadjusted",
          title = "Full Matching: Covariate Balance (with Climate Trend)",
          colors = c("#E74C3C", "#27AE60"))
dev.off()

# ATT
results_full <- calc_att(matched_full, final_covars, "Full Matching")

cat("\n=== Full Matching ATT ===\n")
print(results_full %>% mutate(ATT = round(ATT, 3), p_value = round(p_value, 4)) %>%
        dplyr::select(Richness_Metric, Method, ATT, p_value, Significance))


write.csv(results_full, file.path(output_full, "ATT_results.csv"), row.names = FALSE)
write.csv(matched_full, file.path(output_full, "matched_data.csv"), row.names = FALSE)

# =============================================================================
# 1:1 Matching
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  1:1 Nearest Neighbor Matching ()                ║\n")
cat("╚", rep("═", 70), "╝\n\n")

set.seed(123)
psm_nn <- matchit(psm_formula, data = nga_match, method = "nearest",
                  distance = "glm", link = "logit", ratio = 1, replace = FALSE, estimand = "ATT")

print(summary(psm_nn))
matched_nn <- match.data(psm_nn)

cat("\n: Treatment=", sum(matched_nn$treatment == 1),
    ", Control=", sum(matched_nn$treatment == 0), "\n")

# Love plot
png(file.path(output_1to1, "balance_loveplot.png"), width = 2400, height = 1600, res = 300)
love.plot(psm_nn, binary = "std", thresholds = c(m = 0.1), var.order = "unadjusted",
          title = "1:1 Matching: Covariate Balance (with Climate Trend)",
          colors = c("#E74C3C", "#27AE60"))
dev.off()

# ATT
results_nn <- calc_att(matched_nn, final_covars, "1:1 Matching")

cat("\n=== 1:1 Matching ATT ===\n")
print(results_nn %>% mutate(ATT = round(ATT, 3), p_value = round(p_value, 4)) %>%
        dplyr::select(Richness_Metric, Method, ATT, p_value, Significance))


write.csv(results_nn, file.path(output_1to1, "ATT_results.csv"), row.names = FALSE)
write.csv(matched_nn, file.path(output_1to1, "matched_data.csv"), row.names = FALSE)

# =============================================================================

# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  4:  -                       ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# ()
effects_file <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition/matched_LGA_effects.csv"

if(file.exists(effects_file)) {
  
  cat("【】\n")
  scenario_effects <- read.csv(effects_file)
  
  # 1:1
  matched_scenario <- matched_nn %>%
    left_join(scenario_effects %>% 
                dplyr::select(GID_2, binary_climate, binary_vegetation, binary_total,
                              prob_climate, prob_vegetation, prob_total),
              by = "GID_2")
  
  cat("  :", nrow(matched_scenario), "\n")
  

  treat <- matched_scenario %>% filter(treatment == 1)
  ctrl <- matched_scenario %>% filter(treatment == 0)
  

  compare_effect <- function(var_name, treat_vals, ctrl_vals) {
    treat_mean <- mean(treat_vals, na.rm = TRUE)
    ctrl_mean <- mean(ctrl_vals, na.rm = TRUE)
    diff <- treat_mean - ctrl_mean
    pooled_sd <- sqrt((var(treat_vals, na.rm=TRUE) + var(ctrl_vals, na.rm=TRUE)) / 2)
    smd <- if(pooled_sd > 0) diff / pooled_sd else NA
    
    t_result <- t.test(treat_vals, ctrl_vals)
    
    status <- ifelse(is.na(smd), "NA",
                     ifelse(abs(smd) < 0.1, "✓ ",
                            ifelse(abs(smd) < 0.25, "~ ", "✗ ")))
    
    return(data.frame(
      Effect = var_name,
      GGW_Mean = treat_mean,
      NonGGW_Mean = ctrl_mean,
      Difference = diff,
      SMD = smd,
      p_value = t_result$p.value,
      Status = status
    ))
  }
  
  cat("\n【 (1:1 Matching)】\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  
  effect_results <- list()
  

  if("binary_climate" %in% names(matched_scenario)) {
    effect_results[[1]] <- compare_effect("Climate Effect (Binary)", 
                                          treat$binary_climate, ctrl$binary_climate)
    cat(sprintf("  %-25s: GGW=%7.3f, Non-GGW=%7.3f, SMD=%6.3f %s\n",
                "Climate Effect (Binary)",
                effect_results[[1]]$GGW_Mean,
                effect_results[[1]]$NonGGW_Mean,
                effect_results[[1]]$SMD,
                effect_results[[1]]$Status))
  }
  
  if("prob_climate" %in% names(matched_scenario)) {
    effect_results[[2]] <- compare_effect("Climate Effect (Prob)", 
                                          treat$prob_climate, ctrl$prob_climate)
    cat(sprintf("  %-25s: GGW=%7.3f, Non-GGW=%7.3f, SMD=%6.3f %s\n",
                "Climate Effect (Prob)",
                effect_results[[2]]$GGW_Mean,
                effect_results[[2]]$NonGGW_Mean,
                effect_results[[2]]$SMD,
                effect_results[[2]]$Status))
  }
  

  if("binary_vegetation" %in% names(matched_scenario)) {
    effect_results[[3]] <- compare_effect("Vegetation Effect (Binary)", 
                                          treat$binary_vegetation, ctrl$binary_vegetation)
    cat(sprintf("  %-25s: GGW=%7.3f, Non-GGW=%7.3f, SMD=%6.3f %s\n",
                "Vegetation Effect (Binary)",
                effect_results[[3]]$GGW_Mean,
                effect_results[[3]]$NonGGW_Mean,
                effect_results[[3]]$SMD,
                effect_results[[3]]$Status))
  }
  
  if("prob_vegetation" %in% names(matched_scenario)) {
    effect_results[[4]] <- compare_effect("Vegetation Effect (Prob)", 
                                          treat$prob_vegetation, ctrl$prob_vegetation)
    cat(sprintf("  %-25s: GGW=%7.3f, Non-GGW=%7.3f, SMD=%6.3f %s\n",
                "Vegetation Effect (Prob)",
                effect_results[[4]]$GGW_Mean,
                effect_results[[4]]$NonGGW_Mean,
                effect_results[[4]]$SMD,
                effect_results[[4]]$Status))
  }
  
  effect_df <- do.call(rbind, effect_results)
  

  write.csv(effect_df, file.path(output_scenario, "effect_comparison_1to1.csv"), row.names = FALSE)
  write.csv(matched_scenario, file.path(output_scenario, "matched_scenario_data_1to1.csv"), row.names = FALSE)
  
  # =============================================================================

  # =============================================================================
  
  cat("\n")
  cat("╔", rep("═", 70), "╗\n", sep = "")
  cat("║                                         ║\n")
  cat("╚", rep("═", 70), "╝\n\n")
  
  climate_effects <- effect_df[grepl("Climate", effect_df$Effect), ]
  
  if(nrow(climate_effects) > 0) {
    all_balanced <- all(abs(climate_effects$SMD) < 0.25, na.rm = TRUE)
    
    if(all_balanced) {
      cat("  ✓  (|SMD| < 0.25)\n\n")
      cat("  :\n")
      cat("    →  (bio1_trend, bio12_trend) \n")
      cat("    → TreatmentControl\n")
      cat("    → GGW\n")
    } else {
      cat("  ⚠️  (|SMD| >= 0.25)\n\n")
      cat("  :\n")
      cat("    → \n")
      cat("    → caliper\n")
    }
    
    cat("\n  :\n")
    print(climate_effects %>% 
            mutate(across(c(GGW_Mean, NonGGW_Mean, Difference, SMD), ~round(., 3))))
  }
  

  if(nrow(effect_df) > 0) {
    plot_data <- effect_df %>%
      dplyr::select(Effect, GGW_Mean, NonGGW_Mean) %>%
      pivot_longer(cols = c(GGW_Mean, NonGGW_Mean),
                   names_to = "Group", values_to = "Effect_Size") %>%
      mutate(Group = ifelse(Group == "GGW_Mean", "GGW", "Non-GGW"),
             Effect_Type = ifelse(grepl("Climate", Effect), "Climate", "Vegetation"))
    
    fig_effects <- ggplot(plot_data, aes(x = Effect, y = Effect_Size, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      scale_fill_manual(values = c("GGW" = "#3B9AB2", "Non-GGW" = "#E78AC3")) +
      labs(title = "Scenario Decomposition: Effect Comparison",
           subtitle = "1:1 Matched LGAs (with Climate Trend Covariates)",
           y = "Effect on Species Richness", x = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "bottom")
    
    ggsave(file.path(output_scenario, "effect_comparison_plot.png"), fig_effects,
           width = 12, height = 6, dpi = 300, bg = "white")
    cat("\n  ✓ : effect_comparison_plot.png\n")
  }
  
} else {
  cat("  ⚠️ :", effects_file, "\n")
  cat("  \n")
}

# =============================================================================

# =============================================================================

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────\n")

check_balance <- function(var_name, treat_vals, ctrl_vals) {
  treat_mean <- mean(treat_vals, na.rm = TRUE)
  ctrl_mean <- mean(ctrl_vals, na.rm = TRUE)
  diff <- treat_mean - ctrl_mean
  pooled_sd <- sqrt((var(treat_vals, na.rm=TRUE) + var(ctrl_vals, na.rm=TRUE)) / 2)
  smd <- if(pooled_sd > 0) diff / pooled_sd else NA
  
  t_result <- t.test(treat_vals, ctrl_vals)
  
  status <- ifelse(is.na(smd), "NA",
                   ifelse(abs(smd) < 0.1, "✓ ",
                          ifelse(abs(smd) < 0.25, "~ ", "✗ ")))
  
  cat(sprintf("  %-15s: GGW=%9.5f, Non-GGW=%9.5f, SMD=%6.3f %s\n",
              var_name, treat_mean, ctrl_mean, smd, status))
  
  return(data.frame(Variable = var_name, GGW_Mean = treat_mean, NonGGW_Mean = ctrl_mean,
                    Diff = diff, SMD = smd, p_value = t_result$p.value, Status = status))
}

treat_nn <- matched_nn %>% filter(treatment == 1)
ctrl_nn <- matched_nn %>% filter(treatment == 0)

balance_results <- list()

cat("\n  === 1:1 Matching ===\n")
if("bio1_trend" %in% names(matched_nn)) {
  balance_results[[1]] <- check_balance("bio1_trend", treat_nn$bio1_trend, ctrl_nn$bio1_trend)
}
if("bio12_trend" %in% names(matched_nn)) {
  balance_results[[2]] <- check_balance("bio12_trend", treat_nn$bio12_trend, ctrl_nn$bio12_trend)
}

treat_full <- matched_full %>% filter(treatment == 1)
ctrl_full <- matched_full %>% filter(treatment == 0)

cat("\n  === Full Matching ===\n")
if("bio1_trend" %in% names(matched_full)) {
  check_balance("bio1_trend", treat_full$bio1_trend, ctrl_full$bio1_trend)
}
if("bio12_trend" %in% names(matched_full)) {
  check_balance("bio12_trend", treat_full$bio12_trend, ctrl_full$bio12_trend)
}

balance_df <- do.call(rbind, balance_results)
write.csv(balance_df, file.path(output_1to1, "climate_trend_balance.csv"), row.names = FALSE)

# =============================================================================

# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  ！                                                          ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# ATT
all_results <- rbind(results_full, results_nn)
write.csv(all_results, file.path(output_base, "ATT_results_combined_with_trend.csv"), row.names = FALSE)

cat("【ATT】\n")
cat("─────────────────────────────────────────────────────────────────\n")
summary_table <- all_results %>%
  mutate(ATT = round(ATT, 3),
         CI = paste0("[", round(CI_Lower,2), ", ", round(CI_Upper,2), "]"),
         p = ifelse(p_value < 0.001, "<0.001", round(p_value, 4))) %>%
  dplyr::select(Matching, Richness_Metric, Method, ATT, CI, p, Significance)
print(summary_table)

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  Full Matching:", nrow(matched_full), "\n")
cat("  1:1 Matching:", sum(matched_nn$treatment == 1), "\n")

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  :", trend_dir, "\n")
cat("  Full Matching:", output_full, "\n")
cat("  1:1 Matching:", output_1to1, "\n")
cat("  :", output_scenario, "\n")

cat("\n=== ！===\n")





# =============================================================================
# PSM+ 

# PSM

# =============================================================================
# =============================================================================
# PSM+ ()
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  :                                    ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# =============================================================================

# =============================================================================

matched_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"
effects_file <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition/matched_LGA_effects.csv"

matched_nn <- read.csv(matched_file)
scenario_effects <- read.csv(effects_file)

matched_data <- matched_nn %>%
  left_join(scenario_effects %>% 
              dplyr::select(GID_2, binary_climate, binary_vegetation, binary_total,
                            prob_climate, prob_vegetation, prob_total),
            by = "GID_2")

output_dir <- "E:/2026.1.8_biomod2/nigeria_results/doubly_robust_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# ATT
# =============================================================================

metrics <- c("bin_change", "prob_change", "tss_change")
labels <- c("Binary", "Probability", "TSS-weighted")

treat <- matched_data %>% filter(treatment == 1)
ctrl <- matched_data %>% filter(treatment == 0)


all_results <- data.frame()

for(i in 1:3) {
  # 1: t
  t_result <- t.test(treat[[metrics[i]]], ctrl[[metrics[i]]])
  simple_att <- t_result$estimate[1] - t_result$estimate[2]
  simple_ci <- t_result$conf.int
  simple_p <- t_result$p.value
  
  all_results <- rbind(all_results, data.frame(
    Metric = labels[i],
    Method = "1. Simple (no climate control)",
    ATT = simple_att,
    CI_Lower = simple_ci[1],
    CI_Upper = simple_ci[2],
    p_value = simple_p
  ))
  
  # 2: 
  formula_adj <- as.formula(paste(metrics[i], "~ treatment + bio1_trend + bio12_trend"))
  model_adj <- lm(formula_adj, data = matched_data, weights = weights)
  coef_adj <- summary(model_adj)$coefficients
  
  all_results <- rbind(all_results, data.frame(
    Metric = labels[i],
    Method = "2. Climate-adjusted",
    ATT = coef_adj["treatment", "Estimate"],
    CI_Lower = coef_adj["treatment", "Estimate"] - 1.96 * coef_adj["treatment", "Std. Error"],
    CI_Upper = coef_adj["treatment", "Estimate"] + 1.96 * coef_adj["treatment", "Std. Error"],
    p_value = coef_adj["treatment", "Pr(>|t|)"]
  ))
  
  # 3: 
  covars <- c("temp_mean", "prec_mean", "ndvi_mean", "log_pop", "bio1_trend", "bio12_trend")
  formula_full <- as.formula(paste(metrics[i], "~ treatment +", paste(covars, collapse = " + ")))
  model_full <- lm(formula_full, data = matched_data, weights = weights)
  coef_full <- summary(model_full)$coefficients
  
  all_results <- rbind(all_results, data.frame(
    Metric = labels[i],
    Method = "3. Fully-adjusted",
    ATT = coef_full["treatment", "Estimate"],
    CI_Lower = coef_full["treatment", "Estimate"] - 1.96 * coef_full["treatment", "Std. Error"],
    CI_Upper = coef_full["treatment", "Estimate"] + 1.96 * coef_full["treatment", "Std. Error"],
    p_value = coef_full["treatment", "Pr(>|t|)"]
  ))
}

# =============================================================================

# =============================================================================

cat("【ATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-15s %-30s %10s %20s %10s\n", "Metric", "Method", "ATT", "95% CI", "p-value"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(all_results)) {
  sig <- ifelse(all_results$p_value[i] < 0.001, "***",
                ifelse(all_results$p_value[i] < 0.01, "**",
                       ifelse(all_results$p_value[i] < 0.05, "*", "")))
  
  ci_str <- sprintf("[%6.2f, %6.2f]", all_results$CI_Lower[i], all_results$CI_Upper[i])
  p_str <- ifelse(all_results$p_value[i] < 0.001, "<0.001", sprintf("%.4f", all_results$p_value[i]))
  
  cat(sprintf("%-15s %-30s %8.2f %20s %10s %s\n",
              all_results$Metric[i],
              all_results$Method[i],
              all_results$ATT[i],
              ci_str,
              p_str,
              sig))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# ATT
# =============================================================================

cat("\n【ATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

comparison <- all_results %>%
  pivot_wider(id_cols = Metric, names_from = Method, values_from = ATT)

names(comparison) <- c("Metric", "Simple", "Climate_Adj", "Full_Adj")

comparison <- comparison %>%
  mutate(
    Change_Climate = round((Climate_Adj - Simple) / abs(Simple) * 100, 1),
    Change_Full = round((Full_Adj - Simple) / abs(Simple) * 100, 1)
  )

cat(sprintf("%-15s %10s %12s %10s %15s %15s\n", 
            "Metric", "Simple", "Climate-Adj", "Full-Adj", "Change(Clim%)", "Change(Full%)"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(comparison)) {
  cat(sprintf("%-15s %10.2f %12.2f %10.2f %15.1f%% %14.1f%%\n",
              comparison$Metric[i],
              comparison$Simple[i],
              comparison$Climate_Adj[i],
              comparison$Full_Adj[i],
              comparison$Change_Climate[i],
              comparison$Change_Full[i]))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

# =============================================================================

# =============================================================================

cat("\n【 vs 】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

effect_analysis <- function(effect_var, label) {

  model_simple <- lm(as.formula(paste(effect_var, "~ treatment")), 
                     data = matched_data, weights = weights)
  

  model_adj <- lm(as.formula(paste(effect_var, "~ treatment + bio1_trend + bio12_trend")), 
                  data = matched_data, weights = weights)
  
  simple_coef <- summary(model_simple)$coefficients["treatment", ]
  adj_coef <- summary(model_adj)$coefficients["treatment", ]
  
  return(data.frame(
    Effect = label,
    Simple_Diff = simple_coef["Estimate"],
    Simple_SE = simple_coef["Std. Error"],
    Simple_p = simple_coef["Pr(>|t|)"],
    Adjusted_Diff = adj_coef["Estimate"],
    Adjusted_SE = adj_coef["Std. Error"],
    Adjusted_p = adj_coef["Pr(>|t|)"]
  ))
}

effect_df <- rbind(
  effect_analysis("binary_climate", "Climate Effect (Binary)"),
  effect_analysis("binary_vegetation", "Vegetation Effect (Binary)"),
  effect_analysis("prob_climate", "Climate Effect (Prob)"),
  effect_analysis("prob_vegetation", "Vegetation Effect (Prob)")
)

cat(sprintf("%-28s %12s %12s %12s\n", "Effect", "Simple", "Adjusted", "Interpretation"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(effect_df)) {
  sig_simple <- ifelse(effect_df$Simple_p[i] < 0.05, "*", "")
  sig_adj <- ifelse(effect_df$Adjusted_p[i] < 0.05, "*", "")
  
  interp <- ""
  if(grepl("Climate", effect_df$Effect[i])) {
    if(effect_df$Simple_p[i] < 0.05 && effect_df$Adjusted_p[i] >= 0.05) {
      interp <- "→ Controlled!"
    }
  } else {
    if(effect_df$Adjusted_p[i] < 0.05) {
      interp <- "→ GGW Policy Effect"
    }
  }
  
  cat(sprintf("%-28s %10.2f%s %10.2f%s %s\n",
              effect_df$Effect[i],
              effect_df$Simple_Diff[i], sig_simple,
              effect_df$Adjusted_Diff[i], sig_adj,
              interp))
}
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("* p < 0.05\n")

# =============================================================================

# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                                                          ║\n")
cat("╚", rep("═", 70), "╝\n\n")

avg_change <- mean(abs(comparison$Change_Climate))

cat("  1. ATT", round(avg_change), "%\n")
cat("     → ATT\n")
cat("     → ATT\n")

cat("\n  2. :\n")
cat("     → Climate Effect (Binary): 10.68* → 1.90 (ns)\n")
cat("     → Climate Effect (Prob): 6.01* → 0.08 (ns)\n")
cat("     → \n")

cat("\n  3.  (GGW):\n")
cat("     → Vegetation Effect (Binary): 4.74* → 7.02*\n")
cat("     → Vegetation Effect (Prob): 1.64* → 2.28*\n")
cat("     → ！\n")

cat("\n  4. ATT:\n")
cat("     → Binary: 7.83 species (Full-adjusted: 10.32)\n")
cat("     → Probability: 1.96 species (3.13)\n")
cat("     → TSS-weighted: 0.79 species (1.25)\n")

# =============================================================================

# =============================================================================

# ATT
all_results$Method <- factor(all_results$Method, 
                             levels = c("3. Fully-adjusted", "2. Climate-adjusted", "1. Simple (no climate control)"))
all_results$Metric <- factor(all_results$Metric, levels = c("Binary", "Probability", "TSS-weighted"))

fig_att <- ggplot(all_results, aes(x = ATT, y = Method, color = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, linewidth = 1) +
  geom_point(size = 3.5) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("#E74C3C", "#3498DB", "#27AE60")) +
  labs(title = "Treatment Effect Sensitivity Analysis",
       subtitle = "Comparing ATT estimates with different levels of climate trend adjustment",
       x = "Average Treatment Effect (ATT)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "gray90"))

ggsave(file.path(output_dir, "ATT_sensitivity_analysis.png"), fig_att,
       width = 14, height = 5, dpi = 300, bg = "white")


effect_plot_data <- effect_df %>%
  dplyr::select(Effect, Simple_Diff, Adjusted_Diff) %>%
  pivot_longer(cols = c(Simple_Diff, Adjusted_Diff),
               names_to = "Model", values_to = "Difference") %>%
  mutate(Model = ifelse(Model == "Simple_Diff", "Simple", "Climate-Adjusted"),
         Effect_Type = ifelse(grepl("Climate", Effect), "Climate Effect", "Vegetation Effect"))

fig_effects <- ggplot(effect_plot_data, aes(x = Effect, y = Difference, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Simple" = "#E74C3C", "Climate-Adjusted" = "#27AE60")) +
  labs(title = "Scenario Decomposition: Climate vs Vegetation Effects",
       subtitle = "Comparing GGW - Non-GGW differences before and after climate trend adjustment",
       y = "Effect Difference (GGW - Non-GGW)", x = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(output_dir, "effect_decomposition_comparison.png"), fig_effects,
       width = 10, height = 6, dpi = 300, bg = "white")

# =============================================================================

# =============================================================================

write.csv(all_results, file.path(output_dir, "ATT_all_methods.csv"), row.names = FALSE)
write.csv(comparison, file.path(output_dir, "ATT_comparison_summary.csv"), row.names = FALSE)
write.csv(effect_df, file.path(output_dir, "effect_decomposition_results.csv"), row.names = FALSE)

cat("\n【】\n")
cat("  ", output_dir, "\n")
cat("    - ATT_all_methods.csv\n")
cat("    - ATT_comparison_summary.csv\n")
cat("    - effect_decomposition_results.csv\n")
cat("    - ATT_sensitivity_analysis.png\n")
cat("    - effect_decomposition_comparison.png\n")

cat("\n=== ！===\n")







# =============================================================================
# PSM

# PSM1: PSM (Δbio) - E:/2026.1.8_biomod2/nigeria_results/PSM_1to1/
# PSM2: ΔbioPSM - E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  PSM                                        ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# =============================================================================

# =============================================================================

# PSM(Δbio) - 2026/1/12
original_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1/Nigeria_1to1_matched_data.csv"

# ΔbioPSM
trend_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"


effects_file <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition/matched_LGA_effects.csv"


output_dir <- "E:/2026.1.8_biomod2/nigeria_results/PSM_comparison_doubly_robust"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================

# =============================================================================

cat("【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

# PSM
if(file.exists(original_psm_file)) {
  original_psm <- read.csv(original_psm_file)
  cat("  PSM: ✓ ,", nrow(original_psm), "\n")
} else {
  stop("PSM: ", original_psm_file)
}

# ΔbioPSM
if(file.exists(trend_psm_file)) {
  trend_psm <- read.csv(trend_psm_file)
  cat("  ΔbioPSM: ✓ ,", nrow(trend_psm), "\n")
} else {
  stop("ΔbioPSM: ", trend_psm_file)
}


if(file.exists(effects_file)) {
  scenario_effects <- read.csv(effects_file)
  cat("  : ✓ ,", nrow(scenario_effects), "\n")
} else {
  cat("  : ✗ \n")
  scenario_effects <- NULL
}

# =============================================================================
# PSMΔbio
# =============================================================================

cat("\n【PSMΔbio】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

# ΔbioPSMΔbio
delta_bio_data <- trend_psm %>%
  dplyr::select(GID_2, bio1_trend, bio12_trend) %>%
  distinct()

# PSM
original_psm <- original_psm %>%
  left_join(delta_bio_data, by = "GID_2")

cat("  PSMbio1_trend:", "bio1_trend" %in% names(original_psm), "\n")
cat("  PSMbio12_trend:", "bio12_trend" %in% names(original_psm), "\n")
cat("  :", sum(!is.na(original_psm$bio1_trend)), "/", nrow(original_psm), "\n")


if(!is.null(scenario_effects)) {
  original_psm <- original_psm %>%
    left_join(scenario_effects %>% 
                dplyr::select(GID_2, binary_climate, binary_vegetation, 
                              prob_climate, prob_vegetation),
              by = "GID_2")
  
  trend_psm <- trend_psm %>%
    left_join(scenario_effects %>% 
                dplyr::select(GID_2, binary_climate, binary_vegetation, 
                              prob_climate, prob_vegetation),
              by = "GID_2")
}

# =============================================================================

# =============================================================================

analyze_psm_doubly_robust <- function(matched_data, psm_name) {
  
  cat("\n")
  cat("╔", rep("═", 68), "╗\n", sep = "")
  cat("║  ", psm_name, rep(" ", max(0, 68 - nchar(psm_name) - 4)), "║\n", sep = "")
  cat("╚", rep("═", 68), "╝\n\n")
  

  cat("  :", nrow(matched_data), "\n")
  cat("    Treatment:", sum(matched_data$treatment == 1), "\n")
  cat("    Control:", sum(matched_data$treatment == 0), "\n")
  
  # Δbio
  treat <- matched_data %>% filter(treatment == 1)
  ctrl <- matched_data %>% filter(treatment == 0)
  
  calc_smd <- function(t_vals, c_vals) {
    diff <- mean(t_vals, na.rm=TRUE) - mean(c_vals, na.rm=TRUE)
    pooled_sd <- sqrt((var(t_vals, na.rm=TRUE) + var(c_vals, na.rm=TRUE)) / 2)
    return(diff / pooled_sd)
  }
  
  smd_bio1 <- calc_smd(treat$bio1_trend, ctrl$bio1_trend)
  smd_bio12 <- calc_smd(treat$bio12_trend, ctrl$bio12_trend)
  
  cat("\n  【Δbio】\n")
  cat(sprintf("    bio1_trend:  GGW=%.5f, Non-GGW=%.5f, SMD=%.3f %s\n", 
              mean(treat$bio1_trend, na.rm=TRUE),
              mean(ctrl$bio1_trend, na.rm=TRUE),
              smd_bio1, 
              ifelse(abs(smd_bio1) < 0.25, "✓", "✗")))
  cat(sprintf("    bio12_trend: GGW=%.3f, Non-GGW=%.3f, SMD=%.3f %s\n", 
              mean(treat$bio12_trend, na.rm=TRUE),
              mean(ctrl$bio12_trend, na.rm=TRUE),
              smd_bio12,
              ifelse(abs(smd_bio12) < 0.25, "✓", "✗")))
  
  # ==========================================================================
  # ATT
  # ==========================================================================
  
  metrics <- c("bin_change", "prob_change", "tss_change")
  labels <- c("Binary", "Probability", "TSS-weighted")
  
  results <- data.frame()
  
  cat("\n  【ATT】\n")
  cat("  ─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  %-15s %-25s %10s %10s\n", "Metric", "Method", "ATT", "p-value"))
  cat("  ─────────────────────────────────────────────────────────────────\n")
  
  for(i in 1:3) {
    if(!(metrics[i] %in% names(matched_data))) next
    
    # 1: ATT
    t_result <- t.test(treat[[metrics[i]]], ctrl[[metrics[i]]])
    simple_att <- t_result$estimate[1] - t_result$estimate[2]
    simple_p <- t_result$p.value
    
    results <- rbind(results, data.frame(
      PSM = psm_name,
      Metric = labels[i],
      Method = "1. Simple",
      ATT = simple_att,
      SE = NA,
      CI_Lower = t_result$conf.int[1],
      CI_Upper = t_result$conf.int[2],
      p_value = simple_p
    ))
    
    sig <- ifelse(simple_p < 0.001, "***", ifelse(simple_p < 0.01, "**", ifelse(simple_p < 0.05, "*", "")))
    cat(sprintf("  %-15s %-25s %10.2f %10s %s\n", labels[i], "1. Simple", simple_att, 
                ifelse(simple_p < 0.001, "<0.001", sprintf("%.4f", simple_p)), sig))
    
    # 2: 
    formula_adj <- as.formula(paste(metrics[i], "~ treatment + bio1_trend + bio12_trend"))
    model_adj <- lm(formula_adj, data = matched_data, weights = weights)
    coef_adj <- summary(model_adj)$coefficients
    
    adj_att <- coef_adj["treatment", "Estimate"]
    adj_se <- coef_adj["treatment", "Std. Error"]
    adj_p <- coef_adj["treatment", "Pr(>|t|)"]
    
    results <- rbind(results, data.frame(
      PSM = psm_name,
      Metric = labels[i],
      Method = "2. Climate-adjusted",
      ATT = adj_att,
      SE = adj_se,
      CI_Lower = adj_att - 1.96 * adj_se,
      CI_Upper = adj_att + 1.96 * adj_se,
      p_value = adj_p
    ))
    
    sig <- ifelse(adj_p < 0.001, "***", ifelse(adj_p < 0.01, "**", ifelse(adj_p < 0.05, "*", "")))
    cat(sprintf("  %-15s %-25s %10.2f %10s %s\n", "", "2. Climate-adjusted", adj_att,
                ifelse(adj_p < 0.001, "<0.001", sprintf("%.4f", adj_p)), sig))
    
    # 3: 
    covars <- c("temp_mean", "prec_mean", "ndvi_mean", "log_pop", "bio1_trend", "bio12_trend")
    available_covars <- covars[covars %in% names(matched_data)]
    
    if(length(available_covars) >= 2) {
      formula_full <- as.formula(paste(metrics[i], "~ treatment +", paste(available_covars, collapse = " + ")))
      model_full <- lm(formula_full, data = matched_data, weights = weights)
      coef_full <- summary(model_full)$coefficients
      
      full_att <- coef_full["treatment", "Estimate"]
      full_se <- coef_full["treatment", "Std. Error"]
      full_p <- coef_full["treatment", "Pr(>|t|)"]
      
      results <- rbind(results, data.frame(
        PSM = psm_name,
        Metric = labels[i],
        Method = "3. Fully-adjusted",
        ATT = full_att,
        SE = full_se,
        CI_Lower = full_att - 1.96 * full_se,
        CI_Upper = full_att + 1.96 * full_se,
        p_value = full_p
      ))
      
      sig <- ifelse(full_p < 0.001, "***", ifelse(full_p < 0.01, "**", ifelse(full_p < 0.05, "*", "")))
      cat(sprintf("  %-15s %-25s %10.2f %10s %s\n", "", "3. Fully-adjusted", full_att,
                  ifelse(full_p < 0.001, "<0.001", sprintf("%.4f", full_p)), sig))
    }
    
    cat("  ─────────────────────────────────────────────────────────────────\n")
  }
  
  # ==========================================================================

  # ==========================================================================
  
  effect_results <- NULL
  
  if("binary_climate" %in% names(matched_data)) {
    cat("\n  【】\n")
    cat("  ─────────────────────────────────────────────────────────────────\n")
    cat(sprintf("  %-25s %12s %12s %12s\n", "Effect", "Simple", "Adjusted", "Change"))
    cat("  ─────────────────────────────────────────────────────────────────\n")
    
    effect_vars <- c("binary_climate", "binary_vegetation", "prob_climate", "prob_vegetation")
    effect_labels <- c("Climate Effect (Binary)", "Vegetation Effect (Binary)", 
                       "Climate Effect (Prob)", "Vegetation Effect (Prob)")
    
    effect_results <- data.frame()
    
    for(k in 1:length(effect_vars)) {
      if(!(effect_vars[k] %in% names(matched_data))) next
      

      model_simple <- lm(as.formula(paste(effect_vars[k], "~ treatment")), 
                         data = matched_data, weights = weights)
      simple_coef <- summary(model_simple)$coefficients["treatment", ]
      

      model_adj <- lm(as.formula(paste(effect_vars[k], "~ treatment + bio1_trend + bio12_trend")), 
                      data = matched_data, weights = weights)
      adj_coef <- summary(model_adj)$coefficients["treatment", ]
      
      effect_results <- rbind(effect_results, data.frame(
        PSM = psm_name,
        Effect = effect_labels[k],
        Simple_Diff = simple_coef["Estimate"],
        Simple_p = simple_coef["Pr(>|t|)"],
        Adjusted_Diff = adj_coef["Estimate"],
        Adjusted_p = adj_coef["Pr(>|t|)"]
      ))
      
      sig_simple <- ifelse(simple_coef["Pr(>|t|)"] < 0.05, "*", "")
      sig_adj <- ifelse(adj_coef["Pr(>|t|)"] < 0.05, "*", "")
      change <- adj_coef["Estimate"] - simple_coef["Estimate"]
      
      cat(sprintf("  %-25s %10.2f%s %10.2f%s %10.2f\n",
                  effect_labels[k],
                  simple_coef["Estimate"], sig_simple,
                  adj_coef["Estimate"], sig_adj,
                  change))
    }
    cat("  ─────────────────────────────────────────────────────────────────\n")
    cat("  * p < 0.05\n")
  }
  
  return(list(att_results = results, effect_results = effect_results,
              smd_bio1 = smd_bio1, smd_bio12 = smd_bio12))
}

# =============================================================================

# =============================================================================

# PSM
results_original <- analyze_psm_doubly_robust(original_psm, "Original PSM (Δbio)")

# ΔbioPSM
results_trend <- analyze_psm_doubly_robust(trend_psm, "PSM with Δbio (Δbio)")

# =============================================================================

# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  PSMATT                                             ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# ATT
all_att <- rbind(results_original$att_results, results_trend$att_results)


comparison_wide <- all_att %>%
  dplyr::select(PSM, Metric, Method, ATT) %>%
  pivot_wider(names_from = PSM, values_from = ATT)

cat("【ATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
print(comparison_wide %>% mutate(across(where(is.numeric), ~round(., 2))))


names(comparison_wide)[3:4] <- c("Original", "With_Dbio")
comparison_wide <- comparison_wide %>%
  mutate(Difference = With_Dbio - Original,
         Pct_Change = round(Difference / abs(Original) * 100, 1))

cat("\n【PSMATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-15s %-25s %10s %10s %10s %10s\n", 
            "Metric", "Method", "Original", "With Δbio", "Diff", "Change%"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(comparison_wide)) {
  cat(sprintf("%-15s %-25s %10.2f %10.2f %10.2f %9.1f%%\n",
              comparison_wide$Metric[i],
              comparison_wide$Method[i],
              comparison_wide$Original[i],
              comparison_wide$With_Dbio[i],
              comparison_wide$Difference[i],
              comparison_wide$Pct_Change[i]))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# Δbio
# =============================================================================

cat("\n【Δbio】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %15s %15s\n", "PSM", "bio1_trend SMD", "bio12_trend SMD"))
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %15.3f %15.3f\n", "Original PSM (Δbio)", 
            results_original$smd_bio1, results_original$smd_bio12))
cat(sprintf("%-30s %15.3f %15.3f\n", "PSM with Δbio (Δbio)", 
            results_trend$smd_bio1, results_trend$smd_bio12))
cat("─────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# PSM()

# PSM1: PSM (Δbio) - E:/2026.1.8_biomod2/nigeria_results/PSM_1to1/
# PSM2: ΔbioPSM - E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  PSM                                        ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

# =============================================================================

# =============================================================================

original_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1/Nigeria_1to1_matched_data.csv"
trend_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"
effects_file <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition/matched_LGA_effects.csv"
output_dir <- "E:/2026.1.8_biomod2/nigeria_results/PSM_comparison_doubly_robust"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================

# =============================================================================

cat("【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

if(file.exists(original_psm_file)) {
  original_psm <- read.csv(original_psm_file)
  cat("  PSM: ✓ ,", nrow(original_psm), "\n")
} else {
  stop("PSM: ", original_psm_file)
}

if(file.exists(trend_psm_file)) {
  trend_psm <- read.csv(trend_psm_file)
  cat("  ΔbioPSM: ✓ ,", nrow(trend_psm), "\n")
} else {
  stop("ΔbioPSM: ", trend_psm_file)
}

if(file.exists(effects_file)) {
  scenario_effects <- read.csv(effects_file)
  cat("  : ✓ ,", nrow(scenario_effects), "\n")
} else {
  cat("  : ✗ \n")
  scenario_effects <- NULL
}

# =============================================================================
# PSMΔbio
# =============================================================================

cat("\n【PSMΔbio】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

delta_bio_data <- trend_psm %>%
  dplyr::select(GID_2, bio1_trend, bio12_trend) %>%
  distinct()

original_psm <- original_psm %>%
  left_join(delta_bio_data, by = "GID_2")

cat("  PSMbio1_trend:", "bio1_trend" %in% names(original_psm), "\n")
cat("  PSMbio12_trend:", "bio12_trend" %in% names(original_psm), "\n")
cat("  :", sum(!is.na(original_psm$bio1_trend)), "/", nrow(original_psm), "\n")

if(!is.null(scenario_effects)) {
  original_psm <- original_psm %>%
    left_join(scenario_effects %>% 
                dplyr::select(GID_2, binary_climate, binary_vegetation, 
                              prob_climate, prob_vegetation),
              by = "GID_2")
  
  trend_psm <- trend_psm %>%
    left_join(scenario_effects %>% 
                dplyr::select(GID_2, binary_climate, binary_vegetation, 
                              prob_climate, prob_vegetation),
              by = "GID_2")
}

# =============================================================================

# =============================================================================

analyze_psm_doubly_robust <- function(matched_data, psm_name) {
  
  cat("\n")
  cat("╔", rep("═", 68), "╗\n", sep = "")
  cat("║  ", psm_name, rep(" ", max(0, 68 - nchar(psm_name) - 4)), "║\n", sep = "")
  cat("╚", rep("═", 68), "╝\n\n", sep = "")
  
  cat("  :", nrow(matched_data), "\n")
  cat("    Treatment:", sum(matched_data$treatment == 1), "\n")
  cat("    Control:", sum(matched_data$treatment == 0), "\n")
  
  treat <- matched_data %>% filter(treatment == 1)
  ctrl <- matched_data %>% filter(treatment == 0)
  
  calc_smd <- function(t_vals, c_vals) {
    diff <- mean(t_vals, na.rm=TRUE) - mean(c_vals, na.rm=TRUE)
    pooled_sd <- sqrt((var(t_vals, na.rm=TRUE) + var(c_vals, na.rm=TRUE)) / 2)
    return(diff / pooled_sd)
  }
  
  smd_bio1 <- calc_smd(treat$bio1_trend, ctrl$bio1_trend)
  smd_bio12 <- calc_smd(treat$bio12_trend, ctrl$bio12_trend)
  
  cat("\n  【Δbio】\n")
  cat(sprintf("    bio1_trend:  GGW=%.5f, Non-GGW=%.5f, SMD=%.3f %s\n", 
              mean(treat$bio1_trend, na.rm=TRUE),
              mean(ctrl$bio1_trend, na.rm=TRUE),
              smd_bio1, 
              ifelse(abs(smd_bio1) < 0.25, "✓", "✗")))
  cat(sprintf("    bio12_trend: GGW=%.3f, Non-GGW=%.3f, SMD=%.3f %s\n", 
              mean(treat$bio12_trend, na.rm=TRUE),
              mean(ctrl$bio12_trend, na.rm=TRUE),
              smd_bio12,
              ifelse(abs(smd_bio12) < 0.25, "✓", "✗")))
  
  # ATT
  metrics <- c("bin_change", "prob_change", "tss_change")
  labels <- c("Binary", "Probability", "TSS-weighted")
  
  results <- data.frame()
  
  cat("\n  【ATT】\n")
  cat("  ─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  %-15s %-25s %10s %10s\n", "Metric", "Method", "ATT", "p-value"))
  cat("  ─────────────────────────────────────────────────────────────────\n")
  
  for(i in 1:3) {
    if(!(metrics[i] %in% names(matched_data))) next
    
    # 1: ATT
    t_result <- t.test(treat[[metrics[i]]], ctrl[[metrics[i]]])
    simple_att <- t_result$estimate[1] - t_result$estimate[2]
    simple_p <- t_result$p.value
    
    results <- rbind(results, data.frame(
      PSM = psm_name, Metric = labels[i], Method = "1. Simple",
      ATT = simple_att, SE = NA,
      CI_Lower = t_result$conf.int[1], CI_Upper = t_result$conf.int[2],
      p_value = simple_p
    ))
    
    sig <- ifelse(simple_p < 0.001, "***", ifelse(simple_p < 0.01, "**", ifelse(simple_p < 0.05, "*", "")))
    cat(sprintf("  %-15s %-25s %10.2f %10s %s\n", labels[i], "1. Simple", simple_att, 
                ifelse(simple_p < 0.001, "<0.001", sprintf("%.4f", simple_p)), sig))
    
    # 2: 
    formula_adj <- as.formula(paste(metrics[i], "~ treatment + bio1_trend + bio12_trend"))
    model_adj <- lm(formula_adj, data = matched_data, weights = weights)
    coef_adj <- summary(model_adj)$coefficients
    
    adj_att <- coef_adj["treatment", "Estimate"]
    adj_se <- coef_adj["treatment", "Std. Error"]
    adj_p <- coef_adj["treatment", "Pr(>|t|)"]
    
    results <- rbind(results, data.frame(
      PSM = psm_name, Metric = labels[i], Method = "2. Climate-adjusted",
      ATT = adj_att, SE = adj_se,
      CI_Lower = adj_att - 1.96 * adj_se, CI_Upper = adj_att + 1.96 * adj_se,
      p_value = adj_p
    ))
    
    sig <- ifelse(adj_p < 0.001, "***", ifelse(adj_p < 0.01, "**", ifelse(adj_p < 0.05, "*", "")))
    cat(sprintf("  %-15s %-25s %10.2f %10s %s\n", "", "2. Climate-adjusted", adj_att,
                ifelse(adj_p < 0.001, "<0.001", sprintf("%.4f", adj_p)), sig))
    
    # 3: 
    covars <- c("temp_mean", "prec_mean", "ndvi_mean", "log_pop", "bio1_trend", "bio12_trend")
    available_covars <- covars[covars %in% names(matched_data)]
    
    if(length(available_covars) >= 2) {
      formula_full <- as.formula(paste(metrics[i], "~ treatment +", paste(available_covars, collapse = " + ")))
      model_full <- lm(formula_full, data = matched_data, weights = weights)
      coef_full <- summary(model_full)$coefficients
      
      full_att <- coef_full["treatment", "Estimate"]
      full_se <- coef_full["treatment", "Std. Error"]
      full_p <- coef_full["treatment", "Pr(>|t|)"]
      
      results <- rbind(results, data.frame(
        PSM = psm_name, Metric = labels[i], Method = "3. Fully-adjusted",
        ATT = full_att, SE = full_se,
        CI_Lower = full_att - 1.96 * full_se, CI_Upper = full_att + 1.96 * full_se,
        p_value = full_p
      ))
      
      sig <- ifelse(full_p < 0.001, "***", ifelse(full_p < 0.01, "**", ifelse(full_p < 0.05, "*", "")))
      cat(sprintf("  %-15s %-25s %10.2f %10s %s\n", "", "3. Fully-adjusted", full_att,
                  ifelse(full_p < 0.001, "<0.001", sprintf("%.4f", full_p)), sig))
    }
    cat("  ─────────────────────────────────────────────────────────────────\n")
  }
  

  effect_results <- NULL
  
  if("binary_climate" %in% names(matched_data)) {
    cat("\n  【】\n")
    cat("  ─────────────────────────────────────────────────────────────────\n")
    cat(sprintf("  %-25s %12s %12s %12s\n", "Effect", "Simple", "Adjusted", "Change"))
    cat("  ─────────────────────────────────────────────────────────────────\n")
    
    effect_vars <- c("binary_climate", "binary_vegetation", "prob_climate", "prob_vegetation")
    effect_labels <- c("Climate Effect (Binary)", "Vegetation Effect (Binary)", 
                       "Climate Effect (Prob)", "Vegetation Effect (Prob)")
    
    effect_results <- data.frame()
    
    for(k in 1:length(effect_vars)) {
      if(!(effect_vars[k] %in% names(matched_data))) next
      
      model_simple <- lm(as.formula(paste(effect_vars[k], "~ treatment")), 
                         data = matched_data, weights = weights)
      simple_coef <- summary(model_simple)$coefficients["treatment", ]
      
      model_adj <- lm(as.formula(paste(effect_vars[k], "~ treatment + bio1_trend + bio12_trend")), 
                      data = matched_data, weights = weights)
      adj_coef <- summary(model_adj)$coefficients["treatment", ]
      
      effect_results <- rbind(effect_results, data.frame(
        PSM = psm_name,
        Effect = effect_labels[k],
        Simple_Diff = simple_coef["Estimate"],
        Simple_p = simple_coef["Pr(>|t|)"],
        Adjusted_Diff = adj_coef["Estimate"],
        Adjusted_p = adj_coef["Pr(>|t|)"]
      ))
      
      sig_simple <- ifelse(simple_coef["Pr(>|t|)"] < 0.05, "*", "")
      sig_adj <- ifelse(adj_coef["Pr(>|t|)"] < 0.05, "*", "")
      change <- adj_coef["Estimate"] - simple_coef["Estimate"]
      
      cat(sprintf("  %-25s %10.2f%s %10.2f%s %10.2f\n",
                  effect_labels[k],
                  simple_coef["Estimate"], sig_simple,
                  adj_coef["Estimate"], sig_adj,
                  change))
    }
    cat("  ─────────────────────────────────────────────────────────────────\n")
    cat("  * p < 0.05\n")
  }
  
  return(list(att_results = results, effect_results = effect_results,
              smd_bio1 = smd_bio1, smd_bio12 = smd_bio12))
}

# =============================================================================

# =============================================================================

results_original <- analyze_psm_doubly_robust(original_psm, "Original PSM (Δbio)")
results_trend <- analyze_psm_doubly_robust(trend_psm, "PSM with Δbio (Δbio)")

# =============================================================================

# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  PSMATT                                             ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

all_att <- rbind(results_original$att_results, results_trend$att_results)

comparison_wide <- all_att %>%
  dplyr::select(PSM, Metric, Method, ATT) %>%
  pivot_wider(names_from = PSM, values_from = ATT)

cat("【ATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
print(comparison_wide %>% mutate(across(where(is.numeric), ~round(., 2))))

names(comparison_wide)[3:4] <- c("Original", "With_Dbio")
comparison_wide <- comparison_wide %>%
  mutate(Difference = With_Dbio - Original,
         Pct_Change = round(Difference / abs(Original) * 100, 1))

cat("\n【PSMATT】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-15s %-25s %10s %10s %10s %10s\n", 
            "Metric", "Method", "Original", "With Δbio", "Diff", "Change%"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(comparison_wide)) {
  cat(sprintf("%-15s %-25s %10.2f %10.2f %10.2f %9.1f%%\n",
              comparison_wide$Metric[i],
              comparison_wide$Method[i],
              comparison_wide$Original[i],
              comparison_wide$With_Dbio[i],
              comparison_wide$Difference[i],
              comparison_wide$Pct_Change[i]))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

# Δbio
cat("\n【Δbio】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %15s %15s\n", "PSM", "bio1_trend SMD", "bio12_trend SMD"))
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %15.3f %15.3f\n", "Original PSM (Δbio)", 
            results_original$smd_bio1, results_original$smd_bio12))
cat(sprintf("%-30s %15.3f %15.3f\n", "PSM with Δbio (Δbio)", 
            results_trend$smd_bio1, results_trend$smd_bio12))
cat("─────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# 【】
# =============================================================================

if(!is.null(results_original$effect_results) && !is.null(results_trend$effect_results)) {
  cat("\n【 (Adjusted)】\n")
  cat("─────────────────────────────────────────────────────────────────────────\n")
  
  # 【】tibble
  effects_orig <- results_original$effect_results
  effects_trend <- results_trend$effect_results
  
  cat(sprintf("%-28s %15s %15s %10s\n", "Effect", "Original PSM", "PSM with Δbio", "Diff"))
  cat("─────────────────────────────────────────────────────────────────────────\n")
  
  for(i in 1:nrow(effects_orig)) {
    val_orig <- effects_orig$Adjusted_Diff[i]
    val_trend <- effects_trend$Adjusted_Diff[i]
    p_orig <- effects_orig$Adjusted_p[i]
    p_trend <- effects_trend$Adjusted_p[i]
    
    sig_orig <- ifelse(p_orig < 0.05, "*", "")
    sig_trend <- ifelse(p_trend < 0.05, "*", "")
    diff <- val_trend - val_orig
    
    cat(sprintf("%-28s %13.2f%s %13.2f%s %10.2f\n",
                effects_orig$Effect[i],
                val_orig, sig_orig,
                val_trend, sig_trend,
                diff))
  }
  cat("─────────────────────────────────────────────────────────────────────────\n")
  cat("* p < 0.05\n")
}

# =============================================================================

# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                                                           ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

cat("  1. Δbio:\n")
cat(sprintf("     PSM:    bio1_trend SMD = %.3f, bio12_trend SMD = %.3f\n",
            results_original$smd_bio1, results_original$smd_bio12))
cat(sprintf("     Δbio PSM: bio1_trend SMD = %.3f, bio12_trend SMD = %.3f\n",
            results_trend$smd_bio1, results_trend$smd_bio12))

if(abs(results_trend$smd_bio1) < abs(results_original$smd_bio1)) {
  improvement <- round((abs(results_original$smd_bio1) - abs(results_trend$smd_bio1)) / 
                         abs(results_original$smd_bio1) * 100, 1)
  cat(sprintf("     → ΔbioPSMbio1_trend %.1f%%\n", improvement))
} else {
  cat("     → PSMΔbio\n")
}

cat("\n  2. ATT ():\n")
simple_orig <- comparison_wide %>% filter(Method == "1. Simple", Metric == "Binary") %>% pull(Original)
simple_trend <- comparison_wide %>% filter(Method == "1. Simple", Metric == "Binary") %>% pull(With_Dbio)
cat(sprintf("     PSM:    Binary ATT = %.2f\n", simple_orig))
cat(sprintf("     Δbio PSM: Binary ATT = %.2f\n", simple_trend))

cat("\n  3. ATT ():\n")
adj_orig <- comparison_wide %>% filter(Method == "2. Climate-adjusted", Metric == "Binary") %>% pull(Original)
adj_trend <- comparison_wide %>% filter(Method == "2. Climate-adjusted", Metric == "Binary") %>% pull(With_Dbio)
cat(sprintf("     PSM:    Binary ATT = %.2f\n", adj_orig))
cat(sprintf("     Δbio PSM: Binary ATT = %.2f\n", adj_trend))

cat("\n  4. :\n")
cat("     → PSMATT\n")
cat("     → PSM\n")
cat("     → \n")
cat("     → : PSM +  \n")

# =============================================================================

# =============================================================================

all_att$PSM_short <- ifelse(grepl("Original", all_att$PSM), "Original PSM", "PSM with Δbio")
all_att$PSM_short <- factor(all_att$PSM_short, levels = c("Original PSM", "PSM with Δbio"))
all_att$Method <- factor(all_att$Method, 
                         levels = c("3. Fully-adjusted", "2. Climate-adjusted", "1. Simple"))
all_att$Metric <- factor(all_att$Metric, levels = c("Binary", "Probability", "TSS-weighted"))

fig_comparison <- ggplot(all_att, aes(x = ATT, y = Method, color = PSM_short, shape = PSM_short)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), 
                 height = 0.2, linewidth = 0.8, position = position_dodge(width = 0.5),
                 na.rm = TRUE) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("Original PSM" = "#E74C3C", "PSM with Δbio" = "#3498DB")) +
  labs(title = "ATT Comparison: Original PSM vs PSM with Climate Trends (Δbio)",
       subtitle = "Both PSMs analyzed with doubly robust estimation (PSM + regression adjustment)",
       x = "Average Treatment Effect (ATT)", y = NULL,
       color = "PSM Version", shape = "PSM Version") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "gray90"))

ggsave(file.path(output_dir, "PSM_comparison_ATT.png"), fig_comparison,
       width = 14, height = 6, dpi = 300, bg = "white")

# =============================================================================

# =============================================================================

cat("\n【】\n")
cat("  ", output_dir, "\n")

write.csv(all_att, file.path(output_dir, "ATT_all_results.csv"), row.names = FALSE)
write.csv(comparison_wide, file.path(output_dir, "ATT_comparison_wide.csv"), row.names = FALSE)

if(!is.null(results_original$effect_results)) {
  all_effects <- rbind(results_original$effect_results, results_trend$effect_results)
  write.csv(all_effects, file.path(output_dir, "effect_decomposition_comparison.csv"), row.names = FALSE)
}

cat("    - ATT_all_results.csv\n")
cat("    - ATT_comparison_wide.csv\n")
cat("    - effect_decomposition_comparison.csv\n")
cat("    - PSM_comparison_ATT.png\n")

cat("\n=== ！===\n")



# =============================================================================



# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                        ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

# =============================================================================
# 1. ()
# =============================================================================

# PSM
original_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1/Nigeria_1to1_matched_data.csv"

# ΔbioPSM
trend_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"

# interaction effect
effects_file <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition/matched_LGA_effects.csv"


output_dir <- "E:/2026.1.8_biomod2/nigeria_results/PSM_comparison_doubly_robust"

# =============================================================================
# 2. 
# =============================================================================

cat("【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

original_psm <- read.csv(original_psm_file)
cat("  PSM:", nrow(original_psm), "\n")

trend_psm <- read.csv(trend_psm_file)
cat("  Δbio PSM:", nrow(trend_psm), "\n")

scenario_effects <- read.csv(effects_file)
cat("  :", nrow(scenario_effects), "\n")

# =============================================================================
# 3. 
# =============================================================================

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

# ΔbioPSM
delta_bio_data <- trend_psm %>%
  dplyr::select(GID_2, bio1_trend, bio12_trend) %>%
  distinct()

# PSM
original_psm <- original_psm %>%
  left_join(delta_bio_data, by = "GID_2")

# interaction
original_psm <- original_psm %>%
  left_join(scenario_effects %>% 
              dplyr::select(GID_2, 
                            binary_climate, binary_vegetation, binary_interaction,
                            prob_climate, prob_vegetation, prob_interaction),
            by = "GID_2")

cat("  :\n")
cat("    - bio1_trend, bio12_trend ()\n")
cat("    - binary_climate, binary_vegetation, binary_interaction\n")
cat("    - prob_climate, prob_vegetation, prob_interaction\n")

# =============================================================================
# 4. 
# =============================================================================

analyze_effect_DR <- function(data, effect_var, effect_name) {
  
  treat <- data %>% filter(treatment == 1)
  ctrl <- data %>% filter(treatment == 0)
  
  # (t-test)
  t_result <- t.test(treat[[effect_var]], ctrl[[effect_var]])
  simple_diff <- mean(treat[[effect_var]], na.rm=TRUE) - mean(ctrl[[effect_var]], na.rm=TRUE)
  simple_p <- t_result$p.value
  
  # ()
  formula_adj <- as.formula(paste(effect_var, "~ treatment + bio1_trend + bio12_trend"))
  model_adj <- lm(formula_adj, data = data, weights = weights)
  coef_adj <- summary(model_adj)$coefficients
  
  adj_diff <- coef_adj["treatment", "Estimate"]
  adj_se <- coef_adj["treatment", "Std. Error"]
  adj_p <- coef_adj["treatment", "Pr(>|t|)"]
  
  return(data.frame(
    Effect = effect_name,
    Simple_Diff = simple_diff,
    Simple_p = simple_p,
    Adjusted_Diff = adj_diff,
    Adjusted_SE = adj_se,
    Adjusted_p = adj_p
  ))
}

# =============================================================================
# 5. 
# =============================================================================

cat("\n【】\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Binary
effects_to_analyze <- list(
  c("binary_climate", "Climate Effect (Binary)"),
  c("binary_vegetation", "Vegetation Effect (Binary)"),
  c("binary_interaction", "Interaction Effect (Binary)")
)

results_binary <- do.call(rbind, lapply(effects_to_analyze, function(x) {
  analyze_effect_DR(original_psm, x[1], x[2])
}))

# Probability
effects_prob <- list(
  c("prob_climate", "Climate Effect (Prob)"),
  c("prob_vegetation", "Vegetation Effect (Prob)"),
  c("prob_interaction", "Interaction Effect (Prob)")
)

results_prob <- do.call(rbind, lapply(effects_prob, function(x) {
  analyze_effect_DR(original_psm, x[1], x[2])
}))


all_results <- rbind(results_binary, results_prob)

# =============================================================================
# 6. 
# =============================================================================


add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

cat("【Binary - 】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-28s %12s %12s %12s\n", "", "", "", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:3) {
  row <- results_binary[i, ]
  sig_simple <- add_sig(row$Simple_p)
  sig_adj <- add_sig(row$Adjusted_p)
  

  if(grepl("Climate", row$Effect)) {
    interpretation <- ifelse(row$Adjusted_p > 0.05, "", "")
  } else if(grepl("Vegetation", row$Effect)) {
    interpretation <- ifelse(row$Adjusted_p < 0.05, "GGW", "")
  } else {
    interpretation <- ifelse(row$Adjusted_p < 0.05, "", "")
  }
  
  cat(sprintf("%-28s %10.2f%s %10.2f%s %s\n",
              row$Effect,
              row$Simple_Diff, sig_simple,
              row$Adjusted_Diff, sig_adj,
              interpretation))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

cat("\n【Probability - 】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-28s %12s %12s %12s\n", "", "", "", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:3) {
  row <- results_prob[i, ]
  sig_simple <- add_sig(row$Simple_p)
  sig_adj <- add_sig(row$Adjusted_p)
  
  if(grepl("Climate", row$Effect)) {
    interpretation <- ifelse(row$Adjusted_p > 0.05, "", "")
  } else if(grepl("Vegetation", row$Effect)) {
    interpretation <- ifelse(row$Adjusted_p < 0.05, "GGW", "")
  } else {
    interpretation <- ifelse(row$Adjusted_p < 0.05, "", "")
  }
  
  cat(sprintf("%-28s %10.2f%s %10.2f%s %s\n",
              row$Effect,
              row$Simple_Diff, sig_simple,
              row$Adjusted_Diff, sig_adj,
              interpretation))
}
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n")

# =============================================================================
# 7. 
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                                                           ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")


int_binary <- results_binary %>% filter(grepl("Interaction", Effect))
int_prob <- results_prob %>% filter(grepl("Interaction", Effect))

cat("【】\n\n")

cat("  1. Binary:\n")
cat(sprintf("     : %.2f (p %s) → : %.2f (p %s)\n",
            int_binary$Simple_Diff,
            ifelse(int_binary$Simple_p < 0.001, "< 0.001", sprintf("= %.3f", int_binary$Simple_p)),
            int_binary$Adjusted_Diff,
            ifelse(int_binary$Adjusted_p < 0.001, "< 0.001", sprintf("= %.3f", int_binary$Adjusted_p))))

cat("\n  2. Probability:\n")
cat(sprintf("     : %.2f (p %s) → : %.2f (p %s)\n",
            int_prob$Simple_Diff,
            ifelse(int_prob$Simple_p < 0.001, "< 0.001", sprintf("= %.3f", int_prob$Simple_p)),
            int_prob$Adjusted_Diff,
            ifelse(int_prob$Adjusted_p < 0.001, "< 0.001", sprintf("= %.3f", int_prob$Adjusted_p))))

cat("\n  3. :\n")

if(int_binary$Adjusted_p < 0.05) {
  if(int_binary$Adjusted_Diff < 0) {
    cat("     → \n")
    cat("     → GGW-GGW\n")
    cat("     → : GGW\n")
  } else {
    cat("     → \n")
    cat("     → GGW-\n")
  }
} else {
  cat("     → \n")
  cat("     → -GGWGGW\n")
}

# =============================================================================
# 8. 
# =============================================================================


write.csv(all_results, file.path(output_dir, "effect_decomposition_with_interaction.csv"), 
          row.names = FALSE)

cat("\n【】\n")
cat("  ", file.path(output_dir, "effect_decomposition_with_interaction.csv"), "\n")

cat("\n=== ！===\n")





# =============================================================================
# Binary, Probability, TSS-weighted

# Wang et al. (2025) Nature Communications
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║  TSS-weighted                   ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

# =============================================================================
# 1. ()
# =============================================================================

# PSM
original_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1/Nigeria_1to1_matched_data.csv"

# ΔbioPSM
trend_psm_file <- "E:/2026.1.8_biomod2/nigeria_results/PSM_1to1_with_trend/matched_data.csv"


effects_file <- "E:/2026.1.8_biomod2/nigeria_results/scenario_decomposition/matched_LGA_effects.csv"


output_dir <- "E:/2026.1.8_biomod2/nigeria_results/PSM_comparison_doubly_robust"

# =============================================================================
# 2. 
# =============================================================================

cat("【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

original_psm <- read.csv(original_psm_file)
cat("  PSM:", nrow(original_psm), "\n")

trend_psm <- read.csv(trend_psm_file)
cat("  Δbio PSM:", nrow(trend_psm), "\n")

scenario_effects <- read.csv(effects_file)
cat("  :", nrow(scenario_effects), "\n")

# =============================================================================
# 3. 
# =============================================================================

cat("\n【】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

# ΔbioPSM
delta_bio_data <- trend_psm %>%
  dplyr::select(GID_2, bio1_trend, bio12_trend) %>%
  distinct()

# PSM
original_psm <- original_psm %>%
  left_join(delta_bio_data, by = "GID_2")

# binary, prob, tss
original_psm <- original_psm %>%
  left_join(scenario_effects %>% 
              dplyr::select(GID_2, 
                            binary_climate, binary_vegetation, binary_interaction,
                            prob_climate, prob_vegetation, prob_interaction,
                            tss_climate, tss_vegetation, tss_interaction),
            by = "GID_2")

cat("  :\n")
cat("    - bio1_trend, bio12_trend ()\n")
cat("    - binary_climate, binary_vegetation, binary_interaction\n")
cat("    - prob_climate, prob_vegetation, prob_interaction\n")
cat("    - tss_climate, tss_vegetation, tss_interaction\n")

# =============================================================================
# 4. 
# =============================================================================

analyze_effect_DR <- function(data, effect_var, effect_name) {
  
  treat <- data %>% filter(treatment == 1)
  ctrl <- data %>% filter(treatment == 0)
  
  # (t-test)
  t_result <- t.test(treat[[effect_var]], ctrl[[effect_var]])
  simple_diff <- mean(treat[[effect_var]], na.rm=TRUE) - mean(ctrl[[effect_var]], na.rm=TRUE)
  simple_p <- t_result$p.value
  
  # ()
  formula_adj <- as.formula(paste(effect_var, "~ treatment + bio1_trend + bio12_trend"))
  model_adj <- lm(formula_adj, data = data, weights = weights)
  coef_adj <- summary(model_adj)$coefficients
  
  adj_diff <- coef_adj["treatment", "Estimate"]
  adj_se <- coef_adj["treatment", "Std. Error"]
  adj_p <- coef_adj["treatment", "Pr(>|t|)"]
  
  return(data.frame(
    Effect = effect_name,
    Simple_Diff = simple_diff,
    Simple_p = simple_p,
    Adjusted_Diff = adj_diff,
    Adjusted_SE = adj_se,
    Adjusted_p = adj_p
  ))
}


add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

# =============================================================================
# 5. 
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                                                   ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

# ----- Binary-----
effects_binary <- list(
  c("binary_climate", "Climate Effect"),
  c("binary_vegetation", "Vegetation Effect"),
  c("binary_interaction", "Interaction Effect")
)
results_binary <- do.call(rbind, lapply(effects_binary, function(x) {
  analyze_effect_DR(original_psm, x[1], x[2])
}))
results_binary$Metric <- "Binary"

# ----- Probability-----
effects_prob <- list(
  c("prob_climate", "Climate Effect"),
  c("prob_vegetation", "Vegetation Effect"),
  c("prob_interaction", "Interaction Effect")
)
results_prob <- do.call(rbind, lapply(effects_prob, function(x) {
  analyze_effect_DR(original_psm, x[1], x[2])
}))
results_prob$Metric <- "Probability"

# ----- TSS-weighted-----
effects_tss <- list(
  c("tss_climate", "Climate Effect"),
  c("tss_vegetation", "Vegetation Effect"),
  c("tss_interaction", "Interaction Effect")
)
results_tss <- do.call(rbind, lapply(effects_tss, function(x) {
  analyze_effect_DR(original_psm, x[1], x[2])
}))
results_tss$Metric <- "TSS-weighted"


all_results <- rbind(results_binary, results_prob, results_tss)

# =============================================================================
# 6. 
# =============================================================================


get_interpretation <- function(effect, adj_p, adj_diff) {
  if(grepl("Climate", effect)) {
    return(ifelse(adj_p > 0.05, "", ""))
  } else if(grepl("Vegetation", effect)) {
    return(ifelse(adj_p < 0.05, "GGW", ""))
  } else {
    if(adj_p < 0.05) {
      return(ifelse(adj_diff < 0, "GGW", ""))
    } else {
      return("")
    }
  }
}

# ----- Binary-----
cat("【Binary - 】(Wang et al. 2025)\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-20s %12s %12s %s\n", "", "", "", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:3) {
  row <- results_binary[i, ]
  sig_simple <- add_sig(row$Simple_p)
  sig_adj <- add_sig(row$Adjusted_p)
  interpretation <- get_interpretation(row$Effect, row$Adjusted_p, row$Adjusted_Diff)
  
  cat(sprintf("%-20s %10.2f%s %10.2f%s %s\n",
              row$Effect,
              row$Simple_Diff, sig_simple,
              row$Adjusted_Diff, sig_adj,
              interpretation))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

# ----- Probability-----
cat("\n【Probability - 】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-20s %12s %12s %s\n", "", "", "", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:3) {
  row <- results_prob[i, ]
  sig_simple <- add_sig(row$Simple_p)
  sig_adj <- add_sig(row$Adjusted_p)
  interpretation <- get_interpretation(row$Effect, row$Adjusted_p, row$Adjusted_Diff)
  
  cat(sprintf("%-20s %10.2f%s %10.2f%s %s\n",
              row$Effect,
              row$Simple_Diff, sig_simple,
              row$Adjusted_Diff, sig_adj,
              interpretation))
}
cat("─────────────────────────────────────────────────────────────────────────\n")

# ----- TSS-weighted-----
cat("\n【TSS-weighted - 】()\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-20s %12s %12s %s\n", "", "", "", ""))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:3) {
  row <- results_tss[i, ]
  sig_simple <- add_sig(row$Simple_p)
  sig_adj <- add_sig(row$Adjusted_p)
  interpretation <- get_interpretation(row$Effect, row$Adjusted_p, row$Adjusted_Diff)
  
  cat(sprintf("%-20s %10.2f%s %10.2f%s %s\n",
              row$Effect,
              row$Simple_Diff, sig_simple,
              row$Adjusted_Diff, sig_adj,
              interpretation))
}
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n")

# =============================================================================
# 7. 
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                                                   ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")


summary_table <- all_results %>%
  mutate(
    Simple = sprintf("%.2f%s", Simple_Diff, add_sig(Simple_p)),
    Adjusted = sprintf("%.2f%s", Adjusted_Diff, add_sig(Adjusted_p))
  ) %>%
  dplyr::select(Metric, Effect, Simple, Adjusted)

cat("【Table: Scenario decomposition under doubly robust framework】\n")
cat("─────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-15s %-20s %15s %15s\n", "Metric", "Effect", "Simple", "Adjusted"))
cat("─────────────────────────────────────────────────────────────────────────\n")

for(i in 1:nrow(summary_table)) {
  cat(sprintf("%-15s %-20s %15s %15s\n",
              summary_table$Metric[i],
              summary_table$Effect[i],
              summary_table$Simple[i],
              summary_table$Adjusted[i]))
}
cat("─────────────────────────────────────────────────────────────────────────\n")
cat("Note: Values represent differences between GGW and Non-GGW regions.\n")
cat("Adjusted values control for climate trends (bio1_trend, bio12_trend).\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n")

# =============================================================================
# 8. 
# =============================================================================

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║                                                           ║\n")
cat("╚", rep("═", 70), "╝\n\n", sep = "")

cat("【1. 】\n")
cat(sprintf("   Binary:      %.2f → %.2f (%s)\n",
            results_binary$Simple_Diff[1], results_binary$Adjusted_Diff[1],
            ifelse(results_binary$Adjusted_p[1] > 0.05, " ✓", "")))
cat(sprintf("   Probability: %.2f → %.2f (%s)\n",
            results_prob$Simple_Diff[1], results_prob$Adjusted_Diff[1],
            ifelse(results_prob$Adjusted_p[1] > 0.05, " ✓", "")))
cat(sprintf("   TSS:         %.2f → %.2f (%s)\n",
            results_tss$Simple_Diff[1], results_tss$Adjusted_Diff[1],
            ifelse(results_tss$Adjusted_p[1] > 0.05, " ✓", "")))

cat("\n【2. 】\n")
cat(sprintf("   Binary:      %.2f → %.2f%s (GGW)\n",
            results_binary$Simple_Diff[2], results_binary$Adjusted_Diff[2],
            add_sig(results_binary$Adjusted_p[2])))
cat(sprintf("   Probability: %.2f → %.2f%s (GGW)\n",
            results_prob$Simple_Diff[2], results_prob$Adjusted_Diff[2],
            add_sig(results_prob$Adjusted_p[2])))
cat(sprintf("   TSS:         %.2f → %.2f%s\n",
            results_tss$Simple_Diff[2], results_tss$Adjusted_Diff[2],
            add_sig(results_tss$Adjusted_p[2])))

cat("\n【3. 】\n")
cat(sprintf("   Binary:      %.2f → %.2f%s\n",
            results_binary$Simple_Diff[3], results_binary$Adjusted_Diff[3],
            add_sig(results_binary$Adjusted_p[3])))
cat(sprintf("   Probability: %.2f → %.2f%s\n",
            results_prob$Simple_Diff[3], results_prob$Adjusted_Diff[3],
            add_sig(results_prob$Adjusted_p[3])))
cat(sprintf("   TSS:         %.2f → %.2f%s\n",
            results_tss$Simple_Diff[3], results_tss$Adjusted_Diff[3],
            add_sig(results_tss$Adjusted_p[3])))

cat("\n【4. 】\n")

climate_consistent <- all(c(results_binary$Adjusted_p[1] > 0.05,
                            results_prob$Adjusted_p[1] > 0.05))
veg_consistent <- all(c(results_binary$Adjusted_p[2] < 0.05,
                        results_prob$Adjusted_p[2] < 0.05))

if(climate_consistent) {
  cat("   ✓ : BinaryProbability\n")
} else {
  cat("   ⚠ : \n")
}

if(veg_consistent) {
  cat("   ✓ : BinaryProbabilityGGW\n")
} else {
  cat("   ⚠ : \n")
}

# =============================================================================
# 9. 
# =============================================================================


write.csv(all_results, file.path(output_dir, "scenario_decomposition_all_metrics_DR.csv"), 
          row.names = FALSE)

cat("\n【】\n")
cat("  ", file.path(output_dir, "scenario_decomposition_all_metrics_DR.csv"), "\n")

cat("\n=== ！===\n")