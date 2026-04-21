# =============================================================================
# BIO
# =============================================================================

library(terra)

bio_dir <- "E:/GGW_bird_analysis/ENV/Senegal"

cat("=== BIO ===\n\n")

# BIO_2000
year_dir <- file.path(bio_dir, "BIO_2000")

if(dir.exists(year_dir)) {
  cat("BIO_2000 \n")
  cat(":", year_dir, "\n\n")
  

  all_files <- list.files(year_dir, full.names = FALSE)
  cat(" (", length(all_files), "):\n", sep = "")
  print(all_files)
  
  # tif
  tif_files <- list.files(year_dir, pattern = "\\.tif$", full.names = TRUE)
  cat("\nTIF (", length(tif_files), "):\n", sep = "")
  for(f in tif_files) {
    cat("  -", basename(f), "\n")
  }
  
  # tif
  if(length(tif_files) > 0) {
    cat("\nTIF:\n")
    r <- rast(tif_files[1])
    cat("  :", basename(tif_files[1]), "\n")
    cat("  :", nlyr(r), "\n")
    cat("  :", paste(names(r), collapse = ", "), "\n")
    cat("  :", res(r), "\n")
    cat("  :", as.vector(ext(r)), "\n")
  }
  
} else {
  cat("BIO_2000 !\n")
  cat(":\n")
  print(list.files(bio_dir))
}


# =============================================================================

# =============================================================================

library(terra)

output_dir <- "E:/2026.1.8_biomod2/SEN_covariates"

cat("===  ===\n\n")


temp <- rast(file.path(output_dir, "SEN_temp_mean_2000_2006.tif"))


temp_stats <- global(temp, c("min", "max", "mean"), na.rm = TRUE)
cat(":\n")
cat("  :", round(temp_stats$min, 2), "\n")
cat("  :", round(temp_stats$max, 2), "\n")
cat("  :", round(temp_stats$mean, 2), "\n")

# 10
temp_fixed <- temp * 10


temp_fixed_stats <- global(temp_fixed, c("min", "max", "mean"), na.rm = TRUE)
cat("\n (×10):\n")
cat("  :", round(temp_fixed_stats$min, 2), "°C\n")
cat("  :", round(temp_fixed_stats$max, 2), "°C\n")
cat("  :", round(temp_fixed_stats$mean, 2), "°C\n")


if(temp_fixed_stats$min > 20 && temp_fixed_stats$max < 40) {
  cat("  ✓  (24-30°C)\n")
  

  writeRaster(temp_fixed, file.path(output_dir, "SEN_temp_mean_2000_2006.tif"), overwrite = TRUE)
  cat("\n  ✓ : SEN_temp_mean_2000_2006.tif\n")
  
} else {
  cat("  ⚠ BIO\n")
}

# =============================================================================

# =============================================================================

cat("\n===  ===\n\n")


prec <- rast(file.path(output_dir, "SEN_prec_mean_2000_2006.tif"))

# = / (+ 10)
aridity_fixed <- prec / (temp_fixed + 10)


aridity_stats <- global(aridity_fixed, c("min", "max", "mean"), na.rm = TRUE)
cat(" ():\n")
cat("  :", round(aridity_stats$min, 2), "\n")
cat("  :", round(aridity_stats$max, 2), "\n")
cat("  :", round(aridity_stats$mean, 2), "\n")
cat("  : <20, 20-50, 50-75, >75\n")


writeRaster(aridity_fixed, file.path(output_dir, "SEN_aridity_2000_2006.tif"), overwrite = TRUE)
cat("\n  ✓ : SEN_aridity_2000_2006.tif\n")

# =============================================================================

# =============================================================================

cat("\n===  ===\n\n")


temp_final <- rast(file.path(output_dir, "SEN_temp_mean_2000_2006.tif"))
prec_final <- rast(file.path(output_dir, "SEN_prec_mean_2000_2006.tif"))
cv_final <- rast(file.path(output_dir, "SEN_prec_cv_2000_2006.tif"))
aridity_final <- rast(file.path(output_dir, "SEN_aridity_2000_2006.tif"))

summary_df <- data.frame(
   = c("", "", "CV", ""),
   = c("°C", "mm", "%", "-"),
   = round(c(
    global(temp_final, "min", na.rm=T)$min,
    global(prec_final, "min", na.rm=T)$min,
    global(cv_final, "min", na.rm=T)$min,
    global(aridity_final, "min", na.rm=T)$min
  ), 2),
   = round(c(
    global(temp_final, "max", na.rm=T)$max,
    global(prec_final, "max", na.rm=T)$max,
    global(cv_final, "max", na.rm=T)$max,
    global(aridity_final, "max", na.rm=T)$max
  ), 2),
   = round(c(
    global(temp_final, "mean", na.rm=T)$mean,
    global(prec_final, "mean", na.rm=T)$mean,
    global(cv_final, "mean", na.rm=T)$mean,
    global(aridity_final, "mean", na.rm=T)$mean
  ), 2),
   = c("24-30", "200-1500", "60-150", "5-50")
)

print(summary_df, row.names = FALSE)

cat("\n===  ===\n")
# =============================================================================

# =============================================================================

library(terra)

cat("\n")
cat("╔", rep("═", 65), "╗\n", sep = "")
cat("║                                            ║\n")
cat("╚", rep("═", 65), "╝\n\n")

output_dir <- "E:/2026.1.8_biomod2/SEN_covariates"

# =============================================================================

# =============================================================================

# 1. 
cat("【1】 (temp_mean)\n")
cat(rep("-", 50), "\n", sep = "")
temp <- rast(file.path(output_dir, "SEN_temp_mean_2000_2006.tif"))
temp_stats <- global(temp, c("min", "max", "mean", "sd"), na.rm = TRUE)
cat("  :", round(temp_stats$min, 2), "°C\n")
cat("  :", round(temp_stats$max, 2), "°C\n")
cat("  :", round(temp_stats$mean, 2), "°C\n")
cat("  :", round(temp_stats$sd, 2), "\n")
cat("  :  24-30°C\n")
if(temp_stats$min > 15 && temp_stats$max < 45) {
  cat("  ✓ \n")
} else {
  cat("  ⚠ \n")
}

# 2. 
cat("\n【2】 (prec_mean)\n")
cat(rep("-", 50), "\n", sep = "")
prec <- rast(file.path(output_dir, "SEN_prec_mean_2000_2006.tif"))
prec_stats <- global(prec, c("min", "max", "mean", "sd"), na.rm = TRUE)
cat("  :", round(prec_stats$min, 2), "mm\n")
cat("  :", round(prec_stats$max, 2), "mm\n")
cat("  :", round(prec_stats$mean, 2), "mm\n")
cat("  :", round(prec_stats$sd, 2), "\n")
cat("  :  200-1500mm ()\n")
if(prec_stats$min >= 0 && prec_stats$max < 3000) {
  cat("  ✓ \n")
} else {
  cat("  ⚠ \n")
}

# 3. 
cat("\n【3】 (prec_cv)\n")
cat(rep("-", 50), "\n", sep = "")
cv <- rast(file.path(output_dir, "SEN_prec_cv_2000_2006.tif"))
cv_stats <- global(cv, c("min", "max", "mean", "sd"), na.rm = TRUE)
cat("  :", round(cv_stats$min, 2), "\n")
cat("  :", round(cv_stats$max, 2), "\n")
cat("  :", round(cv_stats$mean, 2), "\n")
cat("  :", round(cv_stats$sd, 2), "\n")
cat("  : CV 60-120%\n")
if(cv_stats$min >= 0 && cv_stats$max < 200) {
  cat("  ✓ \n")
} else {
  cat("  ⚠ \n")
}

# 4. 
cat("\n【4】 (aridity)\n")
cat(rep("-", 50), "\n", sep = "")
aridity <- rast(file.path(output_dir, "SEN_aridity_2000_2006.tif"))
aridity_stats <- global(aridity, c("min", "max", "mean", "sd"), na.rm = TRUE)
cat("  :", round(aridity_stats$min, 2), "\n")
cat("  :", round(aridity_stats$max, 2), "\n")
cat("  :", round(aridity_stats$mean, 2), "\n")
cat("  :", round(aridity_stats$sd, 2), "\n")
cat("  :  <20 , 20-50 , 50-75 , >75 \n")
if(aridity_stats$min >= 0 && aridity_stats$max < 200) {
  cat("  ✓ \n")
} else {
  cat("  ⚠ \n")
}

# =============================================================================

# =============================================================================

cat("\n【5】\n")
cat(rep("-", 50), "\n", sep = "")


png(file.path(output_dir, "SEN_climate_check.png"), width = 2400, height = 2000, res = 300)
par(mfrow = c(2, 2), mar = c(3, 3, 3, 4))


plot(temp, main = " (°C)", col = hcl.colors(50, "YlOrRd"))


plot(prec, main = " (mm)", col = hcl.colors(50, "YlGnBu"))

# CV
plot(cv, main = " (%)", col = hcl.colors(50, "Purples"))


plot(aridity, main = "", col = hcl.colors(50, "RdYlGn"))

dev.off()

cat("  ✓ : SEN_climate_check.png\n")

# =============================================================================

# =============================================================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("\n")
cat(rep("=", 60), "\n", sep = "")

summary_df <- data.frame(
   = c("", "", "CV", ""),
   = c("°C", "mm", "%", "-"),
   = round(c(temp_stats$min, prec_stats$min, cv_stats$min, aridity_stats$min), 2),
   = round(c(temp_stats$max, prec_stats$max, cv_stats$max, aridity_stats$max), 2),
   = round(c(temp_stats$mean, prec_stats$mean, cv_stats$mean, aridity_stats$mean), 2),
   = c("24-30", "200-1500", "60-120", "5-50")
)

print(summary_df, row.names = FALSE)

cat("\n:\n")
cat("  - :  <400mm\n")
cat("  - :  >1000mm\n")
cat("  -  (25-30°C)\n")
cat("  -  (6-10)\n")



# =============================================================================

# Adm4+ lulc_croplandlulc_diversity
# =============================================================================

library(terra)
library(sf)
library(dplyr)

cat("\n")
cat("╔", rep("═", 70), "╗\n", sep = "")
cat("║   | Senegal Covariate Extraction                ║\n")
cat("╚", rep("═", 70), "╝\n\n")

# =============================================================================

# =============================================================================

cov_dir <- "E:/2026.1.8_biomod2/SEN_covariates"
adm_file <- "E:/11.17progress/study_area/SEN_adm4.shp"
ggw_file <- "E:/11.17progress/study_area/SEN_GGW.shp"
output_file <- file.path(cov_dir, "Senegal_Adm4_covariates.csv")

# =============================================================================

# =============================================================================

cat("【 1】\n")
cat(rep("-", 50), "\n", sep = "")


sen_adm4 <- st_read(adm_file, quiet = TRUE)
cat("  :", nrow(sen_adm4), "\n")
cat("  :", paste(names(sen_adm4)[1:5], collapse = ", "), "...\n")

# GGW
sen_ggw <- st_read(ggw_file, quiet = TRUE)
sen_ggw <- st_transform(sen_ggw, st_crs(sen_adm4))
ggw_union <- st_union(sen_ggw)
cat("  GGW: \n")

# =============================================================================
# Treatment
# =============================================================================

cat("\n【 2】Treatment\n")
cat(rep("-", 50), "\n", sep = "")

# GGW
sen_adm4$treatment <- as.integer(st_intersects(sen_adm4, ggw_union, sparse = FALSE))

cat("  GGW (Treatment=1):", sum(sen_adm4$treatment == 1), "\n")
cat("  GGW (Treatment=0):", sum(sen_adm4$treatment == 0), "\n")

# =============================================================================

# =============================================================================

cat("\n【 3】\n")
cat(rep("-", 50), "\n", sep = "")


temp_mean <- rast(file.path(cov_dir, "SEN_temp_mean_2000_2006.tif"))
prec_mean <- rast(file.path(cov_dir, "SEN_prec_mean_2000_2006.tif"))
prec_cv <- rast(file.path(cov_dir, "SEN_prec_cv_2000_2006.tif"))
aridity <- rast(file.path(cov_dir, "SEN_aridity_2000_2006.tif"))
cat("  ✓  (4)\n")


elevation <- rast(file.path(cov_dir, "SEN_elevation.tif"))
slope <- rast(file.path(cov_dir, "SEN_slope.tif"))
roughness <- rast(file.path(cov_dir, "SEN_roughness.tif"))
cat("  ✓  (3)\n")


ndvi <- rast(file.path(cov_dir, "SEN_NDVI_2000_2006_mean.tif"))
pop <- rast(file.path(cov_dir, "SEN_pop_2000_2006_mean.tif"))
cat("  ✓  (2)\n")


lulc <- rast(file.path(cov_dir, "SEN_LULC_2005.tif"))
cat("  ✓  (1)\n")

# =============================================================================

# =============================================================================

cat("\n【 4】\n")
cat(rep("-", 50), "\n", sep = "")

# SpatVector
sen_vect <- vect(sen_adm4)


cat("   temp_mean...")
sen_adm4$temp_mean <- terra::extract(temp_mean, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   prec_mean...")
sen_adm4$prec_mean <- terra::extract(prec_mean, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   prec_cv...")
sen_adm4$prec_cv <- terra::extract(prec_cv, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   aridity_mean...")
sen_adm4$aridity_mean <- terra::extract(aridity, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   elev_mean...")
sen_adm4$elev_mean <- terra::extract(elevation, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   slope_mean...")
sen_adm4$slope_mean <- terra::extract(slope, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   roughness_mean...")
sen_adm4$roughness_mean <- terra::extract(roughness, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   ndvi_mean...")
sen_adm4$ndvi_mean <- terra::extract(ndvi, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

cat("   pop_mean...")
sen_adm4$pop_mean <- terra::extract(pop, sen_vect, fun = mean, na.rm = TRUE)[, 2]
cat(" ✓\n")

# log_pop
sen_adm4$log_pop <- log(sen_adm4$pop_mean + 1)

# =============================================================================
# (lulc_cropland + lulc_diversity)
# =============================================================================

cat("\n【 5】\n")
cat(rep("-", 50), "\n", sep = "")

# MODIS MCD12Q1 LC_Type1 
lulc_classes <- list(
  forest = c(1, 2, 3, 4, 5),
  shrubland = c(6, 7),
  savanna = c(8, 9),
  grassland = c(10),
  wetland = c(11),
  cropland = c(12, 14),
  urban = c(13),
  barren = c(16),
  water = c(17)
)

# Shannon
calc_shannon <- function(props) {
  props <- props[props > 0]
  if(length(props) == 0) return(0)
  -sum(props * log(props))
}


n_adm <- nrow(sen_adm4)
sen_adm4$lulc_cropland <- NA
sen_adm4$lulc_diversity <- NA

cat("  ...\n")
pb <- txtProgressBar(min = 0, max = n_adm, style = 3)

for(i in 1:n_adm) {
  tryCatch({

    adm_i <- sen_vect[i, ]
    lulc_crop <- crop(lulc, adm_i)
    lulc_mask <- mask(lulc_crop, adm_i)
    

    vals <- values(lulc_mask, na.rm = TRUE)
    
    if(length(vals) > 0) {
      total_pixels <- length(vals)
      

      props <- numeric(length(lulc_classes))
      names(props) <- names(lulc_classes)
      
      for(lc_name in names(lulc_classes)) {
        props[lc_name] <- sum(vals %in% lulc_classes[[lc_name]]) / total_pixels
      }
      
      # cropland
      sen_adm4$lulc_cropland[i] <- props["cropland"]
      
      # Shannon
      sen_adm4$lulc_diversity[i] <- calc_shannon(props)
    }
    
  }, error = function(e) {
    # NA
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)
cat("\n  ✓ \n")

# =============================================================================

# =============================================================================

cat("\n【 6】\n")
cat(rep("-", 50), "\n", sep = "")


centroids <- st_coordinates(st_centroid(sen_adm4))
sen_adm4$longitude <- centroids[, 1]
sen_adm4$latitude <- centroids[, 2]


sen_adm4$area_km2 <- as.numeric(st_area(sen_adm4)) / 1e6

cat("  ✓ \n")

# =============================================================================

# =============================================================================

cat("\n【 7】\n")
cat(rep("-", 50), "\n", sep = "")

# geometry
sen_cov <- st_drop_geometry(sen_adm4)



cat("  :\n")
print(names(sen_cov))

# CSV
write.csv(sen_cov, output_file, row.names = FALSE)
cat("\n  ✓ :", output_file, "\n")

# =============================================================================

# =============================================================================

cat("\n【 8】\n")
cat(rep("-", 50), "\n", sep = "")

# PSM
psm_vars <- c("temp_mean", "prec_mean", "prec_cv", "aridity_mean",
              "elev_mean", "slope_mean", "roughness_mean",
              "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")

cat("\n  PSM:\n")
for(v in psm_vars) {
  if(v %in% names(sen_cov)) {
    stats <- c(
      min = min(sen_cov[[v]], na.rm = TRUE),
      max = max(sen_cov[[v]], na.rm = TRUE),
      mean = mean(sen_cov[[v]], na.rm = TRUE),
      na_count = sum(is.na(sen_cov[[v]]))
    )
    cat(sprintf("  %-15s: min=%.2f, max=%.2f, mean=%.2f, NA=%d\n",
                v, stats["min"], stats["max"], stats["mean"], stats["na_count"]))
  } else {
    cat(sprintf("  %-15s: ✗ \n", v))
  }
}

# Treatment
cat("\n  Treatment:\n")
cat("    GGW (1):", sum(sen_cov$treatment == 1, na.rm = TRUE), "\n")
cat("    GGW (0):", sum(sen_cov$treatment == 0, na.rm = TRUE), "\n")

# =============================================================================

# =============================================================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("！\n")
cat(rep("=", 60), "\n", sep = "")

cat("\n:", output_file, "\n")
cat(":", nrow(sen_cov), "\n")
cat(":", ncol(sen_cov), "\n")

cat("\nPSM (11):\n")
cat(paste(psm_vars, collapse = ", "), "\n")

cat("\n: VIF → PSM → DID\n")



# =============================================================================
# NA
# =============================================================================

sen_cov <- read.csv("E:/2026.1.8_biomod2/SEN_covariates/Senegal_Adm4_covariates.csv")

cat("=== NA ===\n\n")

psm_covars <- c("temp_mean", "prec_mean", "prec_cv", "aridity_mean",
                "elev_mean", "slope_mean", "roughness_mean",
                "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")

for(v in psm_covars) {
  na_count <- sum(is.na(sen_cov[[v]]))
  if(na_count > 0) {
    cat(v, ": ", na_count, "NA\n", sep = "")
    # NA
    na_rows <- which(is.na(sen_cov[[v]]))
    cat("  NA:\n")
    print(sen_cov[na_rows, c("GID_4", "NAME_4", "NAME_3", "NAME_2", v)])
  }
}

# NA
cat("\n=== NA ===\n")
na_rows <- which(rowSums(is.na(sen_cov[, psm_covars])) > 0)
cat("", length(na_rows), "NA\n\n")

if(length(na_rows) > 0) {
  print(sen_cov[na_rows, c("GID_4", "NAME_4", "NAME_2", "area_km2", "longitude", "latitude")])
  

  cat("\n:\n")
  print(t(sen_cov[na_rows, psm_covars]))
}


# =============================================================================
# NA
# =============================================================================

sen_cov <- read.csv("E:/2026.1.8_biomod2/SEN_covariates/Senegal_Adm4_covariates.csv")

cat(":", nrow(sen_cov), "\n")

# NA
psm_covars <- c("temp_mean", "prec_mean", "prec_cv", "aridity_mean",
                "elev_mean", "slope_mean", "roughness_mean",
                "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")

na_rows <- which(rowSums(is.na(sen_cov[, psm_covars])) > 0)
cat(":\n")
print(sen_cov[na_rows, c("NAME_4", "NAME_2", "area_km2")])


sen_cov_clean <- sen_cov[-na_rows, ]
cat("\n:", nrow(sen_cov_clean), "\n")

# Treatment
cat("\nTreatment:\n")
cat("  GGW:", sum(sen_cov_clean$treatment == 1), "\n")
cat("  GGW:", sum(sen_cov_clean$treatment == 0), "\n")


write.csv(sen_cov_clean, 
          "E:/2026.1.8_biomod2/SEN_covariates/Senegal_Adm4_covariates_clean.csv",
          row.names = FALSE)
cat("\n✓ : Senegal_Adm4_covariates_clean.csv\n")




# =============================================================================
# VIF
# =============================================================================

library(car)

cat("\n")
cat("╔", rep("═", 65), "╗\n", sep = "")
cat("║   VIF                                               ║\n")
cat("╚", rep("═", 65), "╝\n\n")


sen_cov <- read.csv("E:/2026.1.8_biomod2/SEN_covariates/Senegal_Adm4_covariates_clean.csv")

cat(":", nrow(sen_cov), "\n")
cat("Treatment: GGW=", sum(sen_cov$treatment == 1), 
    ", GGW=", sum(sen_cov$treatment == 0), "\n\n")


psm_covars <- c("temp_mean", "prec_mean", "prec_cv", "aridity_mean",
                "elev_mean", "slope_mean", "roughness_mean",
                "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")

cat(" (", length(psm_covars), "):\n", sep = "")
cat(paste(psm_covars, collapse = ", "), "\n\n")

# =============================================================================
# VIF
# =============================================================================

cat("【VIF】\n")
cat(rep("-", 50), "\n", sep = "")

vif_formula <- as.formula(paste("treatment ~", paste(psm_covars, collapse = " + ")))
vif_model <- glm(vif_formula, data = sen_cov, family = binomial)
vif_values <- car::vif(vif_model)


vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = round(vif_values, 2)
)
vif_df <- vif_df[order(-vif_df$VIF), ]

cat("\nVIF ():\n")
print(vif_df, row.names = FALSE)

# =============================================================================
# VIF
# =============================================================================

cat("\n【VIF】\n")
cat(rep("-", 50), "\n", sep = "")

final_covars <- psm_covars
removed_vars <- c()

while(max(vif_values) > 10) {
  remove_var <- names(which.max(vif_values))
  removed_vars <- c(removed_vars, remove_var)
  cat("  :", remove_var, "(VIF =", round(max(vif_values), 2), ")\n")
  
  final_covars <- setdiff(final_covars, remove_var)
  
  vif_formula <- as.formula(paste("treatment ~", paste(final_covars, collapse = " + ")))
  vif_model <- glm(vif_formula, data = sen_cov, family = binomial)
  vif_values <- car::vif(vif_model)
}

# =============================================================================

# =============================================================================

cat("\n【VIF】\n")
cat(rep("-", 50), "\n", sep = "")

final_vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = round(vif_values, 2)
)
final_vif_df <- final_vif_df[order(-final_vif_df$VIF), ]

print(final_vif_df, row.names = FALSE)

cat("\n (", length(removed_vars), "): ", 
    ifelse(length(removed_vars) > 0, paste(removed_vars, collapse = ", "), ""), "\n", sep = "")

cat("\n (", length(final_covars), "):\n", sep = "")
for(v in final_covars) {
  cat("  -", v, "(VIF =", round(vif_values[v], 2), ")\n")
}

# =============================================================================

# =============================================================================

cat("\n【】\n")
cat(rep("-", 50), "\n", sep = "")

eth_final <- c("prec_mean", "prec_cv", "temp_mean", "slope_mean",
               "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")

cat(" (8):\n")
cat("  ", paste(eth_final, collapse = ", "), "\n")

cat("\n (", length(final_covars), "):\n", sep = "")
cat("  ", paste(final_covars, collapse = ", "), "\n")

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("VIF！\n")
cat(rep("=", 60), "\n", sep = "")

cat("\nPSM:\n")
cat("final_covars <- c(\"", paste(final_covars, collapse = "\", \""), "\")\n", sep = "")


# prec_meanVIF
test_covars <- c("temp_mean", "prec_mean", "elev_mean", "slope_mean",
                 "ndvi_mean", "log_pop", "lulc_cropland", "lulc_diversity")

vif_formula <- as.formula(paste("treatment ~", paste(test_covars, collapse = " + ")))
vif_model <- glm(vif_formula, data = sen_cov, family = binomial)
vif_values <- car::vif(vif_model)

cat("prec_meanaridity_meanVIF:\n")
print(sort(round(vif_values, 2), decreasing = TRUE))