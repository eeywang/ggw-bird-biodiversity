
eth_cov <- read.csv("E:/2026.1.8_biomod2/ETH_covariates/Ethiopia_Adm3_covariates.csv")


cat("===  ===\n")
print(names(eth_cov))


cat("\n===  ===\n")
str(eth_cov)


cat("\n=== 6 ===\n")
print(head(eth_cov))


cat("\n:", nrow(eth_cov), "\n")


2026.1.11

# =============================================================================
# Shannoncropland
# =============================================================================

library(terra)
library(sf)
library(dplyr)

cat("\n")
cat("╔", rep("═", 65), "╗\n", sep = "")
cat("║  ShannonCropland                        ║\n")
cat("╚", rep("═", 65), "╝\n\n")

# =============================================================================

# =============================================================================

cat("【 1】\n")
cat(rep("-", 50), "\n", sep = "")

# LULC
lulc <- rast("E:/2026.1.8_biomod2/ETH_covariates/ETH_LULC_2005.tif")


eth_adm3 <- st_read("E:/11.17progress/study_area/ETH_adm3.shp", quiet = TRUE)
eth_adm3 <- st_transform(eth_adm3, crs(lulc))


eth_cov <- read.csv("E:/2026.1.8_biomod2/ETH_covariates/Ethiopia_Adm3_covariates.csv")

cat("  LULC:", basename(sources(lulc)), "\n")
cat("  :", nrow(eth_adm3), "\n")
cat("  :", nrow(eth_cov), "\n")

# LULC
cat("\n  LULC:\n")
lulc_freq <- freq(lulc)
print(lulc_freq)

# =============================================================================
# LULC
# =============================================================================

cat("\n【 2】LULC\n")
cat(rep("-", 50), "\n", sep = "")

# MODIS MCD12Q1 LC_Type1 :
# 1 = Evergreen Needleleaf Forests
# 2 = Evergreen Broadleaf Forests
# 3 = Deciduous Needleleaf Forests
# 4 = Deciduous Broadleaf Forests
# 5 = Mixed Forests
# 6 = Closed Shrublands
# 7 = Open Shrublands
# 8 = Woody Savannas
# 9 = Savannas
# 10 = Grasslands
# 11 = Permanent Wetlands
# 12 = Croplands
# 13 = Urban and Built-up Lands
# 14 = Cropland/Natural Vegetation Mosaics
# 15 = Permanent Snow and Ice
# 16 = Barren
# 17 = Water Bodies

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

cat("  ", length(lulc_classes), "\n")

# =============================================================================
# Shannon
# =============================================================================

cat("\n【 3】Shannon\n")
cat(rep("-", 50), "\n", sep = "")

# SpatVector
eth_vect <- vect(eth_adm3)


n_adm <- nrow(eth_adm3)
lulc_props <- data.frame(
  GID_3 = eth_adm3$GID_3,
  lulc_forest = numeric(n_adm),
  lulc_shrubland = numeric(n_adm),
  lulc_savanna = numeric(n_adm),
  lulc_grassland = numeric(n_adm),
  lulc_wetland = numeric(n_adm),
  lulc_cropland = numeric(n_adm),
  lulc_urban = numeric(n_adm),
  lulc_barren = numeric(n_adm),
  lulc_water = numeric(n_adm),
  lulc_diversity = numeric(n_adm)
)

# Shannon
calc_shannon <- function(props) {
  props <- props[props > 0]
  if(length(props) == 0) return(0)
  -sum(props * log(props))
}


cat("  ...\n")
pb <- txtProgressBar(min = 0, max = n_adm, style = 3)

for(i in 1:n_adm) {
  
  tryCatch({

    adm_i <- eth_vect[i, ]
    lulc_crop <- crop(lulc, adm_i)
    lulc_mask <- mask(lulc_crop, adm_i)
    

    vals <- values(lulc_mask, na.rm = TRUE)
    
    if(length(vals) > 0) {
      total_pixels <- length(vals)
      

      props <- numeric(length(lulc_classes))
      names(props) <- names(lulc_classes)
      
      for(lc_name in names(lulc_classes)) {
        props[lc_name] <- sum(vals %in% lulc_classes[[lc_name]]) / total_pixels
        lulc_props[i, paste0("lulc_", lc_name)] <- props[lc_name]
      }
      
      # Shannon
      lulc_props$lulc_diversity[i] <- calc_shannon(props)
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

cat("\n【 4】\n")
cat(rep("-", 50), "\n", sep = "")

cat("\n  --- lulc_cropland ---\n")
cat("  Min:", round(min(lulc_props$lulc_cropland, na.rm = TRUE), 4), "\n")
cat("  Max:", round(max(lulc_props$lulc_cropland, na.rm = TRUE), 4), "\n")
cat("  Mean:", round(mean(lulc_props$lulc_cropland, na.rm = TRUE), 4), "\n")
cat("  NA:", sum(is.na(lulc_props$lulc_cropland)), "\n")

cat("\n  --- lulc_diversity (Shannon) ---\n")
cat("  Min:", round(min(lulc_props$lulc_diversity, na.rm = TRUE), 4), "\n")
cat("  Max:", round(max(lulc_props$lulc_diversity, na.rm = TRUE), 4), "\n")
cat("  Mean:", round(mean(lulc_props$lulc_diversity, na.rm = TRUE), 4), "\n")
cat("  NA:", sum(is.na(lulc_props$lulc_diversity)), "\n")


par(mfrow = c(1, 2))
hist(lulc_props$lulc_cropland, breaks = 30, main = "Cropland Proportion", 
     xlab = "Proportion", col = "steelblue")
hist(lulc_props$lulc_diversity, breaks = 30, main = "Shannon Diversity Index", 
     xlab = "Shannon Index", col = "darkgreen")
par(mfrow = c(1, 1))

# =============================================================================

# =============================================================================

cat("\n【 5】\n")
cat(rep("-", 50), "\n", sep = "")


lulc_to_merge <- lulc_props[, c("GID_3", "lulc_cropland", "lulc_diversity")]


if("lulc_cropland" %in% names(eth_cov)) {
  cat("  : lulc_cropland\n")
  eth_cov$lulc_cropland <- NULL
}
if("lulc_diversity" %in% names(eth_cov)) {
  eth_cov$lulc_diversity <- NULL
}


eth_cov_updated <- merge(eth_cov, lulc_to_merge, by = "GID_3", all.x = TRUE)

# log_pop
if(!"log_pop" %in% names(eth_cov_updated)) {
  eth_cov_updated$log_pop <- log(eth_cov_updated$pop_mean + 1)
}

cat("  :", nrow(eth_cov_updated), "\n")
cat("  :", ncol(eth_cov_updated), "\n")


output_file <- "E:/2026.1.8_biomod2/ETH_covariates/Ethiopia_Adm3_covariates_updated.csv"
write.csv(eth_cov_updated, output_file, row.names = FALSE)
cat("\n  ✓ :", output_file, "\n")

# =============================================================================
# VIF
# =============================================================================

cat("\n【 6】VIF\n")
cat(rep("-", 50), "\n", sep = "")

library(car)


psm_covars <- c("prec_mean", "prec_cv", "temp_mean", "aridity_mean",
                "elev_mean", "slope_mean", "roughness_mean",
                "ndvi_mean", "log_pop", 
                "lulc_cropland", "lulc_diversity")


available_covars <- psm_covars[psm_covars %in% names(eth_cov_updated)]
missing_covars <- psm_covars[!psm_covars %in% names(eth_cov_updated)]

if(length(missing_covars) > 0) {
  cat("  : :", paste(missing_covars, collapse = ", "), "\n")
}

cat("  :", paste(available_covars, collapse = ", "), "\n\n")

# VIF
vif_data <- eth_cov_updated[, c("treatment", available_covars)]
vif_data <- na.omit(vif_data)

cat("  VIF:", nrow(vif_data), "\n")

# VIF
vif_formula <- as.formula(paste("treatment ~", paste(available_covars, collapse = " + ")))
vif_model <- glm(vif_formula, data = vif_data, family = binomial)

# VIF
vif_values <- car::vif(vif_model)

cat("\n  VIF:\n")
vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = round(vif_values, 2)
)
vif_df <- vif_df[order(-vif_df$VIF), ]
print(vif_df, row.names = FALSE)

# VIF > 10
if(max(vif_values) > 10) {
  cat("\n  ⚠ VIF > 10:\n")
  high_vif <- names(vif_values[vif_values > 10])
  for(v in high_vif) {
    cat("    -", v, "(VIF =", round(vif_values[v], 2), ")\n")
  }
} else {
  cat("\n  ✓ VIF < 10\n")
}

# VIF
final_covars <- available_covars
while(max(vif_values) > 10) {
  remove_var <- names(which.max(vif_values))
  cat("\n  :", remove_var, "(VIF =", round(max(vif_values), 2), ")\n")
  final_covars <- setdiff(final_covars, remove_var)
  
  vif_formula <- as.formula(paste("treatment ~", paste(final_covars, collapse = " + ")))
  vif_model <- glm(vif_formula, data = vif_data, family = binomial)
  vif_values <- car::vif(vif_model)
}

cat("\n   (VIF < 10):\n")
for(v in final_covars) {
  cat("    -", v, "(VIF =", round(vif_values[v], 2), ")\n")
}

# =============================================================================

# =============================================================================

cat("\n")
cat(rep("=", 65), "\n", sep = "")
cat("！\n")
cat(rep("=", 65), "\n", sep = "")

cat("\n:", output_file, "\n")
cat("\nPSM:\n")
cat(paste(final_covars, collapse = ", "), "\n")