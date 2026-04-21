# ==============================================================================
# Nature/NEE-Compliant Figure Code — COMPLETE VERSION
# 
# Nature Branded Research Journals Artwork Standards:
#   - Font: Helvetica / Arial (sans-serif)
#   - Max text size: 7pt   Min text size: 5pt
#   - 2-column width: 180mm = 7.087 inches
#   - 1-column width: 88mm  = 3.465 inches
#   - Format: Vector (PDF via cairo_pdf)
#   - Color: RGB for original research
#
# Output:
#   *_180mm: full-width version

if (!require("pacman")) install.packages("pacman")
pacman::p_load(sf, ggplot2, dplyr, terra, scales,
               exactextractr, rnaturalearth, rnaturalearthdata,
               cowplot, grid, gridExtra)

# ==============================================================================
# 0. NATURE-COMPLIANT CONSTANTS
# ==============================================================================

# --- ---
WIDTH_2COL   <- 180 / 25.4   # 180mm → 7.087 inches
MAP_H_FULL   <- 60  / 25.4   # 180mm ≈ 60mm

# --- 180mm font sizes
SIZE_AXIS_TEXT    <- 6
SIZE_AXIS_TITLE   <- 7
SIZE_LEGEND_TEXT  <- 5.5
SIZE_LEGEND_TITLE <- 6


# --- 180mm line widths
LW_BORDER     <- 0.5
LW_COUNTRY    <- 0.6
LW_GGW        <- 0.8
LW_ADMIN      <- 0.05
LW_AFRICA     <- 0.15
LW_AXIS_TICK  <- 0.3


# --- 180mm ---
LEGEND_KEY_W  <- 1.5    # cm
LEGEND_KEY_H  <- 0.2    # cm


cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  Nature/NEE Artwork Standards Applied\n")
cat("  Max text: 7pt | Min text: 5pt\n")
cat("  180mm full-width + \n")
cat("══════════════════════════════════════════════════════════════\n\n")


# ==============================================================================
# 1. FILE PATHS
# ==============================================================================

ggw_all_path   <- "D:/OneDrive - University of Leeds/Biodiversity/GGWAREA.shp"
nga_shp_path   <- "E:/11.17progress/study_area/Nigeria_adm2.shp"
sen_shp_path   <- "E:/11.17progress/study_area/SEN_adm4.shp"
eth_shp_path   <- "E:/11.17progress/study_area/ETH_adm3.shp"

nga_effect_dir <- "E:/2026.1.8_biomod2/NGA_scenario_analysis"
sen_effect_dir <- "E:/2026.1.8_biomod2/SEN_scenario_analysis/ggw_senegal/scenario_analysis_50threshold"
eth_effect_dir <- "E:/2026.1.8_biomod2/ETH_scenario_analysis/ggw_ethiopia/scenario_analysis_50threshold"

nga_guild_dir  <- "D:/OneDrive - University of Leeds/biodiversity_manscript/results/Nigeria/guild_richness"
sen_guild_dir  <- "D:/OneDrive - University of Leeds/biodiversity_manscript/results/Senegal/guild_richness"
eth_guild_dir  <- "D:/OneDrive - University of Leeds/biodiversity_manscript/results/Ethiopia/guild_richness"

output_dir     <- "D:/OneDrive - University of Leeds/biodiversity_manscript/figures/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

cat("========== Loading Data ==========\n")

target_crs <- st_crs(4326)

ggw_all     <- st_read(ggw_all_path, quiet = TRUE)
ggw_outline <- ggw_all %>% st_union() %>% st_sf() %>% st_transform(target_crs)

nga_admin <- st_read(nga_shp_path, quiet = TRUE) %>% st_transform(target_crs)
sen_admin <- st_read(sen_shp_path, quiet = TRUE) %>% st_transform(target_crs)
eth_admin <- st_read(eth_shp_path, quiet = TRUE) %>% st_transform(target_crs)

nga_outline <- nga_admin %>% st_union() %>% st_sf() %>% mutate(country = "Nigeria")
sen_outline <- sen_admin %>% st_union() %>% st_sf() %>% mutate(country = "Senegal")
eth_outline <- eth_admin %>% st_union() %>% st_sf() %>% mutate(country = "Ethiopia")

africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf") %>%
  st_transform(target_crs)

# --- Scenario effects ---
load_effects <- function(effect_dir) {
  list(
    total       = rast(file.path(effect_dir, "effect_probability_total.tif")),
    climate     = rast(file.path(effect_dir, "effect_probability_climate.tif")),
    vegetation  = rast(file.path(effect_dir, "effect_probability_vegetation.tif")),
    interaction = rast(file.path(effect_dir, "effect_probability_interaction.tif"))
  )
}

nga_effects <- load_effects(nga_effect_dir)
sen_effects <- load_effects(sen_effect_dir)
eth_effects <- load_effects(eth_effect_dir)

calc_admin_means <- function(admin_shp, effects) {
  admin_shp$effect_total       <- exact_extract(effects$total,       admin_shp, 'mean')
  admin_shp$effect_climate     <- exact_extract(effects$climate,     admin_shp, 'mean')
  admin_shp$effect_vegetation  <- exact_extract(effects$vegetation,  admin_shp, 'mean')
  admin_shp$effect_interaction <- exact_extract(effects$interaction, admin_shp, 'mean')
  return(admin_shp)
}

nga_admin <- calc_admin_means(nga_admin, nga_effects)
sen_admin <- calc_admin_means(sen_admin, sen_effects)
eth_admin <- calc_admin_means(eth_admin, eth_effects)

# --- Guild richness ---
load_guild_rasters <- function(guild_dir) {
  list(
    Woodland    = rast(file.path(guild_dir, "richness_prob_change_Woodland.tif")),
    Openhabitat = rast(file.path(guild_dir, "richness_prob_change_Openhabitat.tif")),
    Wetland     = rast(file.path(guild_dir, "richness_prob_change_Wetland.tif"))
  )
}

nga_guilds <- load_guild_rasters(nga_guild_dir)
sen_guilds <- load_guild_rasters(sen_guild_dir)
eth_guilds <- load_guild_rasters(eth_guild_dir)

calc_guild_means <- function(admin_shp, guild_rasters) {
  for (guild in names(guild_rasters)) {
    r <- guild_rasters[[guild]]
    if (nlyr(r) > 1) r <- r[[1]]
    col_name <- paste0("guild_", guild)
    admin_shp[[col_name]] <- as.numeric(exact_extract(r, admin_shp, 'mean'))
  }
  return(admin_shp)
}

nga_admin <- calc_guild_means(nga_admin, nga_guilds)
sen_admin <- calc_guild_means(sen_admin, sen_guilds)
eth_admin <- calc_guild_means(eth_admin, eth_guilds)

cat("All data loaded successfully!\n")


# ==============================================================================
# 3. SCALE RANGES
# ==============================================================================

all_main_vals <- c(
  nga_admin$effect_total, nga_admin$effect_climate, nga_admin$effect_vegetation,
  sen_admin$effect_total, sen_admin$effect_climate, sen_admin$effect_vegetation,
  eth_admin$effect_total, eth_admin$effect_climate, eth_admin$effect_vegetation
)
main_max <- max(ceiling(max(abs(all_main_vals), na.rm = TRUE)), 25)

all_int_vals <- c(
  nga_admin$effect_interaction,
  sen_admin$effect_interaction,
  eth_admin$effect_interaction
)
int_max <- max(ceiling(max(abs(all_int_vals), na.rm = TRUE)), 3)

wo_vals <- c(
  nga_admin$guild_Woodland, nga_admin$guild_Openhabitat,
  sen_admin$guild_Woodland, sen_admin$guild_Openhabitat,
  eth_admin$guild_Woodland, eth_admin$guild_Openhabitat
)
wo_max <- max(ceiling(max(abs(wo_vals), na.rm = TRUE)), 10)

wet_vals <- c(
  nga_admin$guild_Wetland,
  sen_admin$guild_Wetland,
  eth_admin$guild_Wetland
)
wet_max <- max(ceiling(max(abs(wet_vals), na.rm = TRUE)), 5)

cat("Scenario main: \u00B1", main_max, " | Interaction: \u00B1", int_max, "\n")
cat("Guild WO: \u00B1", wo_max, " | Wetland: \u00B1", wet_max, "\n")


# ==============================================================================
# 4. BOUNDING BOX
# ==============================================================================

all_countries <- bind_rows(nga_outline, sen_outline, eth_outline)
all_geoms     <- bind_rows(all_countries, ggw_outline %>% mutate(country = "GGW"))
bbox_all      <- st_bbox(all_geoms)

x_range <- bbox_all["xmax"] - bbox_all["xmin"]
y_range <- bbox_all["ymax"] - bbox_all["ymin"]

xlim_tight <- c(bbox_all["xmin"] - x_range * 0.03, bbox_all["xmax"] + x_range * 0.03)
ylim_tight <- c(bbox_all["ymin"] - y_range * 0.05, bbox_all["ymax"] + y_range * 0.05)


# ==============================================================================
# 5. COLOR PALETTES
# ==============================================================================

rdylgn_colors <- c(
  "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b",
  "#ffffbf",
  "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"
)

rdbu_colors <- c(
  "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
  "#f7f7f7",
  "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
)


# ==============================================================================
# 6A. THEME — 180mm full-width version
# ==============================================================================

theme_nature <- theme(
  panel.grid.major   = element_blank(),
  panel.grid.minor   = element_blank(),
  panel.background   = element_rect(fill = "white", color = NA),
  panel.border       = element_rect(fill = NA, color = "black", linewidth = LW_BORDER),
  plot.background    = element_rect(fill = "white", color = NA),
  text               = element_text(family = "Helvetica"),
  axis.text          = element_text(size = SIZE_AXIS_TEXT,   color = "black"),
  axis.title         = element_text(size = SIZE_AXIS_TITLE,  color = "black"),
  axis.title.x       = element_text(margin = margin(3, 0, 0, 0)),
  axis.title.y       = element_text(margin = margin(0, 3, 0, 0)),
  axis.ticks         = element_line(color = "black", linewidth = LW_AXIS_TICK),
  axis.ticks.length  = unit(0.08, "cm"),
  legend.position    = "bottom",
  legend.key.width   = unit(LEGEND_KEY_W, "cm"),
  legend.key.height  = unit(LEGEND_KEY_H, "cm"),
  legend.title       = element_text(size = SIZE_LEGEND_TITLE, family = "Helvetica"),
  legend.text        = element_text(size = SIZE_LEGEND_TEXT,  family = "Helvetica"),
  legend.margin      = margin(0, 0, 0, 0),
  legend.box.margin  = margin(-3, 0, 0, 0),
  plot.margin        = margin(2, 3, 2, 2),
  plot.title         = element_blank(),
  plot.subtitle      = element_blank()
)


# ==============================================================================
theme_nature_half <- theme(
  panel.grid.major   = element_blank(),
  panel.grid.minor   = element_blank(),
  panel.background   = element_rect(fill = "white", color = NA),
  panel.border       = element_rect(fill = NA, color = "black", linewidth = LW_BORDER_H),
  plot.background    = element_rect(fill = "white", color = NA),
  text               = element_text(family = "Helvetica"),
  axis.text          = element_text(size = SIZE_AXIS_TEXT_H,   color = "black"),
  axis.title         = element_text(size = SIZE_AXIS_TITLE_H,  color = "black"),
  axis.title.x       = element_text(margin = margin(2, 0, 0, 0)),
  axis.title.y       = element_text(margin = margin(0, 2, 0, 0)),
  axis.ticks         = element_line(color = "black", linewidth = LW_AXIS_TICK_H),
  axis.ticks.length  = unit(0.05, "cm"),
  legend.position    = "bottom",
  legend.key.width   = unit(LEGEND_KEY_W_H, "cm"),
  legend.key.height  = unit(LEGEND_KEY_H_H, "cm"),
  legend.title       = element_text(size = SIZE_LEGEND_TITLE_H, family = "Helvetica"),
  legend.text        = element_text(size = SIZE_LEGEND_TEXT_H,  family = "Helvetica"),
  legend.margin      = margin(0, 0, 0, 0),
  legend.box.margin  = margin(-3, 0, 0, 0),
  plot.margin        = margin(1, 2, 1, 1),
  plot.title         = element_blank(),
  plot.subtitle      = element_blank()
)


# ==============================================================================
# 7A. Scenario map — 180mm # ==============================================================================

make_scenario_map <- function(nga_data, sen_data, eth_data,
                              var_name, vmax, is_interaction = FALSE,
                              show_x_axis = TRUE, show_y_axis = TRUE) {
  
  brks <- if (is_interaction) seq(-vmax, vmax, length.out = 7) else seq(-vmax, vmax, length.out = 9)
  scale_label <- paste0("\u0394 species richness (scale: \u00B1", vmax, ")")
  
  p <- ggplot() +
    geom_sf(data = africa, fill = "#f0f0f0", color = "#cccccc", linewidth = LW_AFRICA) +
    geom_sf(data = sen_data, aes(fill = .data[[var_name]]), color = "gray60", linewidth = LW_ADMIN) +
    geom_sf(data = nga_data, aes(fill = .data[[var_name]]), color = "gray60", linewidth = LW_ADMIN) +
    geom_sf(data = eth_data, aes(fill = .data[[var_name]]), color = "gray60", linewidth = LW_ADMIN) +
    geom_sf(data = sen_outline, fill = NA, color = "black", linewidth = LW_COUNTRY) +
    geom_sf(data = nga_outline, fill = NA, color = "black", linewidth = LW_COUNTRY) +
    geom_sf(data = eth_outline, fill = NA, color = "black", linewidth = LW_COUNTRY) +
    geom_sf(data = ggw_outline, fill = NA, color = "#2e7d32", linewidth = LW_GGW, linetype = "dashed") +
    scale_fill_gradientn(
      colors = rdylgn_colors, limits = c(-vmax, vmax), oob = scales::squish,
      na.value = "gray90", name = scale_label, breaks = brks,
      labels = function(x) ifelse(x > 0, paste0("+", round(x, 1)), round(x, 1))
    ) +
    coord_sf(xlim = xlim_tight, ylim = ylim_tight, expand = FALSE) +
    labs(
      x = if (show_x_axis) expression("Longitude ("*degree*"E)") else NULL,
      y = if (show_y_axis) expression("Latitude ("*degree*"N)")  else NULL
    ) +
    theme_nature +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
  
  if (!show_x_axis) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  if (!show_y_axis) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  return(p)
}


# ==============================================================================
# 7B. Guild map — 180mm # ==============================================================================

make_guild_map <- function(nga_data, sen_data, eth_data,
                           var_name, vmax,
                           show_x_axis = TRUE, show_y_axis = TRUE) {
  
  brks        <- seq(-vmax, vmax, length.out = 9)
  scale_label <- paste0("\u0394 species richness (scale: \u00B1", vmax, ")")
  
  p <- ggplot() +
    geom_sf(data = africa, fill = "#f0f0f0", color = "#cccccc", linewidth = LW_AFRICA) +
    geom_sf(data = sen_data, aes(fill = .data[[var_name]]), color = "gray60", linewidth = LW_ADMIN) +
    geom_sf(data = nga_data, aes(fill = .data[[var_name]]), color = "gray60", linewidth = LW_ADMIN) +
    geom_sf(data = eth_data, aes(fill = .data[[var_name]]), color = "gray60", linewidth = LW_ADMIN) +
    geom_sf(data = sen_outline, fill = NA, color = "black", linewidth = LW_COUNTRY) +
    geom_sf(data = nga_outline, fill = NA, color = "black", linewidth = LW_COUNTRY) +
    geom_sf(data = eth_outline, fill = NA, color = "black", linewidth = LW_COUNTRY) +
    geom_sf(data = ggw_outline, fill = NA, color = "#2e7d32", linewidth = LW_GGW, linetype = "dashed") +
    scale_fill_gradientn(
      colors = rdbu_colors, limits = c(-vmax, vmax), oob = scales::squish,
      na.value = "gray90", name = scale_label, breaks = brks,
      labels = function(x) ifelse(x > 0, paste0("+", round(x, 1)), round(x, 1))
    ) +
    coord_sf(xlim = xlim_tight, ylim = ylim_tight, expand = FALSE) +
    labs(
      x = if (show_x_axis) expression("Longitude ("*degree*"E)") else NULL,
      y = if (show_y_axis) expression("Latitude ("*degree*"N)")  else NULL
    ) +
    theme_nature +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
  
  if (!show_x_axis) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  if (!show_y_axis) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  return(p)
}


# ==============================================================================
# ==============================================================================
# ==============================================================================
# 8. SAVE FUNCTION
# ==============================================================================

save_nature <- function(plot, name, w = WIDTH_2COL, h = MAP_H_FULL) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")),
         plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(name, ".png")),
         plot, width = w, height = h, dpi = 300, bg = "white")
  cat("\u2713", name, sprintf("(%.0f \u00d7 %.0f mm)\n", w * 25.4, h * 25.4))
}


# ==============================================================================
# PART A: Fig 2 — Scenario Decomposition Maps
# ==============================================================================

cat("\n========== PART A: Fig 2 Scenario Maps ==========\n")

# --- 180mm full-width version ---
cat("  [180mm version]\n")

fig2a_180 <- make_scenario_map(nga_admin, sen_admin, eth_admin, "effect_total",       main_max)
fig2b_180 <- make_scenario_map(nga_admin, sen_admin, eth_admin, "effect_climate",     main_max)
fig2c_180 <- make_scenario_map(nga_admin, sen_admin, eth_admin, "effect_vegetation",  main_max)
fig2d_180 <- make_scenario_map(nga_admin, sen_admin, eth_admin, "effect_interaction", int_max, is_interaction = TRUE)

save_nature(fig2a_180, "Fig2a_total_180mm")
save_nature(fig2b_180, "Fig2b_climate_180mm")
save_nature(fig2c_180, "Fig2c_vegetation_180mm")
save_nature(fig2d_180, "Fig2d_interaction_180mm")



# ==============================================================================
# PART B: Fig 4 — Guild Spatial Maps
# ==============================================================================

cat("\n========== PART B: Fig 4 Guild Maps ==========\n")

# --- 180mm full-width version ---
cat("  [180mm version]\n")

fig4a_180 <- make_guild_map(nga_admin, sen_admin, eth_admin, "guild_Woodland",    wo_max)
fig4b_180 <- make_guild_map(nga_admin, sen_admin, eth_admin, "guild_Openhabitat", wo_max)
fig4c_180 <- make_guild_map(nga_admin, sen_admin, eth_admin, "guild_Wetland",     wet_max)

save_nature(fig4a_180, "Fig4a_woodland_180mm")
save_nature(fig4b_180, "Fig4b_openhabitat_180mm")
save_nature(fig4c_180, "Fig4c_wetland_180mm")



# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  All Done!\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  Output:", output_dir, "\n\n")
cat("  180mm versions: *_180mm.pdf / .png\n")
cat("    Font: 5.5-7pt | Lines: 0.05-0.8 | Legend: 1.5cm bar\n")
cat("    Use for: standalone figures, single-column full width\n\n")
cat("  \n")
cat("    Font: 5-5.5pt | Lines: 0.03-0.5 | Legend: 0.8cm bar\n")
cat("    Use for: 2-column composite in Adobe Illustrator\n\n")
cat("  Panel labels (a)(b)(c)(d)(e) → add in AI\n")
cat("  Python boxplots → see separate scripts\n")
cat("══════════════════════════════════════════════════════════════\n")
