# ==============================================================================
# Guild spatial maps — matching scenario decomposition style
# colour scheme: RdBu (), : guild richness change
# Output: three separate plots (14 × 8 inches), no panel labels (added in Illustrator)
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(sf, ggplot2, dplyr, cowplot, scales,
               rnaturalearth, rnaturalearthdata)

# ==============================================================================
# 1. FILE PATHS
# ==============================================================================

ggw_all_path <- "D:/OneDrive - University of Leeds/Biodiversity/GGWAREA.shp"
nga_shp_path <- "E:/11.17progress/study_area/Nigeria_adm2.shp"
sen_shp_path <- "E:/11.17progress/study_area/SEN_adm4.shp"
eth_shp_path <- "E:/11.17progress/study_area/ETH_adm3.shp"

output_dir <- "D:/OneDrive - University of Leeds/biodiversity_manscript/figures/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

cat("\n========== Loading Data ==========\n")

ggw_all <- st_read(ggw_all_path, quiet = TRUE)
ggw_outline <- ggw_all %>% st_union() %>% st_sf()

nga_admin <- st_read(nga_shp_path, quiet = TRUE)
sen_admin <- st_read(sen_shp_path, quiet = TRUE)
eth_admin <- st_read(eth_shp_path, quiet = TRUE)

nga_outline <- nga_admin %>% st_union() %>% st_sf() %>% mutate(country = "Nigeria")
sen_outline <- sen_admin %>% st_union() %>% st_sf() %>% mutate(country = "Senegal")
eth_outline <- eth_admin %>% st_union() %>% st_sf() %>% mutate(country = "Ethiopia")

# Guild richness CSVs
ng_csv <- read.csv("D:/OneDrive - University of Leeds/biodiversity_manscript/results/Nigeria/guild_richness/Nigeria_guild_richness_by_LGA.csv")
sn_csv <- read.csv("D:/OneDrive - University of Leeds/biodiversity_manscript/results/Senegal/guild_richness/Senegal_guild_richness_by_admin.csv")
et_csv <- read.csv("D:/OneDrive - University of Leeds/biodiversity_manscript/results/Ethiopia/guild_richness/Ethiopia_guild_richness_by_admin.csv")

# Join guild data to admin shapefiles
nga_admin <- nga_admin %>% left_join(ng_csv, by = "GID_2")
sen_admin <- sen_admin %>% left_join(sn_csv, by = "GID_4")
eth_admin <- eth_admin %>% left_join(et_csv, by = "GID_3")

africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")

cat("Data loaded successfully!\n")

# ==============================================================================
# 3. COLOR PALETTE — RdBu (11 colors, )
# ==============================================================================

rdbu_colors <- c(
  "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
  "#f7f7f7",
  "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
)

scale_fill_rdbu <- function(limits, ...) {
  scale_fill_gradientn(
    colors = rdbu_colors,
    limits = limits,
    oob = scales::squish,
    na.value = "gray90",
    ...
  )
}

# ==============================================================================
# 4. DETERMINE GLOBAL SCALE RANGE
# ==============================================================================

all_guild_vals <- c(
  nga_admin$Woodland_change, nga_admin$Openhabitat_change, nga_admin$Wetland_change,
  sen_admin$Woodland_change, sen_admin$Openhabitat_change, sen_admin$Wetland_change,
  eth_admin$Woodland_change, eth_admin$Openhabitat_change, eth_admin$Wetland_change
)
guild_max <- ceiling(max(abs(all_guild_vals), na.rm = TRUE))
guild_max <- max(guild_max, 10)  # At least ±10

cat("Scale for guild maps: ±", guild_max, "\n")

# ==============================================================================
# 5. BOUNDING BOX — Plan B # ==============================================================================

all_countries <- bind_rows(nga_outline, sen_outline, eth_outline)
bbox_all <- st_bbox(all_countries)

x_range <- bbox_all["xmax"] - bbox_all["xmin"]
y_range <- bbox_all["ymax"] - bbox_all["ymin"]
buffer_x <- x_range * 0.08
buffer_y <- y_range * 0.08

xlim_tight <- c(bbox_all["xmin"] - buffer_x, bbox_all["xmax"] + buffer_x)
ylim_tight <- c(bbox_all["ymin"] - buffer_y, bbox_all["ymax"] + buffer_y)

cat("Map extent: X [", xlim_tight[1], ",", xlim_tight[2],
    "] Y [", ylim_tight[1], ",", ylim_tight[2], "]\n")

# ==============================================================================
# 6. THEME — Plan B # ==============================================================================

map_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    legend.key.width = unit(2.5, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid = element_line(color = "gray85", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# ==============================================================================
# 7. MAP FUNCTION — Plan B make_effect_map # ==============================================================================

make_guild_map <- function(nga_data, sen_data, eth_data,
                           var_name, title, subtitle, vmax) {

  nga_mean <- mean(nga_data[[var_name]], na.rm = TRUE)
  sen_mean <- mean(sen_data[[var_name]], na.rm = TRUE)
  eth_mean <- mean(eth_data[[var_name]], na.rm = TRUE)

  nga_centroid <- st_centroid(st_union(nga_data))
  sen_centroid <- st_centroid(st_union(sen_data))
  eth_centroid <- st_centroid(st_union(eth_data))

  p <- ggplot() +
    # Africa background
    geom_sf(data = africa, fill = "#f5f5f5", color = "#cccccc", linewidth = 0.2) +

    # Three countries with guild values
    geom_sf(data = sen_data, aes(fill = .data[[var_name]]), color = "gray50", linewidth = 0.1) +
    geom_sf(data = nga_data, aes(fill = .data[[var_name]]), color = "gray50", linewidth = 0.1) +
    geom_sf(data = eth_data, aes(fill = .data[[var_name]]), color = "gray50", linewidth = 0.1) +

    # Country borders (thick)
    geom_sf(data = sen_outline, fill = NA, color = "#333333", linewidth = 0.8) +
    geom_sf(data = nga_outline, fill = NA, color = "#333333", linewidth = 0.8) +
    geom_sf(data = eth_outline, fill = NA, color = "#333333", linewidth = 0.8) +

    # GGW boundary (green dashed)
    geom_sf(data = ggw_outline, fill = NA, color = "#2e7d32",
            linewidth = 1.2, linetype = "dashed") +

    # Color scale — RdBu
    scale_fill_rdbu(
      limits = c(-vmax, vmax),
      name = expression(Delta ~ "species richness"),
      breaks = seq(-vmax, vmax, length.out = 9),
      labels = function(x) ifelse(x > 0, paste0("+", round(x, 1)), round(x, 1))
    ) +

    # Map extent
    coord_sf(xlim = xlim_tight, ylim = ylim_tight, expand = FALSE) +

    # Labels
    labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude",
      y = "Latitude"
    ) +

    # Country name annotations
    annotate("label",
             x = st_coordinates(sen_centroid)[1],
             y = st_coordinates(sen_centroid)[2],
             label = sprintf("Senegal\n(%+.2f)", sen_mean),
             size = 3.2, fontface = "bold",
             fill = "white", alpha = 0.9, label.size = 0.3) +
    annotate("label",
             x = st_coordinates(nga_centroid)[1],
             y = st_coordinates(nga_centroid)[2],
             label = sprintf("Nigeria\n(%+.2f)", nga_mean),
             size = 3.2, fontface = "bold",
             fill = "white", alpha = 0.9, label.size = 0.3) +
    annotate("label",
             x = st_coordinates(eth_centroid)[1],
             y = st_coordinates(eth_centroid)[2],
             label = sprintf("Ethiopia\n(%+.2f)", eth_mean),
             size = 3.2, fontface = "bold",
             fill = "white", alpha = 0.9, label.size = 0.3) +

    # Theme
    map_theme +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

  return(p)
}

# ==============================================================================
# 8. GENERATE THREE GUILD MAPS
# ==============================================================================

cat("\n========== Generating Guild Maps ==========\n")

# --- Woodland ---
fig_woodland <- make_guild_map(
  nga_admin, sen_admin, eth_admin,
  "Woodland_change",
  "Woodland",
  "Change in woodland bird species richness between early and late periods",
  guild_max
)
ggsave(file.path(output_dir, "Fig4_guild_Woodland.png"),
       fig_woodland, width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig4_guild_Woodland.pdf"),
       fig_woodland, width = 14, height = 8, device = cairo_pdf)
cat("Woodland map saved!\n")

# --- Open-habitat ---
fig_openhabitat <- make_guild_map(
  nga_admin, sen_admin, eth_admin,
  "Openhabitat_change",
  "Open-habitat",
  "Change in open-habitat bird species richness between early and late periods",
  guild_max
)
ggsave(file.path(output_dir, "Fig4_guild_Openhabitat.png"),
       fig_openhabitat, width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig4_guild_Openhabitat.pdf"),
       fig_openhabitat, width = 14, height = 8, device = cairo_pdf)
cat("Open-habitat map saved!\n")

# --- Wetland ---
fig_wetland <- make_guild_map(
  nga_admin, sen_admin, eth_admin,
  "Wetland_change",
  "Wetland",
  "Change in wetland bird species richness between early and late periods",
  guild_max
)
ggsave(file.path(output_dir, "Fig4_guild_Wetland.png"),
       fig_wetland, width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig4_guild_Wetland.pdf"),
       fig_wetland, width = 14, height = 8, device = cairo_pdf)
cat("Wetland map saved!\n")

cat("\n\u2713 All three guild maps saved to:", output_dir, "\n")
