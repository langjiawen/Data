#地图


#采用
install.packages("sf")
install.packages("sp")
install.packages("ggOceanMaps", 
                 repos = c("https://mikkovihtakari.r-universe.dev", 
                           "https://cloud.r-project.org"))
install.packages("ggOceanMapsData", 
                 repos = c("https://mikkovihtakari.r-universe.dev", 
                           "https://cloud.r-project.org"))


library(ggOceanMaps)
library(ggspatial)
library(ggplot2)

options(ggOceanMaps.userdir = "~/ggOceanMapsData_cache")

stations <- data.frame(
  Site = c("A","B","C","D"),
  Longitude = c(-3.747567, -5.389552, -4.920561, -5.600900),
  Latitude = c(54.323133, 53.886122, 55.435981, 53.752522)
)
stations$Depth <- c(68,108,72,80)

p <- basemap(
  limits = c(-8, -2, 53, 56.5),
  bathymetry = TRUE,
  bathy.style = "rcb"
) +
  ggspatial::geom_spatial_point(
    data = stations,
    aes(x = Longitude, y = Latitude),
    shape = 17, colour = "red3", size = 4,
    crs = 4326
  ) +
  ggspatial::geom_spatial_label(
    data = stations,
    aes(x = Longitude, y = Latitude, label = Site),
    nudge_y = 0.08, fontface = "bold", size = 5,
    label.size = 0, fill = NA,
    crs = 4326
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering(fill = c("black", "white"), line_col = "black")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.25,
    pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

print(p)



#最新采用 大加小
install.packages("patchwork")
library(ggOceanMaps)
library(ggspatial)
library(ggplot2)
library(patchwork)   # 用于拼图，没装先 install.packages("patchwork")

options(ggOceanMaps.userdir = "~/ggOceanMapsData_cache")

stations <- data.frame(
  Site = c("A","B","C","D"),
  Longitude = c(-3.747567, -5.389552, -4.920561, -5.600900),
  Latitude = c(54.323133, 53.886122, 55.435981, 53.752522)
)
stations$Depth <- c(68,108,72,80)

##################################################
## 主图：局部站点图（就是你刚画好的那张）
##################################################
main_map <- basemap(
  limits = c(-8, -2, 53, 56.5),
  bathymetry = TRUE,
  bathy.style = "rcb"
) +
  ggspatial::geom_spatial_point(
    data = stations,
    aes(x = Longitude, y = Latitude),
    shape = 17, colour = "red3", size = 4,
    crs = 4326
  ) +
  ggspatial::geom_spatial_label(
    data = stations,
    aes(x = Longitude, y = Latitude, label = Site),
    nudge_y = 0.08, fontface = "bold", size = 5,
    label.size = 0, fill = NA,
    crs = 4326
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering(fill = c("black", "white"), line_col = "black")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.25,
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

##################################################
## 小图：英国+周边海域全图，画一个矩形标出研究区域范围
##################################################
uk_box <- data.frame(
  lon = c(-8, -2, -2, -8, -8),
  lat = c(53, 53, 56.5, 56.5, 53)
)

overview_map <- basemap(
  limits = c(-11, 2, 49, 61),   # 覆盖整个英国+爱尔兰海+北海一部分
  bathymetry = TRUE,
  bathy.style = "rcb"
) +
  ggspatial::geom_spatial_path(
    data = uk_box,
    aes(x = lon, y = lat),
    colour = "red", linewidth = 0.8,
    crs = 4326
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    legend.position = "none"
  )

##################################################
## 拼图：小图嵌进大图右下角
##################################################
final_plot <- main_map +
  inset_element(
    overview_map,
    left = 0.68, bottom = 0.02, right = 0.99, top = 0.31
  )

final_plot

ggsave("stations_map_with_inset.png", plot = final_plot, width = 9, height = 8, dpi = 300)





R.version.string







#site C 5-star
library(ggOceanMaps)
library(ggspatial)
library(ggplot2)

options(ggOceanMaps.userdir = "~/ggOceanMapsData_cache")

# 五个CTD坐标（Site C站位）
ctd_points <- data.frame(
  CTD = c("CTD 1", "CTD 2", "CTD 3", "CTD 4", "CTD 5"),
  Longitude = c(-4.920561, -4.917820, -4.915167, -4.919733, -4.915767),
  Latitude  = c(55.435981, 55.435303, 55.434500, 55.433833, 55.436783)
)

# 星形连线：以CTD1为中心，向其余四点连线（可根据实际布放逻辑调整连线顺序）
star_lines <- data.frame(
  x    = rep(ctd_points$Longitude[1], 4),
  y    = rep(ctd_points$Latitude[1], 4),
  xend = ctd_points$Longitude[-1],
  yend = ctd_points$Latitude[-1]
)

p_ctd <- basemap(
  limits = c(-4.93, -4.905, 55.430, 55.440),
  bathymetry = TRUE,
  bathy.style = "rcb"
) +
  # 星形连线
  geom_segment(
    data = star_lines,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey40", linetype = "dashed", linewidth = 0.4
  ) +
  # CTD站位点
  ggspatial::geom_spatial_point(
    data = ctd_points,
    aes(x = Longitude, y = Latitude),
    shape = 21, fill = "dodgerblue2", colour = "black",
    size = 4, stroke = 0.6,
    crs = 4326
  ) +
  # CTD标签
  ggspatial::geom_spatial_label(
    data = ctd_points,
    aes(x = Longitude, y = Latitude, label = CTD),
    nudge_y = 0.0006, fontface = "bold", size = 4,
    label.size = 0, fill = NA,
    crs = 4326
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
    style = north_arrow_fancy_orienteering(fill = c("black", "white"), line_col = "black")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.3,
    pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in")
  ) +
  ggtitle("Site C: five-star CTD deployment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(p_ctd)

ggsave("SiteC_CTD_deployment.png", plot = p_ctd, width = 8, height = 6, dpi = 300)
