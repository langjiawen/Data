#ecotaxa 新

# ==============================================================================
# 1. 环境准备与数据导入
# ==============================================================================
library(data.table)
library(dplyr)
library(ggplot2)

# 设置工作目录
setwd("/Users/langjiawen/Desktop/Ecotaxa-export-data")

# 读取 TSV 导出文件
data <- fread("TSV_21783_20260713_1410.tsv")

# ==============================================================================
# 2. 筛选有效 CTD 站位（排除错误的 dy206_ctd_063）
# ==============================================================================
selected_ctd <- c(
  "dy206_ctd_002",
  "dy206_ctd_007",
  "dy206_ctd_032",
  "dy206_ctd_038",
  # "dy206_ctd_063", # 已排除错误的 063 数据
  "dy206_ctd_074",
  "dy206_ctd_087",
  "dy206_ctd_092"
)

# 提取筛选站位的数据
particle <- subset(data, sample_stationid %in% selected_ctd)

# 分配 Site 与 Tow 属性
particle$Site <- NA
particle$Tow  <- NA

particle$Site[particle$sample_stationid %in% c("dy206_ctd_002","dy206_ctd_007")] <- "C"
particle$Site[particle$sample_stationid %in% c("dy206_ctd_032","dy206_ctd_038")] <- "B"
particle$Site[particle$sample_stationid %in% c("dy206_ctd_074")]               <- "A"
particle$Site[particle$sample_stationid %in% c("dy206_ctd_087","dy206_ctd_092")] <- "D"

particle$Tow[particle$sample_stationid %in% c("dy206_ctd_002", "dy206_ctd_032", "dy206_ctd_087")] <- "Pre"
particle$Tow[particle$sample_stationid %in% c("dy206_ctd_007", "dy206_ctd_038", "dy206_ctd_074", "dy206_ctd_092")] <- "Post"

# ==============================================================================
# 3. 定义元数据表 (Metadata)
# ==============================================================================
metadata <- data.frame(
  Site = c(
    "A","A",
    "B","B","B","B",
    "C","C","C","C",
    "D","D","D","D"
  ),
  Tow = c(
    "Post","Post",
    "Pre","Post","Pre","Post",
    "Pre","Post","Pre","Post",
    "Pre","Post","Pre","Post"
  ),
  CTD = c(
    "dy206_ctd_074","dy206_ctd_074",
    "dy206_ctd_032","dy206_ctd_038","dy206_ctd_032","dy206_ctd_038",
    "dy206_ctd_002","dy206_ctd_007","dy206_ctd_002","dy206_ctd_007",
    "dy206_ctd_087","dy206_ctd_092","dy206_ctd_087","dy206_ctd_092"
  ),
  Layer = c(
    "Low","High",
    "Low","Low","High","High",
    "Low","Low","High","High",
    "Low","Low","High","High"
  ),
  Depth = c(
    12.5, 32.5,
    32.5, 32.5, 62.5, 62.5,
    32.5, 32.5, 67.5, 67.5,
    32.5, 32.5, 62.5, 62.5
  )
)

# ==============================================================================
# 4. 根据深度切片筛选颗粒数据
# ==============================================================================
particle_selected <- data.frame()

for(i in 1:nrow(metadata)){
  temp <- subset(
    particle,
    sample_stationid == metadata$CTD[i] &
      object_depth_min >= metadata$Depth[i] - 2.5 &
      object_depth_min <  metadata$Depth[i] + 2.5
  )
  
  temp$Site  <- metadata$Site[i]
  temp$Tow   <- metadata$Tow[i]
  temp$Layer <- metadata$Layer[i]
  
  particle_selected <- rbind(particle_selected, temp)
}

# ==============================================================================
# 5. 计算颗粒物数量与采样体积
# ==============================================================================
particle_summary <- particle_selected %>%
  filter(object_annotation_category %in% c(
    "aggregates",
    "detritus<detritus",
    "faecal pellet"
  )) %>%
  group_by(Site, Tow, Layer, object_annotation_category) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  )

volume_summary <- particle_selected %>%
  group_by(Site, Tow, Layer) %>%
  summarise(
    n_images = n(),
    total_volume_L = n() * 0.57,
    .groups = "drop"
  )

particle_summary <- left_join(
  particle_summary,
  volume_summary,
  by = c("Site", "Tow", "Layer")
)

particle_summary$Abundance <- particle_summary$Count / particle_summary$total_volume_L

# ==============================================================================
# 6. 手动追加正确的 063 数据 (Site A, Tow Pre)
# ==============================================================================
siteA_pre_fixed <- data.frame(
  Site = "A",
  Tow = "Pre",
  Layer = c("Low", "Low", "Low", "High", "High", "High"),
  object_annotation_category = c(
    "aggregates", "detritus<detritus", "faecal pellet",
    "aggregates", "detritus<detritus", "faecal pellet"
  ),
  Count = c(
    122, 17, 258,   # Low 水层: Aggregates, Detritus, Faecal pellet
    376, 56, 1069   # High 水层: Aggregates, Detritus, Faecal pellet
  ),
  n_images = c(
    1310, 1310, 1310, # Low 水层图片张数
    2289, 2289, 2289  # High 水层图片张数
  )
) %>%
  mutate(
    total_volume_L = n_images * 0.57,
    Abundance = Count / total_volume_L
  )

particle_summary_all <- bind_rows(particle_summary, siteA_pre_fixed)

# ==============================================================================
# 7. 计算比例 (Proportion) 并调整因子顺序（关键修改）
# ==============================================================================
particle_prop <- particle_summary_all %>%
  group_by(Site, Tow, Layer) %>%
  mutate(
    Proportion = Count / sum(Count)
  ) %>%
  ungroup() %>%
  mutate(
    Category = case_when(
      object_annotation_category == "aggregates" ~ "Aggregates",
      object_annotation_category == "detritus<detritus" ~ "Detritus",
      object_annotation_category == "faecal pellet" ~ "Faecal pellet",
      TRUE ~ object_annotation_category
    )
  )

# --- 因子顺序设置（决定画图时的先后顺序）---
particle_prop$Category <- factor(
  particle_prop$Category,
  levels = c("Aggregates", "Detritus", "Faecal pellet")
)

particle_prop$Site <- factor(
  particle_prop$Site,
  levels = c("A", "B", "C", "D")
)

# 1. 强制 Tow 的顺序为 Pre 在前，Post 在后
particle_prop$Tow <- factor(
  particle_prop$Tow,
  levels = c("Pre", "Post")
)

# 2. 强制 Layer 的分版顺序为 Low 在前，High 在后
particle_prop$Layer <- factor(
  particle_prop$Layer,
  levels = c("Low", "High")
)

# 3. 显式生成按 A.Pre -> A.Post -> B.Pre -> B.Post ... 排序的横轴因子
site_tow_levels <- as.vector(t(outer(c("A", "B", "C", "D"), c("Pre", "Post"), paste, sep = ".")))

particle_prop$Site_Tow <- factor(
  paste(particle_prop$Site, particle_prop$Tow, sep = "."),
  levels = site_tow_levels
)

# ==============================================================================
# 8. 绘制堆叠柱状图并显式打印
# ==============================================================================
p <- ggplot(particle_prop,
            aes(x = Site_Tow,
                y = Proportion,
                fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~Layer) +
  labs(
    x = "Site & Tow Stage",
    y = "Proportion",
    fill = "Particle type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 强制输出图像到 Plots 窗口
print(p)



# ==============================================================================
# 1. 加载所需软件包
# ==============================================================================
library(dplyr)
library(tidyr)
library(vegan)     # 用于 NMDS 和 PERMANOVA (adonis2)
library(ggplot2)

# ==============================================================================
# 2. 数据整理：转换为“样方-物种”矩阵 (Site/Sample × Particle Category)
# ==============================================================================
# 将数据整理为每个样本一行的频数表（或丰度表）
nmds_wide <- particle_summary_all %>%
  select(Site, Tow, Layer, object_annotation_category, Count) %>%
  pivot_wider(
    names_from = object_annotation_category,
    values_from = Count,
    values_fill = 0
  ) %>%
  mutate(Sample_ID = paste(Site, Tow, Layer, sep = "_"))

# 提取环境/分组元数据 (Metadata)
metadata_nmds <- nmds_wide %>% select(Sample_ID, Site, Tow, Layer)

# 提取矩阵数据并转换为纯数值矩阵 (Community Matrix)
count_matrix <- nmds_wide %>%
  select(-Sample_ID, -Site, -Tow, -Layer) %>%
  as.matrix()

# ==============================================================================
# 3. 运行 NMDS 分析 (使用 Bray-Curtis 距离)
# ==============================================================================
set.seed(123) # 设置随机种子以保证结果可重复

# k=2 表示降到 2 维空间；trymax 为最大迭代次数
nmds_res <- metaMDS(count_matrix, distance = "bray", k = 2, trymax = 100)

# 查看 Stress 值（压力值）：<0.2 通常表示拟合良好，<0.1 表示极好
print(paste("NMDS Stress value:", round(nmds_res$stress, 4)))

# ==============================================================================
# 4. 运行 PERMANOVA 统计检验 (adonis2)
# ==============================================================================
# 检验拖网处理 (Tow) 和水层 (Layer) 对颗粒物整体构成的显著性影响
permanova_res <- adonis2(count_matrix ~ Tow * Layer, data = metadata_nmds, method = "bray", permutations = 999)
print("--- PERMANOVA 检验结果 ---")
print(permanova_res)

# ==============================================================================
# 5. 提取 NMDS 坐标并准备画图
# ==============================================================================
# 提取每个样本在 NMDS1 和 NMDS2 轴上的坐标
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$Sample_ID <- metadata_nmds$Sample_ID
nmds_scores$Site <- metadata_nmds$Site
nmds_scores$Tow <- metadata_nmds$Tow
nmds_scores$Layer <- metadata_nmds$Layer

# 因子化，保持图例顺序一致
nmds_scores$Tow <- factor(nmds_scores$Tow, levels = c("Pre", "Post"))
nmds_scores$Layer <- factor(nmds_scores$Layer, levels = c("Low", "High"))

# ==============================================================================
# 6. 使用 ggplot2 绘制高质量 NMDS 散点图
# ==============================================================================
stress_text <- paste("Stress =", round(nmds_res$stress, 3))

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Tow, shape = Layer)) +
  geom_point(size = 4, alpha = 0.8) +
  # 按照 Tow 分组绘制 95% 置信椭圆 (若每个组样本数较少，ellipse 可能会有警告，属正常现象)
  stat_ellipse(aes(group = Tow, fill = Tow), geom = "polygon", alpha = 0.15, level = 0.95, color = NA) +
  labs(
    title = "NMDS Quantification of Particle Composition Shift",
    subtitle = stress_text,
    x = "NMDS 1",
    y = "NMDS 2",
    color = "Tow Stage",
    fill = "Tow Stage",
    shape = "Turbidity Layer"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey92", linetype = "dashed"),
    legend.position = "right"
  )

# 强制在 Plots 窗口渲染图像
print(p_nmds)
