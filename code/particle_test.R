#ecotaxa 新 —— Total Particle Abundance (整根水柱) —— 逐站位 Pre vs Post 检验

# ==============================================================================
# 1. 环境准备与数据导入
# ==============================================================================
library(data.table)
library(dplyr)

# 设置工作目录
setwd("/Users/langjiawen/Desktop/Ecotaxa-export-data")

# 读取 TSV 导出文件
data <- fread("TSV_21783_20260713_1410.tsv")

# ==============================================================================
# 2. 筛选有效 CTD 站位（排除错误的 dy206_ctd_063，即 Site A - Pre）
# ==============================================================================
selected_ctd <- c(
  "dy206_ctd_002",
  "dy206_ctd_007",
  "dy206_ctd_032",
  "dy206_ctd_038",
  # "dy206_ctd_063", # Site A - Pre，原始数据有误，已排除，下面手动补充总数
  "dy206_ctd_074",
  "dy206_ctd_087",
  "dy206_ctd_092"
)

# 提取筛选站位的数据（整根水柱，不做深度切片）
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
# 3. 计算每个 Site x Tow 整根水柱的总颗粒数（不限类别，不分层）
# ==============================================================================
total_summary <- particle %>%
  filter(!is.na(Site) & !is.na(Tow)) %>%
  group_by(Site, Tow) %>%
  summarise(
    Total_Count = n(),
    .groups = "drop"
  )

# ==============================================================================
# 4. 手动追加缺失的 Site A - Pre 整根水柱总颗粒数
# ==============================================================================
# TODO: 把下面的 NA 换成 dy206_ctd_069 整根水柱（所有深度、所有类别）的真实总颗粒数
siteA_pre_total <- data.frame(
  Site = "A",
  Tow = "Pre",
  Total_Count = NA_integer_   # <- 在这里填入正确的总数
)

total_abundance_all <- bind_rows(total_summary, siteA_pre_total)

total_abundance_all$Site <- factor(total_abundance_all$Site, levels = c("A", "B", "C", "D"))
total_abundance_all$Tow  <- factor(total_abundance_all$Tow, levels = c("Pre", "Post"))

print(total_abundance_all)
stopifnot(nrow(total_abundance_all) == 8)

if (any(is.na(total_abundance_all$Total_Count))) {
  stop("Site A - Pre 的 Total_Count 还是 NA，请先在第 4 步填入真实总数再做检验。")
}

# ==============================================================================
# 5. 逐站位检验：每个 Site 自身 Pre vs Post 是否显著
# ==============================================================================
# 说明：
# 每个 Site 每个 Tow 阶段只有一个总数，没有重复观测，无法做经典 t 检验
# （t 检验需要组内多个观测值来估计方差）。
# 这里把 Total_Count 当作计数数据（Poisson 过程），用 poisson.test 逐站位
# 比较 Pre 和 Post 两个计数是否有显著差异（默认假设两个阶段的采样量/曝光相同；
# 如果 Pre、Post 实际采样体积或时长不同，需要把对应的 exposure 值填入下面的
# exposure_pre / exposure_post 向量）。

site_test_results <- lapply(levels(total_abundance_all$Site), function(s) {
  
  sub <- total_abundance_all %>% filter(Site == s)
  
  pre_count  <- sub %>% filter(Tow == "Pre")  %>% pull(Total_Count)
  post_count <- sub %>% filter(Tow == "Post") %>% pull(Total_Count)
  
  # 如果 Pre / Post 的采样体积或时长不同，把 1 换成实际的 exposure（如采样体积、时长）
  exposure_pre  <- 1
  exposure_post <- 1
  
  test_res <- poisson.test(
    x = c(pre_count, post_count),
    T = c(exposure_pre, exposure_post)
  )
  
  data.frame(
    Site = s,
    Pre_Count = pre_count,
    Post_Count = post_count,
    Rate_Ratio_Post_over_Pre = as.numeric(test_res$estimate),
    p_value = test_res$p.value
  )
})

site_test_results <- bind_rows(site_test_results)
print(site_test_results)
