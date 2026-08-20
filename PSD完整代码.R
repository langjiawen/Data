#PSD完整代码

#Site B
library(tidyverse)

############################################################
## 1. Read Function (补齐了头部定义)
############################################################
read_uvp <- function(file) {
  lines <- readLines(file)
  dat <- read.delim(
    file,
    skip = 7,
    sep = ";",
    header = FALSE
  )
  colnames(dat) <- strsplit(lines[7], ";")[[1]]
  return(dat)
}

# 读取数据
pre  <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_032_odv.txt")
post <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_037_odv.txt")

############################################################
## 2. Settings
############################################################
depth_col <- "Depth [m]:PRIMARYVAR:DOUBLE"

abundance_cols <- grep("\\[# l-1\\]", names(pre))

mid_size <- c(57.4, 72.3, 91.3, 115, 144.5, 182, 229.5, 289.5, 364.5, 459, 578.5)

selected_pattern <- paste(
  c("50.8-64", "64-80.6", "80.6-102", "102-128", "128-161", 
    "161-203", "203-256", "256-323", "323-406", "406-512", "512-645"),
  collapse = "|"
)

############################################################
## 3. Extract PSD for Selected Depths
############################################################
extract_psd <- function(dat, depth, tow) {
  x <- dat %>% filter(.data[[depth_col]] == depth)
  
  if(nrow(x) == 0) return(NULL) # 增加容错判断
  
  out <- data.frame(
    SizeBin = names(dat)[abundance_cols],
    Abundance = as.numeric(x[1, abundance_cols])
  ) %>%
    filter(grepl(selected_pattern, SizeBin))
  
  out$MidSize <- mid_size
  out$Depth   <- paste0(depth, " m")
  out$Tow     <- tow
  
  return(out)
}

# 组合数据
psd <- bind_rows(
  extract_psd(pre, 47.5, "Pre"),
  extract_psd(pre, 102.5, "Pre"),
  extract_psd(post, 47.5, "Post"),
  extract_psd(post, 102.5, "Post")
) %>%
  mutate(
    Station = "B",
    Layer   = ifelse(Depth == "47.5 m", "Shallow", "Deep"),
    Group   = paste(Depth, Tow),
    Abundance = ifelse(Abundance <= 0, NA, Abundance)
  )

############################################################
## 4. PSD Figure
############################################################
ggplot(psd, aes(x = MidSize, y = Abundance, colour = Depth, linetype = Tow, shape = Tow, group = interaction(Depth, Tow))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c("47.5 m" = "black", "102.5 m" = "grey50")) +
  scale_linetype_manual(values = c("Pre" = "solid", "Post" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Site B",
    x = expression("Particle diameter (" * mu * "m)"),
    y = expression("Particle abundance (# " * L^{-1} * ")")
  )

############################################################
## 5. Calculate Slopes Across All Depths (重构优化版)
############################################################
# 定义一个计算单份数据所有深度斜率的辅助函数
calc_depth_slopes <- function(dat, tow_label) {
  depths_list <- sort(unique(dat[[depth_col]]))
  res <- data.frame()
  
  for (d in depths_list) {
    psd_depth <- dat %>% filter(.data[[depth_col]] == d)
    
    psd_tmp <- data.frame(
      SizeBin   = names(dat)[abundance_cols],
      Abundance = as.numeric(psd_depth[1, abundance_cols])
    ) %>%
      filter(grepl(selected_pattern, SizeBin)) %>%
      mutate(
        MidSize   = mid_size,
        Abundance = ifelse(Abundance <= 0, NA, Abundance)
      ) %>%
      na.omit()
    
    # 至少需要 6 个点才进行回归拟合
    if (nrow(psd_tmp) >= 6) {
      model <- lm(log10(Abundance) ~ log10(MidSize), data = psd_tmp)
      res <- rbind(res, data.frame(
        Station = "B",
        Tow     = tow_label,
        Depth   = d,
        Slope   = coef(model)[2]
      ))
    }
  }
  return(res)
}

# 修正变量名，一次性算完 pre 和 post
slope_df <- rbind(
  calc_depth_slopes(pre, "Pre"),
  calc_depth_slopes(post, "Post")
)

# 打印结果
print(slope_df)

############################################################
## 6. Slope Boxplot
############################################################
slope_df$Tow <- factor(slope_df$Tow, levels = c("Pre", "Post"))
ggplot(
  slope_df,
  aes(
    x = Tow,
    y = Slope,
    fill = Tow
  )
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.6
  ) +
  # =========================================================
# 关键修改：加上 show.legend = FALSE，点就不会混进图例中了
# =========================================================
geom_jitter(
  width = 0.1,
  size = 2,
  show.legend = FALSE
) +
  scale_fill_manual(
    values = c("Pre" = "grey70", "Post" = "grey30")
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Site B",
    x = "",
    y = "PSD slope"
  )

# 默认的 Welch's t-test (不假设两组方差相等，论文中最推荐的做法)
t_res <- t.test(Slope ~ Tow, data = slope_df)
print(t_res)


############################################################
## 7. Statistical Analysis (T-test & Assumption Checks)
############################################################
library(dplyr)

# ----------------------------------------------------
# A. 检查正态性 (Shapiro-Wilk Test)
# ----------------------------------------------------
normality_check <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    N = n(),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE),
    Shapiro_p = shapiro.test(Slope)$p.value
  )

cat("--- 1. 正态性检验结果 (p > 0.05 说明符合正态分布) ---\n")
print(normality_check)

# ----------------------------------------------------
# B. 执行 T 检验 (Welch's T-test)
# ----------------------------------------------------
# 注：Welch's t-test 比传统学生 t 检验更稳健，无需强求两组方差相等
t_results <- t.test(Slope ~ Tow, data = slope_df, var.equal = FALSE)

cat("\n--- 2. T 检验详细结果 ---\n")
print(t_results)

# 提取关键论文指标
t_stat <- round(t_results$statistic, 3)
df_val <- round(t_results$parameter, 2)
p_val  <- t_results$p.value

cat("\n--- 3. 论文格式描述 (APA Style) ---\n")
if (p_val < 0.001) {
  p_text <- "p < 0.001"
} else {
  p_text <- paste0("p = ", round(p_val, 4))
}

cat(paste0("A Welch's two-sample t-test showed ",
           ifelse(p_val < 0.05, "a significant change", "no significant difference"),
           " in PSD slope between Pre-trawling and Post-trawling (t(", 
           df_val, ") = ", t_stat, ", ", p_text, ").\n"))

# ----------------------------------------------------
# C. 若数据不符合正态分布，可补充非参数检验 (Wilcoxon Test)
# ----------------------------------------------------
# wilcox_res <- wilcox.test(Slope ~ Tow, data = slope_df)
# print(wilcox_res)


install.packages("ggpubr") # 如果没安装过需要运行此行
library(ggpubr)

ggplot(
  slope_df,
  aes(
    x = Tow,
    y = Slope,
    fill = Tow
  )
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.6,
    outlier.shape = NA # 隐藏默认离群点，避免与 geom_jitter 重复
  ) +
  geom_jitter(
    width = 0.1,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("Pre" = "grey70", "Post" = "grey30")
  ) +
  # 自动计算 T 检验并在图上方标注显著性标记 (如 ns, *, **, ***)
  stat_compare_means(
    method = "t.test",
    label = "p.format",    # 显示星号/ns，若想直接显示数值可改为 label = "p.format"
    label.x = 1.5,         # 标记横坐标位置 (Pre和Post中间)
    label.y = max(slope_df$Slope, na.rm = TRUE) + 0.1, # 标记高度
    size = 5
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") + # 已经有 x 轴了，可隐藏颜色图例
  labs(
    title = "Site B",
    x = "",
    y = "PSD slope"
  )

print(t_results)

summary_values <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    n = sum(!is.na(Slope)),
    Median = median(Slope, na.rm = TRUE),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE)
  )

print(summary_values)





#site C

library(tidyverse)
library(ggpubr) # 确保已安装 ggpubr 用于在图上标注 T 检验结果

############################################################
## 1. Read Function
############################################################
read_uvp <- function(file) {
  lines <- readLines(file)
  dat <- read.delim(
    file,
    skip = 7,
    sep = ";",
    header = FALSE
  )
  colnames(dat) <- strsplit(lines[7], ";")[[1]]
  return(dat)
}

# 读取 Site C 的 Pre/Post 数据
pre  <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_002_odv.txt")
post <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_007_odv.txt")

############################################################
## 2. Settings
############################################################
depth_col <- "Depth [m]:PRIMARYVAR:DOUBLE"

abundance_cols <- grep("\\[# l-1\\]", names(pre))

mid_size <- c(57.4, 72.3, 91.3, 115, 144.5, 182, 229.5, 289.5, 364.5, 459, 578.5)

selected_pattern <- paste(
  c("50.8-64", "64-80.6", "80.6-102", "102-128", "128-161", 
    "161-203", "203-256", "256-323", "323-406", "406-512", "512-645"),
  collapse = "|"
)

############################################################
## 3. Extract PSD for Selected Depths (Site C: 32.5 m & 67.5 m)
############################################################
extract_psd <- function(dat, depth, tow) {
  x <- dat %>% filter(.data[[depth_col]] == depth)
  
  if(nrow(x) == 0) return(NULL) # 容错判断
  
  out <- data.frame(
    SizeBin = names(dat)[abundance_cols],
    Abundance = as.numeric(x[1, abundance_cols])
  ) %>%
    filter(grepl(selected_pattern, SizeBin))
  
  out$MidSize <- mid_size
  out$Depth   <- paste0(depth, " m")
  out$Tow     <- tow
  
  return(out)
}

# 组合 Site C 数据
psd <- bind_rows(
  extract_psd(pre, 32.5, "Pre"),
  extract_psd(pre, 67.5, "Pre"),
  extract_psd(post, 32.5, "Post"),
  extract_psd(post, 67.5, "Post")
) %>%
  mutate(
    Station   = "C",
    Layer     = ifelse(Depth == "32.5 m", "Shallow", "Deep"),
    Group     = paste(Depth, Tow),
    Abundance = ifelse(Abundance <= 0, NA, Abundance)
  )

############################################################
## 4. PSD Figure (Site C)
############################################################
ggplot(psd, aes(x = MidSize, y = Abundance, colour = Depth, linetype = Tow, shape = Tow, group = interaction(Depth, Tow))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c("32.5 m" = "black", "67.5 m" = "grey50")) +
  scale_linetype_manual(values = c("Pre" = "solid", "Post" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Site C",
    x = expression("Particle diameter (" * mu * "m)"),
    y = expression("Particle abundance (# " * L^{-1} * ")")
  )

############################################################
## 5. Calculate Slopes Across All Depths (Site C)
############################################################
calc_depth_slopes <- function(dat, tow_label) {
  depths_list <- sort(unique(dat[[depth_col]]))
  res <- data.frame()
  
  for (d in depths_list) {
    psd_depth <- dat %>% filter(.data[[depth_col]] == d)
    
    psd_tmp <- data.frame(
      SizeBin   = names(dat)[abundance_cols],
      Abundance = as.numeric(psd_depth[1, abundance_cols])
    ) %>%
      filter(grepl(selected_pattern, SizeBin)) %>%
      mutate(
        MidSize   = mid_size,
        Abundance = ifelse(Abundance <= 0, NA, Abundance)
      ) %>%
      na.omit()
    
    if (nrow(psd_tmp) >= 6) {
      model <- lm(log10(Abundance) ~ log10(MidSize), data = psd_tmp)
      res <- rbind(res, data.frame(
        Station = "C",
        Tow     = tow_label,
        Depth   = d,
        Slope   = coef(model)[2]
      ))
    }
  }
  return(res)
}

slope_df <- rbind(
  calc_depth_slopes(pre, "Pre"),
  calc_depth_slopes(post, "Post")
)

slope_df$Tow <- factor(slope_df$Tow, levels = c("Pre", "Post"))
print(slope_df)

############################################################
## 6. Statistical Analysis for Site C (T-test & Assumptions)
############################################################
# A. 正态性检验
normality_check <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    N = n(),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE),
    Shapiro_p = shapiro.test(Slope)$p.value
  )

cat("\n====================================================\n")
cat("--- 1. Site C 正态性检验 (p > 0.05 符合正态分布) ---\n")
print(normality_check)

# B. 执行 Welch's T 检验
t_results <- t.test(Slope ~ Tow, data = slope_df, var.equal = FALSE)

cat("\n--- 2. Site C T 检验详细结果 ---\n")
print(t_results)

# C. 自动生成论文描述 (APA 格式)
t_stat <- round(t_results$statistic, 3)
df_val <- round(t_results$parameter, 2)
p_val  <- t_results$p.value
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 4)))

cat("\n--- 3. 论文格式描述 (APA Style) ---\n")
cat(paste0("For Site C, a Welch's two-sample t-test revealed ",
           ifelse(p_val < 0.05, "a significant difference", "no significant difference"),
           " in PSD slope between Pre- and Post-trawling (t(", 
           df_val, ") = ", t_stat, ", ", p_text, ").\n"))
cat("====================================================\n\n")

############################################################
## 7. Slope Boxplot with T-test Significance (Site C)
############################################################
ggplot(
  slope_df,
  aes(
    x = Tow,
    y = Slope,
    fill = Tow
  )
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("Pre" = "grey70", "Post" = "grey30")
  ) +
  # 自动在箱形图顶部标注 T 检验 p 值
  stat_compare_means(
    method = "t.test",
    label = "p.format",                                  # 若想显示星号/ns可改为 "p.signif"
    label.x = 1.5,                                       # 标注放置在横轴中间
    label.y = max(slope_df$Slope, na.rm = TRUE) + 0.1,   # 标记放置在顶部高度
    size = 5
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Site C",
    x = "",
    y = "PSD slope"
  )

print(t_results)

summary_values <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    n = sum(!is.na(Slope)),
    Median = median(Slope, na.rm = TRUE),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE)
  )

print(summary_values)






#site D
library(tidyverse)
library(ggpubr) # 确保已安装 ggpubr 用于在图上标注 T 检验结果

############################################################
## 1. Read Function
############################################################
read_uvp <- function(file) {
  lines <- readLines(file)
  dat <- read.delim(
    file,
    skip = 7,
    sep = ";",
    header = FALSE
  )
  colnames(dat) <- strsplit(lines[7], ";")[[1]]
  return(dat)
}

# 读取 Site D 的 Pre/Post 数据
pre  <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_087_odv.txt")
post <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_092_odv.txt")

############################################################
## 2. Settings
############################################################
depth_col <- "Depth [m]:PRIMARYVAR:DOUBLE"

abundance_cols <- grep("\\[# l-1\\]", names(pre))

mid_size <- c(57.4, 72.3, 91.3, 115, 144.5, 182, 229.5, 289.5, 364.5, 459, 578.5)

selected_pattern <- paste(
  c("50.8-64", "64-80.6", "80.6-102", "102-128", "128-161", 
    "161-203", "203-256", "256-323", "323-406", "406-512", "512-645"),
  collapse = "|"
)

############################################################
## 3. Extract PSD for Selected Depths (Site D: 32.5 m & 67.5 m)
############################################################
extract_psd <- function(dat, depth, tow) {
  x <- dat %>% filter(.data[[depth_col]] == depth)
  
  if(nrow(x) == 0) return(NULL) # 容错判断
  
  out <- data.frame(
    SizeBin = names(dat)[abundance_cols],
    Abundance = as.numeric(x[1, abundance_cols])
  ) %>%
    filter(grepl(selected_pattern, SizeBin))
  
  out$MidSize <- mid_size
  out$Depth   <- paste0(depth, " m")
  out$Tow     <- tow
  
  return(out)
}

# 组合 Site D 数据
psd <- bind_rows(
  extract_psd(pre, 32.5, "Pre"),
  extract_psd(pre, 67.5, "Pre"),
  extract_psd(post, 32.5, "Post"),
  extract_psd(post, 67.5, "Post")
) %>%
  mutate(
    Station   = "D",
    Layer     = ifelse(Depth == "32.5 m", "Shallow", "Deep"),
    Group     = paste(Depth, Tow),
    Abundance = ifelse(Abundance <= 0, NA, Abundance)
  )

############################################################
## 4. PSD Figure (Site D)
############################################################
ggplot(psd, aes(x = MidSize, y = Abundance, colour = Depth, linetype = Tow, shape = Tow, group = interaction(Depth, Tow))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c("32.5 m" = "black", "67.5 m" = "grey50")) +
  scale_linetype_manual(values = c("Pre" = "solid", "Post" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Site D",
    x = expression("Particle diameter (" * mu * "m)"),
    y = expression("Particle abundance (# " * L^{-1} * ")")
  )

############################################################
## 5. Calculate Slopes Across All Depths (Site D)
############################################################
calc_depth_slopes <- function(dat, tow_label) {
  depths_list <- sort(unique(dat[[depth_col]]))
  res <- data.frame()
  
  for (d in depths_list) {
    psd_depth <- dat %>% filter(.data[[depth_col]] == d)
    
    psd_tmp <- data.frame(
      SizeBin   = names(dat)[abundance_cols],
      Abundance = as.numeric(psd_depth[1, abundance_cols])
    ) %>%
      filter(grepl(selected_pattern, SizeBin)) %>%
      mutate(
        MidSize   = mid_size,
        Abundance = ifelse(Abundance <= 0, NA, Abundance)
      ) %>%
      na.omit()
    
    if (nrow(psd_tmp) >= 6) {
      model <- lm(log10(Abundance) ~ log10(MidSize), data = psd_tmp)
      res <- rbind(res, data.frame(
        Station = "D",
        Tow     = tow_label,
        Depth   = d,
        Slope   = coef(model)[2]
      ))
    }
  }
  return(res)
}

slope_df <- rbind(
  calc_depth_slopes(pre, "Pre"),
  calc_depth_slopes(post, "Post")
)

slope_df$Tow <- factor(slope_df$Tow, levels = c("Pre", "Post"))
print(slope_df)

############################################################
## 6. Statistical Analysis for Site D (T-test & Assumptions)
############################################################
# A. 正态性检验
normality_check <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    N = n(),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE),
    Shapiro_p = shapiro.test(Slope)$p.value
  )

cat("\n====================================================\n")
cat("--- 1. Site D 正态性检验 (p > 0.05 符合正态分布) ---\n")
print(normality_check)

# B. 执行 Welch's T 检验
t_results <- t.test(Slope ~ Tow, data = slope_df, var.equal = FALSE)

cat("\n--- 2. Site D T 检验详细结果 ---\n")
print(t_results)

# C. 自动生成论文描述 (APA 格式)
t_stat <- round(t_results$statistic, 3)
df_val <- round(t_results$parameter, 2)
p_val  <- t_results$p.value
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 4)))

cat("\n--- 3. 论文格式描述 (APA Style) ---\n")
cat(paste0("For Site D, a Welch's two-sample t-test revealed ",
           ifelse(p_val < 0.05, "a significant difference", "no significant difference"),
           " in PSD slope between Pre- and Post-trawling (t(", 
           df_val, ") = ", t_stat, ", ", p_text, ").\n"))
cat("====================================================\n\n")

############################################################
## 7. Slope Boxplot with T-test Significance (Site D)
############################################################
ggplot(
  slope_df,
  aes(
    x = Tow,
    y = Slope,
    fill = Tow
  )
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("Pre" = "grey70", "Post" = "grey30")
  ) +
  # 自动在箱形图顶部标注 T 检验 p 值
  stat_compare_means(
    method = "t.test",
    label = "p.format",                                  # 若想显示星号/ns可改为 "p.signif"
    label.x = 1.5,                                       # 标注放置在横轴中间
    label.y = max(slope_df$Slope, na.rm = TRUE) + 0.1,   # 标记放置在顶部高度
    size = 5
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Site D",
    x = "",
    y = "PSD slope"
  )


print(t_results)

summary_values <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    n = sum(!is.na(Slope)),
    Median = median(Slope, na.rm = TRUE),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE)
  )

print(summary_values)








#site A
library(tidyverse)
library(ggpubr) # 确保已安装 ggpubr 用于在图上标注 T 检验结果

############################################################
## 1. Read Function
############################################################
read_uvp <- function(file) {
  lines <- readLines(file)
  dat <- read.delim(
    file,
    skip = 7,
    sep = ";",
    header = FALSE
  )
  colnames(dat) <- strsplit(lines[7], ";")[[1]]
  return(dat)
}

# 读取 Site A 的 Pre/Post 数据
pre  <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_069_odv.txt")
post <- read_uvp("/Users/langjiawen/Desktop/uvp/dy206_ctd_074_odv.txt")

############################################################
## 2. Settings
############################################################
depth_col <- "Depth [m]:PRIMARYVAR:DOUBLE"

abundance_cols <- grep("\\[# l-1\\]", names(pre))

mid_size <- c(57.4, 72.3, 91.3, 115, 144.5, 182, 229.5, 289.5, 364.5, 459, 578.5)

selected_pattern <- paste(
  c("50.8-64", "64-80.6", "80.6-102", "102-128", "128-161", 
    "161-203", "203-256", "256-323", "323-406", "406-512", "512-645"),
  collapse = "|"
)

############################################################
## 3. Extract PSD for Selected Depths (Site A: 17.5 m & 32.5 m)
############################################################
extract_psd <- function(dat, depth, tow) {
  x <- dat %>% filter(.data[[depth_col]] == depth)
  
  if(nrow(x) == 0) return(NULL) # 容错判断
  
  out <- data.frame(
    SizeBin = names(dat)[abundance_cols],
    Abundance = as.numeric(x[1, abundance_cols])
  ) %>%
    filter(grepl(selected_pattern, SizeBin))
  
  out$MidSize <- mid_size
  out$Depth   <- paste0(depth, " m")
  out$Tow     <- tow
  
  return(out)
}

# 组合 Site A 数据
psd <- bind_rows(
  extract_psd(pre, 17.5, "Pre"),
  extract_psd(pre, 32.5, "Pre"),
  extract_psd(post, 17.5, "Post"),
  extract_psd(post, 32.5, "Post")
) %>%
  mutate(
    Station   = "A",
    Layer     = ifelse(Depth == "17.5 m", "Shallow", "Deep"),
    Group     = paste(Depth, Tow),
    Abundance = ifelse(Abundance <= 0, NA, Abundance)
  )

############################################################
## 4. PSD Figure (Site A)
############################################################
ggplot(psd, aes(x = MidSize, y = Abundance, colour = Depth, linetype = Tow, shape = Tow, group = interaction(Depth, Tow))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c("17.5 m" = "black", "32.5 m" = "grey50")) +
  scale_linetype_manual(values = c("Pre" = "solid", "Post" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Site A",
    x = expression("Particle diameter (" * mu * "m)"),
    y = expression("Particle abundance (# " * L^{-1} * ")")
  )

############################################################
## 5. Calculate Slopes Across All Depths (Site A)
############################################################
calc_depth_slopes <- function(dat, tow_label) {
  depths_list <- sort(unique(dat[[depth_col]]))
  res <- data.frame()
  
  for (d in depths_list) {
    psd_depth <- dat %>% filter(.data[[depth_col]] == d)
    
    psd_tmp <- data.frame(
      SizeBin   = names(dat)[abundance_cols],
      Abundance = as.numeric(psd_depth[1, abundance_cols])
    ) %>%
      filter(grepl(selected_pattern, SizeBin)) %>%
      mutate(
        MidSize   = mid_size,
        Abundance = ifelse(Abundance <= 0, NA, Abundance)
      ) %>%
      na.omit()
    
    if (nrow(psd_tmp) >= 6) {
      model <- lm(log10(Abundance) ~ log10(MidSize), data = psd_tmp)
      res <- rbind(res, data.frame(
        Station = "A",
        Tow     = tow_label,
        Depth   = d,
        Slope   = coef(model)[2]
      ))
    }
  }
  return(res)
}

slope_df <- rbind(
  calc_depth_slopes(pre, "Pre"),
  calc_depth_slopes(post, "Post")
)

slope_df$Tow <- factor(slope_df$Tow, levels = c("Pre", "Post"))
print(slope_df)

############################################################
## 6. Statistical Analysis for Site A (T-test & Assumptions)
############################################################
# A. 正态性检验
normality_check <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    N = n(),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE),
    Shapiro_p = shapiro.test(Slope)$p.value
  )

cat("\n====================================================\n")
cat("--- 1. Site A 正态性检验 (p > 0.05 符合正态分布) ---\n")
print(normality_check)

# B. 执行 Welch's T 检验
t_results <- t.test(Slope ~ Tow, data = slope_df, var.equal = FALSE)

cat("\n--- 2. Site A T 检验详细结果 ---\n")
print(t_results)

# C. 自动生成论文描述 (APA 格式)
t_stat <- round(t_results$statistic, 3)
df_val <- round(t_results$parameter, 2)
p_val  <- t_results$p.value
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 4)))

cat("\n--- 3. 论文格式描述 (APA Style) ---\n")
cat(paste0("For Site A, a Welch's two-sample t-test revealed ",
           ifelse(p_val < 0.05, "a significant difference", "no significant difference"),
           " in PSD slope between Pre- and Post-trawling (t(", 
           df_val, ") = ", t_stat, ", ", p_text, ").\n"))
cat("====================================================\n\n")

############################################################
## 7. Slope Boxplot with T-test Significance (Site A)
############################################################
ggplot(
  slope_df,
  aes(
    x = Tow,
    y = Slope,
    fill = Tow
  )
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("Pre" = "grey70", "Post" = "grey30")
  ) +
  # 自动在箱形图顶部标注 T 检验 p 值
  stat_compare_means(
    method = "t.test",
    label = "p.format",                                  # 若想显示星号/ns可改为 "p.signif"
    label.x = 1.5,                                       # 标注放置在横轴中间
    label.y = max(slope_df$Slope, na.rm = TRUE) + 0.1,   # 标记放置在顶部高度
    size = 5
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Site A",
    x = "",
    y = "PSD slope"
  )


print(t_results)

summary_values <- slope_df %>%
  group_by(Tow) %>%
  summarise(
    n = sum(!is.na(Slope)),
    Median = median(Slope, na.rm = TRUE),
    Mean = mean(Slope, na.rm = TRUE),
    SD = sd(Slope, na.rm = TRUE)
  )

print(summary_values)
