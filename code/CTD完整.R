#ctd pre post对比图

########################################################
## CTD Profile Analysis
## Imperial College London
########################################################

library(oce)
library(dplyr)

setwd("/Users/langjiawen/Desktop/DY206_CTD008_align_CTM_Derive_2Hz_Strip")

###############################
## Function to process CTDs
###############################

process_ctd_group <- function(files){
  all_data <- list()
  
  for(i in seq_along(files)){
    ctd <- read.ctd(files[i])
    
    df <- data.frame(
      depth = ctd[["depth"]],
      turb  = ctd[["turbidity"]],
      temp  = ctd[["temperature"]],
      sal   = ctd[["salinity"]],
      chl   = ctd[["fluorescence"]]
    )
    
    df <- df %>%
      filter(complete.cases(.)) %>%
      filter(depth > 2)
    
    df$depth_bin <- floor(df$depth/2)*2
    
    df_bin <- df %>%
      group_by(depth_bin) %>%
      summarise(
        turb = mean(turb),
        chl  = mean(chl),
        temp = mean(temp),
        sal  = mean(sal),
        .groups = "drop"
      )
    
    all_data[[i]] <- df_bin
  }
  
  combined <- bind_rows(all_data)
  
  summary_profile <- combined %>%
    group_by(depth_bin) %>%
    summarise(
      turb_mean = mean(turb),
      turb_sd   = sd(turb),
      chl_mean  = mean(chl),
      chl_sd    = sd(chl),
      temp_mean = mean(temp),
      temp_sd   = sd(temp),
      sal_mean  = mean(sal),
      sal_sd    = sd(sal),
      .groups   = "drop"
    )
  
  return(summary_profile)
}

###############################
## Load Data
###############################

files_pre <- c(
  "DY206_CTD087_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD088_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD089_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD090_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD091_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_pre <- process_ctd_group(files_pre)

files_post <- c(
  "DY206_CTD092_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD093_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD094_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD095_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD096_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_post <- process_ctd_group(files_post)

########################################################
## Plot CTD profiles
########################################################

plot_site <- function(mean_profile_pre,
                      mean_profile_post,
                      low_depth,
                      high_depth,
                      site_name){
  
  # Axis ranges
  turb_range <- range(
    c(
      mean_profile_pre$turb_mean - mean_profile_pre$turb_sd,
      mean_profile_pre$turb_mean + mean_profile_pre$turb_sd,
      mean_profile_post$turb_mean - mean_profile_post$turb_sd,
      mean_profile_post$turb_mean + mean_profile_post$turb_sd
    ),
    na.rm = TRUE
  )
  
  chl_range <- range(
    c(
      mean_profile_pre$chl_mean - mean_profile_pre$chl_sd,
      mean_profile_pre$chl_mean + mean_profile_pre$chl_sd,
      mean_profile_post$chl_mean - mean_profile_post$chl_sd,
      mean_profile_post$chl_mean + mean_profile_post$chl_sd
    ),
    na.rm = TRUE
  )
  
  depth_range <- range(
    c(mean_profile_pre$depth_bin, mean_profile_post$depth_bin)
  )
  
  # 预留右侧充足边距（13行字符宽度）
  par(mar = c(5, 5, 6, 13))
  
  # ----------------------------------------------------
  # 1. Turbidity Axis & Layers
  # ----------------------------------------------------
  plot(
    mean_profile_pre$turb_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    ylim = rev(depth_range),
    xlim = turb_range,
    xlab = expression("Mean turbidity ("*m^{-1}~sr^{-1}*")"),
    ylab = "Depth (m)"
  )
  title(main = site_name, line = 4)
  
  # Turbidity SD Shading
  polygon(
    c(mean_profile_pre$turb_mean - mean_profile_pre$turb_sd,
      rev(mean_profile_pre$turb_mean + mean_profile_pre$turb_sd)),
    c(mean_profile_pre$depth_bin, rev(mean_profile_pre$depth_bin)),
    col = rgb(0, 0, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(mean_profile_post$turb_mean - mean_profile_post$turb_sd,
      rev(mean_profile_post$turb_mean + mean_profile_post$turb_sd)),
    c(mean_profile_post$depth_bin, rev(mean_profile_post$depth_bin)),
    col = rgb(0, 0, 0, 0.08),
    border = NA
  )
  
  # Mean Curves
  lines(mean_profile_pre$turb_mean, mean_profile_pre$depth_bin, lwd = 2, lty = 1, col = "black")
  lines(mean_profile_post$turb_mean, mean_profile_post$depth_bin, lwd = 2, lty = 2, col = "black")
  
  # Selected Depths
  low_pre   <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = low_depth)
  high_pre  <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = high_depth)
  low_post  <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = low_depth)
  high_post <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = high_depth)
  
  points(low_pre$y, low_pre$x, pch = 16, col = "blue", cex = 1.4)
  points(high_pre$y, high_pre$x, pch = 16, col = "red", cex = 1.4)
  points(low_post$y, low_post$x, pch = 17, col = "blue", cex = 1.4)
  points(high_post$y, high_post$x, pch = 17, col = "red", cex = 1.4)
  
  text(low_pre$y, low_pre$x, "Low", pos = 4, col = "blue", cex = 0.85)
  text(high_pre$y, high_pre$x, "High", pos = 4, col = "red", cex = 0.85)
  
  # ----------------------------------------------------
  # 2. Chlorophyll Axis & Layers
  # ----------------------------------------------------
  par(new = TRUE)
  
  plot(
    mean_profile_pre$chl_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = rev(depth_range),
    xlim = chl_range
  )
  
  # Chl SD Shading
  polygon(
    c(mean_profile_pre$chl_mean - mean_profile_pre$chl_sd,
      rev(mean_profile_pre$chl_mean + mean_profile_pre$chl_sd)),
    c(mean_profile_pre$depth_bin, rev(mean_profile_pre$depth_bin)),
    col = rgb(0, 0.39, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(mean_profile_post$chl_mean - mean_profile_post$chl_sd,
      rev(mean_profile_post$chl_mean + mean_profile_post$chl_sd)),
    c(mean_profile_post$depth_bin, rev(mean_profile_post$depth_bin)),
    col = rgb(0.13, 0.55, 0.13, 0.08),
    border = NA
  )
  
  # Chl Mean Curves
  lines(mean_profile_pre$chl_mean, mean_profile_pre$depth_bin, col = "darkgreen", lwd = 2, lty = 1)
  lines(mean_profile_post$chl_mean, mean_profile_post$depth_bin, col = "darkgreen", lwd = 2, lty = 2)
  
  axis(3)
  mtext("Mean chlorophyll (μg/L)", side = 3, line = 2)
  
  # ----------------------------------------------------
  # 3. Legend (改用精确绝对坐标与 XPD 控制)
  # ----------------------------------------------------
  par(xpd = TRUE) # 允许绘制到外边距区域
  
  # 计算图例的绝对 X, Y 坐标：
  # x_pos: 位于当前坐标系 (Chl) 右界外侧一点
  # y_pos: 位于纵轴顶部 (最小深度)
  x_pos <- chl_range[2] + diff(chl_range) * 0.08
  y_pos <- depth_range[1]
  
  legend(
    x = x_pos,
    y = y_pos,
    legend = c(
      "Turbidity Pre",
      "Turbidity Post",
      "Turbidity SD",
      "Chl Pre",
      "Chl Post",
      "Chl SD",
      "Low turbidity",
      "High turbidity"
    ),
    col = c("black", "black", NA, "darkgreen", "darkgreen", NA, "blue", "red"),
    fill = c(NA, NA, rgb(0, 0, 0, 0.15), NA, NA, rgb(0, 0.39, 0, 0.15), NA, NA),
    border = c(NA, NA, NA, NA, NA, NA, NA, NA),
    lty = c(1, 2, NA, 1, 2, NA, NA, NA),
    pch = c(NA, NA, NA, NA, NA, NA, 16, 16),
    lwd = c(2, 2, NA, 2, 2, NA, NA, NA),
    pt.cex = 1.3,
    bty = "n",
    xjust = 0 # 左对齐，确保文字向右自然延伸且留有足够 Margins
  )
}

# ----------------------------------------------------
# 执行保存（导出为 PNG 图片）
# ----------------------------------------------------

png(
  filename = "Site_D_CTD_Profile.png", 
  width = 10,       # 10x8 比例配合 mar=13 能完美呈现图例
  height = 8, 
  units = "in", 
  res = 300
)

plot_site(
  mean_profile_pre,
  mean_profile_post,
  low_depth = 32.0,
  high_depth = 70.5,
  site_name = "Site D (Sand)"
)

dev.off()

library(dplyr)

target_low  <- 32.5
target_high <- 70.5

# 合并数据
merged_df <- inner_join(
  mean_profile_pre, mean_profile_post,
  by = "depth_bin", suffix = c("_pre", "_post")
)

# 动态找到数据中离 17.5 和 30.5 最近的实际 depth_bin 数值
actual_low_depth  <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_low))]
actual_high_depth <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_high))]

target_p_values <- merged_df %>%
  filter(depth_bin %in% c(actual_low_depth, actual_high_depth)) %>%
  mutate(
    Layer = ifelse(depth_bin == actual_low_depth, "Low turbidity layer", "High turbidity layer"),
    n_pre = 5, n_post = 5,
    
    # Turbidity T-test
    se_turb = sqrt((turb_sd_pre^2 / n_pre) + (turb_sd_post^2 / n_post)),
    t_turb  = (turb_mean_post - turb_mean_pre) / se_turb,
    p_turb  = 2 * pt(abs(t_turb), df = n_pre + n_post - 2, lower.tail = FALSE),
    
    # Chlorophyll T-test
    se_chl  = sqrt((chl_sd_pre^2 / n_pre) + (chl_sd_post^2 / n_post)),
    t_chl   = (chl_mean_post - chl_mean_pre) / se_chl,
    p_chl   = 2 * pt(abs(t_chl), df = n_pre + n_post - 2, lower.tail = FALSE)
  ) %>%
  select(Layer, depth_bin, t_turb, p_turb, t_chl, p_chl) %>%
  mutate(
    turb_sig = case_when(p_turb < 0.001 ~ "***", p_turb < 0.01 ~ "**", p_turb < 0.05 ~ "*", TRUE ~ "ns"),
    chl_sig  = case_when(p_chl < 0.001 ~ "***", p_chl < 0.01 ~ "**", p_chl < 0.05 ~ "*", TRUE ~ "ns")
  )

# 打印结果
print("=== Low & High 深度 Turbidity 与 Chlorophyll T检验结果 ===")
print(target_p_values)

range(mean_profile_pre$turb_mean, na.rm = TRUE)
range(mean_profile_post$turb_mean, na.rm = TRUE)

range(mean_profile_pre$chl_mean, na.rm = TRUE)
range(mean_profile_post$chl_mean, na.rm = TRUE)







#siteC

########################################################
## Plot CTD profiles
########################################################

files_pre <- c(
  "DY206_CTD002_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD003_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD004_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD005_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD006_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_pre <- process_ctd_group(files_pre)

files_post <- c(
  "DY206_CTD007_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD008_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD009_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD010_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD011_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_post <- process_ctd_group(files_post)


plot_site <- function(mean_profile_pre,
                      mean_profile_post,
                      low_depth,
                      high_depth,
                      site_name){
  
  # ----------------------------------------------------
  # 0. 数据清洗：剔除 NA 并按深度排序（彻底解决 SD 阴影乱线问题）
  # ----------------------------------------------------
  clean_data <- function(df, mean_col, sd_col) {
    df %>%
      filter(!is.na(.data[[mean_col]]), !is.na(.data[[sd_col]])) %>%
      arrange(depth_bin)
  }
  
  pre_turb_clean  <- clean_data(mean_profile_pre, "turb_mean", "turb_sd")
  post_turb_clean <- clean_data(mean_profile_post, "turb_mean", "turb_sd")
  pre_chl_clean   <- clean_data(mean_profile_pre, "chl_mean", "chl_sd")
  post_chl_clean  <- clean_data(mean_profile_post, "chl_mean", "chl_sd")
  
  # Axis ranges
  turb_range <- range(
    c(
      pre_turb_clean$turb_mean - pre_turb_clean$turb_sd,
      pre_turb_clean$turb_mean + pre_turb_clean$turb_sd,
      post_turb_clean$turb_mean - post_turb_clean$turb_sd,
      post_turb_clean$turb_mean + post_turb_clean$turb_sd
    ),
    na.rm = TRUE
  )
  
  chl_range <- range(
    c(
      pre_chl_clean$chl_mean - pre_chl_clean$chl_sd,
      pre_chl_clean$chl_mean + pre_chl_clean$chl_sd,
      post_chl_clean$chl_mean - post_chl_clean$chl_sd,
      post_chl_clean$chl_mean + post_chl_clean$chl_sd
    ),
    na.rm = TRUE
  )
  
  depth_range <- range(
    c(mean_profile_pre$depth_bin, mean_profile_post$depth_bin),
    na.rm = TRUE
  )
  
  # 预留右侧充足边距存放图例
  par(mar = c(5, 5, 6, 13))
  
  # ----------------------------------------------------
  # 1. Turbidity Axis & Layers
  # ----------------------------------------------------
  plot(
    mean_profile_pre$turb_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    ylim = rev(depth_range),
    xlim = turb_range,
    xlab = expression("Mean turbidity ("*m^{-1}~sr^{-1}*")"),
    ylab = "Depth (m)"
  )
  title(main = site_name, line = 4)
  
  # Turbidity SD Shading
  polygon(
    c(pre_turb_clean$turb_mean - pre_turb_clean$turb_sd,
      rev(pre_turb_clean$turb_mean + pre_turb_clean$turb_sd)),
    c(pre_turb_clean$depth_bin, rev(pre_turb_clean$depth_bin)),
    col = rgb(0, 0, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(post_turb_clean$turb_mean - post_turb_clean$turb_sd,
      rev(post_turb_clean$turb_mean + post_turb_clean$turb_sd)),
    c(post_turb_clean$depth_bin, rev(post_turb_clean$depth_bin)),
    col = rgb(0, 0, 0, 0.08),
    border = NA
  )
  
  # Mean Curves
  lines(mean_profile_pre$turb_mean, mean_profile_pre$depth_bin, lwd = 2, lty = 1, col = "black")
  lines(mean_profile_post$turb_mean, mean_profile_post$depth_bin, lwd = 2, lty = 2, col = "black")
  
  # Selected Depths
  low_pre   <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = low_depth)
  high_pre  <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = high_depth)
  low_post  <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = low_depth)
  high_post <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = high_depth)
  
  points(low_pre$y, low_pre$x, pch = 16, col = "blue", cex = 1.4)
  points(high_pre$y, high_pre$x, pch = 16, col = "red", cex = 1.4)
  points(low_post$y, low_post$x, pch = 17, col = "blue", cex = 1.4)
  points(high_post$y, high_post$x, pch = 17, col = "red", cex = 1.4)
  
  text(low_pre$y, low_pre$x, "Low", pos = 4, col = "blue", cex = 0.85)
  text(high_pre$y, high_pre$x, "High", pos = 4, col = "red", cex = 0.85)
  
  # ----------------------------------------------------
  # 2. Chlorophyll Axis & Layers
  # ----------------------------------------------------
  par(new = TRUE)
  
  plot(
    mean_profile_pre$chl_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = rev(depth_range),
    xlim = chl_range
  )
  
  # Chl SD Shading
  polygon(
    c(pre_chl_clean$chl_mean - pre_chl_clean$chl_sd,
      rev(pre_chl_clean$chl_mean + pre_chl_clean$chl_sd)),
    c(pre_chl_clean$depth_bin, rev(pre_chl_clean$depth_bin)),
    col = rgb(0, 0.39, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(post_chl_clean$chl_mean - post_chl_clean$chl_sd,
      rev(post_chl_clean$chl_mean + post_chl_clean$chl_sd)),
    c(post_chl_clean$depth_bin, rev(post_chl_clean$depth_bin)),
    col = rgb(0.13, 0.55, 0.13, 0.08),
    border = NA
  )
  
  # Chl Mean Curves
  lines(mean_profile_pre$chl_mean, mean_profile_pre$depth_bin, col = "darkgreen", lwd = 2, lty = 1)
  lines(mean_profile_post$chl_mean, mean_profile_post$depth_bin, col = "darkgreen", lwd = 2, lty = 2)
  
  axis(3)
  mtext("Mean chlorophyll (μg/L)", side = 3, line = 2)
  
  # ----------------------------------------------------
  # 3. Legend (绝对坐标准确定位，防截断)
  # ----------------------------------------------------
  par(xpd = TRUE)
  
  x_pos <- chl_range[2] + diff(chl_range) * 0.08
  y_pos <- depth_range[1]
  
  legend(
    x = x_pos,
    y = y_pos,
    legend = c(
      "Turbidity Pre",
      "Turbidity Post",
      "Turbidity SD",
      "Chl Pre",
      "Chl Post",
      "Chl SD",
      "Low turbidity",
      "High turbidity"
    ),
    col = c("black", "black", NA, "darkgreen", "darkgreen", NA, "blue", "red"),
    fill = c(NA, NA, rgb(0, 0, 0, 0.15), NA, NA, rgb(0, 0.39, 0, 0.15), NA, NA),
    border = c(NA, NA, NA, NA, NA, NA, NA, NA),
    lty = c(1, 2, NA, 1, 2, NA, NA, NA),
    pch = c(NA, NA, NA, NA, NA, NA, 16, 16),
    lwd = c(2, 2, NA, 2, 2, NA, NA, NA),
    pt.cex = 1.3,
    bty = "n",
    xjust = 0
  )
}

# ----------------------------------------------------
# 4. 执行直接保存到本地图片文件
# ----------------------------------------------------
png(
  filename = "Site_C_CTD_Profile.png", 
  width = 10, 
  height = 8, 
  units = "in", 
  res = 300
)

plot_site(
  mean_profile_pre,
  mean_profile_post,
  low_depth = 32.5,
  high_depth = 65,
  site_name = "Site C (Mud)"
)

dev.off()

library(dplyr)

target_low  <- 32.5
target_high <- 65

# 合并数据
merged_df <- inner_join(
  mean_profile_pre, mean_profile_post,
  by = "depth_bin", suffix = c("_pre", "_post")
)

# 动态找到数据中离 17.5 和 30.5 最近的实际 depth_bin 数值
actual_low_depth  <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_low))]
actual_high_depth <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_high))]

target_p_values <- merged_df %>%
  filter(depth_bin %in% c(actual_low_depth, actual_high_depth)) %>%
  mutate(
    Layer = ifelse(depth_bin == actual_low_depth, "Low turbidity layer", "High turbidity layer"),
    n_pre = 5, n_post = 5,
    
    # Turbidity T-test
    se_turb = sqrt((turb_sd_pre^2 / n_pre) + (turb_sd_post^2 / n_post)),
    t_turb  = (turb_mean_post - turb_mean_pre) / se_turb,
    p_turb  = 2 * pt(abs(t_turb), df = n_pre + n_post - 2, lower.tail = FALSE),
    
    # Chlorophyll T-test
    se_chl  = sqrt((chl_sd_pre^2 / n_pre) + (chl_sd_post^2 / n_post)),
    t_chl   = (chl_mean_post - chl_mean_pre) / se_chl,
    p_chl   = 2 * pt(abs(t_chl), df = n_pre + n_post - 2, lower.tail = FALSE)
  ) %>%
  select(Layer, depth_bin, t_turb, p_turb, t_chl, p_chl) %>%
  mutate(
    turb_sig = case_when(p_turb < 0.001 ~ "***", p_turb < 0.01 ~ "**", p_turb < 0.05 ~ "*", TRUE ~ "ns"),
    chl_sig  = case_when(p_chl < 0.001 ~ "***", p_chl < 0.01 ~ "**", p_chl < 0.05 ~ "*", TRUE ~ "ns")
  )

# 打印结果
print("=== Low & High 深度 Turbidity 与 Chlorophyll T检验结果 ===")
print(target_p_values)

range(mean_profile_pre$turb_mean, na.rm = TRUE)
range(mean_profile_post$turb_mean, na.rm = TRUE)

range(mean_profile_pre$chl_mean, na.rm = TRUE)
range(mean_profile_post$chl_mean, na.rm = TRUE)






#SITE B

files_pre <- c(
  "DY206_CTD032_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD033_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD034_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD035_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD036_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_pre <- process_ctd_group(files_pre)

files_post <- c(
  "DY206_CTD037_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD038_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD039_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD040_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD041_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_post <- process_ctd_group(files_post)


plot_site <- function(mean_profile_pre,
                      mean_profile_post,
                      low_depth,
                      high_depth,
                      site_name){
  
  # ----------------------------------------------------
  # 0. 数据清洗：剔除 NA 并按深度排序（彻底解决 SD 阴影乱线问题）
  # ----------------------------------------------------
  clean_data <- function(df, mean_col, sd_col) {
    df %>%
      filter(!is.na(.data[[mean_col]]), !is.na(.data[[sd_col]])) %>%
      arrange(depth_bin)
  }
  
  pre_turb_clean  <- clean_data(mean_profile_pre, "turb_mean", "turb_sd")
  post_turb_clean <- clean_data(mean_profile_post, "turb_mean", "turb_sd")
  pre_chl_clean   <- clean_data(mean_profile_pre, "chl_mean", "chl_sd")
  post_chl_clean  <- clean_data(mean_profile_post, "chl_mean", "chl_sd")
  
  # Axis ranges
  turb_range <- range(
    c(
      pre_turb_clean$turb_mean - pre_turb_clean$turb_sd,
      pre_turb_clean$turb_mean + pre_turb_clean$turb_sd,
      post_turb_clean$turb_mean - post_turb_clean$turb_sd,
      post_turb_clean$turb_mean + post_turb_clean$turb_sd
    ),
    na.rm = TRUE
  )
  
  chl_range <- range(
    c(
      pre_chl_clean$chl_mean - pre_chl_clean$chl_sd,
      pre_chl_clean$chl_mean + pre_chl_clean$chl_sd,
      post_chl_clean$chl_mean - post_chl_clean$chl_sd,
      post_chl_clean$chl_mean + post_chl_clean$chl_sd
    ),
    na.rm = TRUE
  )
  
  depth_range <- range(
    c(mean_profile_pre$depth_bin, mean_profile_post$depth_bin),
    na.rm = TRUE
  )
  
  # 预留右侧充足边距存放图例
  par(mar = c(5, 5, 6, 13))
  
  # ----------------------------------------------------
  # 1. Turbidity Axis & Layers
  # ----------------------------------------------------
  plot(
    mean_profile_pre$turb_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    ylim = rev(depth_range),
    xlim = turb_range,
    xlab = expression("Mean turbidity ("*m^{-1}~sr^{-1}*")"),
    ylab = "Depth (m)"
  )
  title(main = site_name, line = 4)
  
  # Turbidity SD Shading
  polygon(
    c(pre_turb_clean$turb_mean - pre_turb_clean$turb_sd,
      rev(pre_turb_clean$turb_mean + pre_turb_clean$turb_sd)),
    c(pre_turb_clean$depth_bin, rev(pre_turb_clean$depth_bin)),
    col = rgb(0, 0, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(post_turb_clean$turb_mean - post_turb_clean$turb_sd,
      rev(post_turb_clean$turb_mean + post_turb_clean$turb_sd)),
    c(post_turb_clean$depth_bin, rev(post_turb_clean$depth_bin)),
    col = rgb(0, 0, 0, 0.08),
    border = NA
  )
  
  # Mean Curves
  lines(mean_profile_pre$turb_mean, mean_profile_pre$depth_bin, lwd = 2, lty = 1, col = "black")
  lines(mean_profile_post$turb_mean, mean_profile_post$depth_bin, lwd = 2, lty = 2, col = "black")
  
  # Selected Depths
  low_pre   <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = low_depth)
  high_pre  <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = high_depth)
  low_post  <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = low_depth)
  high_post <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = high_depth)
  
  points(low_pre$y, low_pre$x, pch = 16, col = "blue", cex = 1.4)
  points(high_pre$y, high_pre$x, pch = 16, col = "red", cex = 1.4)
  points(low_post$y, low_post$x, pch = 17, col = "blue", cex = 1.4)
  points(high_post$y, high_post$x, pch = 17, col = "red", cex = 1.4)
  
  text(low_pre$y, low_pre$x, "Low", pos = 4, col = "blue", cex = 0.85)
  text(high_pre$y, high_pre$x, "High", pos = 4, col = "red", cex = 0.85)
  
  # ----------------------------------------------------
  # 2. Chlorophyll Axis & Layers
  # ----------------------------------------------------
  par(new = TRUE)
  
  plot(
    mean_profile_pre$chl_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = rev(depth_range),
    xlim = chl_range
  )
  
  # Chl SD Shading
  polygon(
    c(pre_chl_clean$chl_mean - pre_chl_clean$chl_sd,
      rev(pre_chl_clean$chl_mean + pre_chl_clean$chl_sd)),
    c(pre_chl_clean$depth_bin, rev(pre_chl_clean$depth_bin)),
    col = rgb(0, 0.39, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(post_chl_clean$chl_mean - post_chl_clean$chl_sd,
      rev(post_chl_clean$chl_mean + post_chl_clean$chl_sd)),
    c(post_chl_clean$depth_bin, rev(post_chl_clean$depth_bin)),
    col = rgb(0.13, 0.55, 0.13, 0.08),
    border = NA
  )
  
  # Chl Mean Curves
  lines(mean_profile_pre$chl_mean, mean_profile_pre$depth_bin, col = "darkgreen", lwd = 2, lty = 1)
  lines(mean_profile_post$chl_mean, mean_profile_post$depth_bin, col = "darkgreen", lwd = 2, lty = 2)
  
  axis(3)
  mtext("Mean chlorophyll (μg/L)", side = 3, line = 2)
  
  # ----------------------------------------------------
  # 3. Legend (绝对坐标准确定位，防截断)
  # ----------------------------------------------------
  par(xpd = TRUE)
  
  x_pos <- chl_range[2] + diff(chl_range) * 0.08
  y_pos <- depth_range[1]
  
  legend(
    x = x_pos,
    y = y_pos,
    legend = c(
      "Turbidity Pre",
      "Turbidity Post",
      "Turbidity SD",
      "Chl Pre",
      "Chl Post",
      "Chl SD",
      "Low turbidity",
      "High turbidity"
    ),
    col = c("black", "black", NA, "darkgreen", "darkgreen", NA, "blue", "red"),
    fill = c(NA, NA, rgb(0, 0, 0, 0.15), NA, NA, rgb(0, 0.39, 0, 0.15), NA, NA),
    border = c(NA, NA, NA, NA, NA, NA, NA, NA),
    lty = c(1, 2, NA, 1, 2, NA, NA, NA),
    pch = c(NA, NA, NA, NA, NA, NA, 16, 16),
    lwd = c(2, 2, NA, 2, 2, NA, NA, NA),
    pt.cex = 1.3,
    bty = "n",
    xjust = 0
  )
}

# ----------------------------------------------------
# 4. 执行直接保存到本地图片文件
# ----------------------------------------------------
png(
  filename = "Site_B_CTD_Profile.png", 
  width = 10, 
  height = 8, 
  units = "in", 
  res = 300
)

plot_site(
  mean_profile_pre,
  mean_profile_post,
  low_depth = 47.5,
  high_depth = 105,
  site_name = "Site B (Mud-muddy sand)"
)

dev.off()

library(dplyr)

target_low  <- 47.5
target_high <- 105

# 合并数据
merged_df <- inner_join(
  mean_profile_pre, mean_profile_post,
  by = "depth_bin", suffix = c("_pre", "_post")
)

# 动态找到数据中离 17.5 和 30.5 最近的实际 depth_bin 数值
actual_low_depth  <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_low))]
actual_high_depth <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_high))]

target_p_values <- merged_df %>%
  filter(depth_bin %in% c(actual_low_depth, actual_high_depth)) %>%
  mutate(
    Layer = ifelse(depth_bin == actual_low_depth, "Low turbidity layer", "High turbidity layer"),
    n_pre = 5, n_post = 5,
    
    # Turbidity T-test
    se_turb = sqrt((turb_sd_pre^2 / n_pre) + (turb_sd_post^2 / n_post)),
    t_turb  = (turb_mean_post - turb_mean_pre) / se_turb,
    p_turb  = 2 * pt(abs(t_turb), df = n_pre + n_post - 2, lower.tail = FALSE),
    
    # Chlorophyll T-test
    se_chl  = sqrt((chl_sd_pre^2 / n_pre) + (chl_sd_post^2 / n_post)),
    t_chl   = (chl_mean_post - chl_mean_pre) / se_chl,
    p_chl   = 2 * pt(abs(t_chl), df = n_pre + n_post - 2, lower.tail = FALSE)
  ) %>%
  select(Layer, depth_bin, t_turb, p_turb, t_chl, p_chl) %>%
  mutate(
    turb_sig = case_when(p_turb < 0.001 ~ "***", p_turb < 0.01 ~ "**", p_turb < 0.05 ~ "*", TRUE ~ "ns"),
    chl_sig  = case_when(p_chl < 0.001 ~ "***", p_chl < 0.01 ~ "**", p_chl < 0.05 ~ "*", TRUE ~ "ns")
  )

# 打印结果
print("=== Low & High 深度 Turbidity 与 Chlorophyll T检验结果 ===")
print(target_p_values)

range(mean_profile_pre$turb_mean, na.rm = TRUE)
range(mean_profile_post$turb_mean, na.rm = TRUE)

range(mean_profile_pre$chl_mean, na.rm = TRUE)
range(mean_profile_post$chl_mean, na.rm = TRUE)










#site A
files_pre <- c(
  "DY206_CTD073_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD068_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD070_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD071_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD072_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_pre <- process_ctd_group(files_pre)

files_post <- c(
  "DY206_CTD074_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD075_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD076_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD077_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD078_align_CTM_Derive_2Hz_Strip.cnv"
)

mean_profile_post <- process_ctd_group(files_post)


plot_site <- function(mean_profile_pre,
                      mean_profile_post,
                      low_depth,
                      high_depth,
                      site_name){
  
  # ----------------------------------------------------
  # 0. 数据清洗：剔除 NA 并按深度排序（彻底解决 SD 阴影乱线问题）
  # ----------------------------------------------------
  clean_data <- function(df, mean_col, sd_col) {
    df %>%
      filter(!is.na(.data[[mean_col]]), !is.na(.data[[sd_col]])) %>%
      arrange(depth_bin)
  }
  
  pre_turb_clean  <- clean_data(mean_profile_pre, "turb_mean", "turb_sd")
  post_turb_clean <- clean_data(mean_profile_post, "turb_mean", "turb_sd")
  pre_chl_clean   <- clean_data(mean_profile_pre, "chl_mean", "chl_sd")
  post_chl_clean  <- clean_data(mean_profile_post, "chl_mean", "chl_sd")
  
  # Axis ranges
  turb_range <- range(
    c(
      pre_turb_clean$turb_mean - pre_turb_clean$turb_sd,
      pre_turb_clean$turb_mean + pre_turb_clean$turb_sd,
      post_turb_clean$turb_mean - post_turb_clean$turb_sd,
      post_turb_clean$turb_mean + post_turb_clean$turb_sd
    ),
    na.rm = TRUE
  )
  
  chl_range <- range(
    c(
      pre_chl_clean$chl_mean - pre_chl_clean$chl_sd,
      pre_chl_clean$chl_mean + pre_chl_clean$chl_sd,
      post_chl_clean$chl_mean - post_chl_clean$chl_sd,
      post_chl_clean$chl_mean + post_chl_clean$chl_sd
    ),
    na.rm = TRUE
  )
  
  depth_range <- range(
    c(mean_profile_pre$depth_bin, mean_profile_post$depth_bin),
    na.rm = TRUE
  )
  
  # 预留右侧充足边距存放图例
  par(mar = c(5, 5, 6, 13))
  
  # ----------------------------------------------------
  # 1. Turbidity Axis & Layers
  # ----------------------------------------------------
  plot(
    mean_profile_pre$turb_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    ylim = rev(depth_range),
    xlim = turb_range,
    xlab = expression("Mean turbidity ("*m^{-1}~sr^{-1}*")"),
    ylab = "Depth (m)"
  )
  title(main = site_name, line = 4)
  
  # Turbidity SD Shading
  polygon(
    c(pre_turb_clean$turb_mean - pre_turb_clean$turb_sd,
      rev(pre_turb_clean$turb_mean + pre_turb_clean$turb_sd)),
    c(pre_turb_clean$depth_bin, rev(pre_turb_clean$depth_bin)),
    col = rgb(0, 0, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(post_turb_clean$turb_mean - post_turb_clean$turb_sd,
      rev(post_turb_clean$turb_mean + post_turb_clean$turb_sd)),
    c(post_turb_clean$depth_bin, rev(post_turb_clean$depth_bin)),
    col = rgb(0, 0, 0, 0.08),
    border = NA
  )
  
  # Mean Curves
  lines(mean_profile_pre$turb_mean, mean_profile_pre$depth_bin, lwd = 2, lty = 1, col = "black")
  lines(mean_profile_post$turb_mean, mean_profile_post$depth_bin, lwd = 2, lty = 2, col = "black")
  
  # Selected Depths
  low_pre   <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = low_depth)
  high_pre  <- approx(x = mean_profile_pre$depth_bin, y = mean_profile_pre$turb_mean, xout = high_depth)
  low_post  <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = low_depth)
  high_post <- approx(x = mean_profile_post$depth_bin, y = mean_profile_post$turb_mean, xout = high_depth)
  
  points(low_pre$y, low_pre$x, pch = 16, col = "blue", cex = 1.4)
  points(high_pre$y, high_pre$x, pch = 16, col = "red", cex = 1.4)
  points(low_post$y, low_post$x, pch = 17, col = "blue", cex = 1.4)
  points(high_post$y, high_post$x, pch = 17, col = "red", cex = 1.4)
  
  text(low_pre$y, low_pre$x, "Low", pos = 4, col = "blue", cex = 0.85)
  text(high_pre$y, high_pre$x, "High", pos = 4, col = "red", cex = 0.85)
  
  # ----------------------------------------------------
  # 2. Chlorophyll Axis & Layers
  # ----------------------------------------------------
  par(new = TRUE)
  
  plot(
    mean_profile_pre$chl_mean,
    mean_profile_pre$depth_bin,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = rev(depth_range),
    xlim = chl_range
  )
  
  # Chl SD Shading
  polygon(
    c(pre_chl_clean$chl_mean - pre_chl_clean$chl_sd,
      rev(pre_chl_clean$chl_mean + pre_chl_clean$chl_sd)),
    c(pre_chl_clean$depth_bin, rev(pre_chl_clean$depth_bin)),
    col = rgb(0, 0.39, 0, 0.12),
    border = NA
  )
  
  polygon(
    c(post_chl_clean$chl_mean - post_chl_clean$chl_sd,
      rev(post_chl_clean$chl_mean + post_chl_clean$chl_sd)),
    c(post_chl_clean$depth_bin, rev(post_chl_clean$depth_bin)),
    col = rgb(0.13, 0.55, 0.13, 0.08),
    border = NA
  )
  
  # Chl Mean Curves
  lines(mean_profile_pre$chl_mean, mean_profile_pre$depth_bin, col = "darkgreen", lwd = 2, lty = 1)
  lines(mean_profile_post$chl_mean, mean_profile_post$depth_bin, col = "darkgreen", lwd = 2, lty = 2)
  
  axis(3)
  mtext("Mean chlorophyll (μg/L)", side = 3, line = 2)
  
  # ----------------------------------------------------
  # 3. Legend (绝对坐标准确定位，防截断)
  # ----------------------------------------------------
  par(xpd = TRUE)
  
  x_pos <- chl_range[2] + diff(chl_range) * 0.08
  y_pos <- depth_range[1]
  
  legend(
    x = x_pos,
    y = y_pos,
    legend = c(
      "Turbidity Pre",
      "Turbidity Post",
      "Turbidity SD",
      "Chl Pre",
      "Chl Post",
      "Chl SD",
      "Low turbidity",
      "High turbidity"
    ),
    col = c("black", "black", NA, "darkgreen", "darkgreen", NA, "blue", "red"),
    fill = c(NA, NA, rgb(0, 0, 0, 0.15), NA, NA, rgb(0, 0.39, 0, 0.15), NA, NA),
    border = c(NA, NA, NA, NA, NA, NA, NA, NA),
    lty = c(1, 2, NA, 1, 2, NA, NA, NA),
    pch = c(NA, NA, NA, NA, NA, NA, 16, 16),
    lwd = c(2, 2, NA, 2, 2, NA, NA, NA),
    pt.cex = 1.3,
    bty = "n",
    xjust = 0
  )
}

# ----------------------------------------------------
# 4. 执行直接保存到本地图片文件
# ----------------------------------------------------
png(
  filename = "Site_A_CTD_Profile.png", 
  width = 10, 
  height = 8, 
  units = "in", 
  res = 300
)

plot_site(
  mean_profile_pre,
  mean_profile_post,
  low_depth = 17.5,
  high_depth = 30.5,
  site_name = "Site A (Muddy sand)"
)

dev.off()

# ====================================================
# 4. 【更稳定版】计算 Low & High 深度的 p 值表格
# ====================================================
library(dplyr)

target_low  <- 17.5
target_high <- 30.5

# 合并数据
merged_df <- inner_join(
  mean_profile_pre, mean_profile_post,
  by = "depth_bin", suffix = c("_pre", "_post")
)

# 动态找到数据中离 17.5 和 30.5 最近的实际 depth_bin 数值
actual_low_depth  <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_low))]
actual_high_depth <- merged_df$depth_bin[which.min(abs(merged_df$depth_bin - target_high))]

target_p_values <- merged_df %>%
  filter(depth_bin %in% c(actual_low_depth, actual_high_depth)) %>%
  mutate(
    Layer = ifelse(depth_bin == actual_low_depth, "Low turbidity layer", "High turbidity layer"),
    n_pre = 5, n_post = 5,
    
    # Turbidity T-test
    se_turb = sqrt((turb_sd_pre^2 / n_pre) + (turb_sd_post^2 / n_post)),
    t_turb  = (turb_mean_post - turb_mean_pre) / se_turb,
    p_turb  = 2 * pt(abs(t_turb), df = n_pre + n_post - 2, lower.tail = FALSE),
    
    # Chlorophyll T-test
    se_chl  = sqrt((chl_sd_pre^2 / n_pre) + (chl_sd_post^2 / n_post)),
    t_chl   = (chl_mean_post - chl_mean_pre) / se_chl,
    p_chl   = 2 * pt(abs(t_chl), df = n_pre + n_post - 2, lower.tail = FALSE)
  ) %>%
  select(Layer, depth_bin, t_turb, p_turb, t_chl, p_chl) %>%
  mutate(
    turb_sig = case_when(p_turb < 0.001 ~ "***", p_turb < 0.01 ~ "**", p_turb < 0.05 ~ "*", TRUE ~ "ns"),
    chl_sig  = case_when(p_chl < 0.001 ~ "***", p_chl < 0.01 ~ "**", p_chl < 0.05 ~ "*", TRUE ~ "ns")
  )

# 打印结果
print("=== Low & High 深度 Turbidity 与 Chlorophyll T检验结果 ===")
print(target_p_values)

range(mean_profile_pre$turb_mean, na.rm = TRUE)
range(mean_profile_post$turb_mean, na.rm = TRUE)

range(mean_profile_pre$chl_mean, na.rm = TRUE)
range(mean_profile_post$chl_mean, na.rm = TRUE)

