#背景指标图
library(oce)
library(dplyr)

png(
  filename = "SiteB_Temperature_Salinity.png",
  width = 1800,
  height = 2600,
  res = 300
)

setwd("/Users/langjiawen/Desktop/DY206_CTD008_align_CTM_Derive_2Hz_Strip")

############################################################
## Function
############################################################

process_ctd_group <- function(files){
  
  all_data <- list()
  
  for(i in seq_along(files)){
    
    ctd <- read.ctd(files[i])
    
    df <- data.frame(
      
      depth = ctd[["depth"]],
      temp  = ctd[["temperature"]],
      sal   = ctd[["salinity"]]
      
    )
    
    df <- df %>%
      filter(complete.cases(.)) %>%
      filter(depth > 2)
    
    df$depth_bin <- floor(df$depth/2)*2
    
    df_bin <- df %>%
      group_by(depth_bin) %>%
      summarise(
        
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
      
      temp_mean = mean(temp),
      sal_mean  = mean(sal),
      
      .groups = "drop"
      
    )
  
  return(summary_profile)
  
}

############################################################
## Site B Pre-trawl CTDs
############################################################

files_pre <- c(
  
  "DY206_CTD032_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD033_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD034_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD035_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD036_align_CTM_Derive_2Hz_Strip.cnv"
  
)

############################################################
## Process
############################################################

profile <- process_ctd_group(files_pre)

profile %>%
  summarise(
    temp_min = min(temp_mean, na.rm = TRUE),
    temp_max = max(temp_mean, na.rm = TRUE),
    sal_min  = min(sal_mean, na.rm = TRUE),
    sal_max  = max(sal_mean, na.rm = TRUE)
  )


############################################################
## Axis ranges
############################################################

depth_range <- range(profile$depth_bin)

temp_range <- range(profile$temp_mean)

sal_range <- range(profile$sal_mean)

############################################################
## Plot
############################################################

par(mar = c(5,5,7,5))


plot(
  
  profile$temp_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "red",
  
  ylim = rev(depth_range),
  xlim = temp_range,
  
  xlab = "Temperature (°C)",
  ylab = "Depth (m)",
  
  main = "Site B"
  
)

## Salinity on top axis

par(new = TRUE)

plot(
  
  profile$sal_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "blue",
  
  axes = FALSE,
  xlab = "",
  ylab = "",
  
  ylim = rev(depth_range),
  xlim = sal_range
  
)

axis(3)

mtext(
  
  "Salinity (PSU)",
  
  side = 3,
  
  line = 2
  
)

############################################################
## Legend
############################################################

legend(
  
  "bottomleft",
  
  legend = c(
    "Temperature",
    "Salinity"
  ),
  
  col = c(
    "red",
    "blue"
  ),
  
  lwd = 2,
  
  bty = "n"
  
)

dev.off()




#Sice c
png(
  filename = "SiteC_Temperature_Salinity.png",
  width = 1800,
  height = 2600,
  res = 300
)

setwd("/Users/langjiawen/Desktop/DY206_CTD008_align_CTM_Derive_2Hz_Strip")

############################################################
## Function
############################################################

process_ctd_group <- function(files){
  
  all_data <- list()
  
  for(i in seq_along(files)){
    
    ctd <- read.ctd(files[i])
    
    df <- data.frame(
      
      depth = ctd[["depth"]],
      temp  = ctd[["temperature"]],
      sal   = ctd[["salinity"]]
      
    )
    
    df <- df %>%
      filter(complete.cases(.)) %>%
      filter(depth > 2)
    
    df$depth_bin <- floor(df$depth/2)*2
    
    df_bin <- df %>%
      group_by(depth_bin) %>%
      summarise(
        
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
      
      temp_mean = mean(temp),
      sal_mean  = mean(sal),
      
      .groups = "drop"
      
    )
  
  return(summary_profile)
  
}

############################################################

files_pre <- c(
  
  "DY206_CTD002_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD003_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD004_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD005_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD006_align_CTM_Derive_2Hz_Strip.cnv"
  
)

############################################################

profile <- process_ctd_group(files_pre)

profile %>%
  summarise(
    temp_min = min(temp_mean, na.rm = TRUE),
    temp_max = max(temp_mean, na.rm = TRUE),
    sal_min  = min(sal_mean, na.rm = TRUE),
    sal_max  = max(sal_mean, na.rm = TRUE)
  )

############################################################
## Axis ranges
############################################################

depth_range <- range(profile$depth_bin)

temp_range <- range(profile$temp_mean)

sal_range <- range(profile$sal_mean)

############################################################
## Plot
############################################################

par(mar = c(5,5,7,5))


plot(
  
  profile$temp_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "red",
  
  ylim = rev(depth_range),
  xlim = temp_range,
  
  xlab = "Temperature (°C)",
  ylab = "Depth (m)",
  
  main = "Site C"
  
)

## Salinity on top axis

par(new = TRUE)

plot(
  
  profile$sal_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "blue",
  
  axes = FALSE,
  xlab = "",
  ylab = "",
  
  ylim = rev(depth_range),
  xlim = sal_range
  
)

axis(3)

mtext(
  
  "Salinity (PSU)",
  
  side = 3,
  
  line = 2
  
)

############################################################
## Legend
############################################################

legend(
  
  "bottomleft",
  
  legend = c(
    "Temperature",
    "Salinity"
  ),
  
  col = c(
    "red",
    "blue"
  ),
  
  lwd = 2,
  
  bty = "n"
  
)

dev.off()


#Sice D
png(
  filename = "SiteD_Temperature_Salinity.png",
  width = 1800,
  height = 2600,
  res = 300
)

setwd("/Users/langjiawen/Desktop/DY206_CTD008_align_CTM_Derive_2Hz_Strip")

############################################################
## Function
############################################################

process_ctd_group <- function(files){
  
  all_data <- list()
  
  for(i in seq_along(files)){
    
    ctd <- read.ctd(files[i])
    
    df <- data.frame(
      
      depth = ctd[["depth"]],
      temp  = ctd[["temperature"]],
      sal   = ctd[["salinity"]]
      
    )
    
    df <- df %>%
      filter(complete.cases(.)) %>%
      filter(depth > 2)
    
    df$depth_bin <- floor(df$depth/2)*2
    
    df_bin <- df %>%
      group_by(depth_bin) %>%
      summarise(
        
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
      
      temp_mean = mean(temp),
      sal_mean  = mean(sal),
      
      .groups = "drop"
      
    )
  
  return(summary_profile)
  
}

############################################################

files_pre <- c(
  
  "DY206_CTD087_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD088_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD089_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD090_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD091_align_CTM_Derive_2Hz_Strip.cnv"
  
)

############################################################

profile <- process_ctd_group(files_pre)

profile %>%
  summarise(
    temp_min = min(temp_mean, na.rm = TRUE),
    temp_max = max(temp_mean, na.rm = TRUE),
    sal_min  = min(sal_mean, na.rm = TRUE),
    sal_max  = max(sal_mean, na.rm = TRUE)
  )

############################################################
## Axis ranges
############################################################

depth_range <- range(profile$depth_bin)

temp_range <- range(profile$temp_mean)

sal_range <- range(profile$sal_mean)

############################################################
## Plot
############################################################

par(mar = c(5,5,7,5))


plot(
  
  profile$temp_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "red",
  
  ylim = rev(depth_range),
  xlim = temp_range,
  
  xlab = "Temperature (°C)",
  ylab = "Depth (m)",
  
  main = "Site D"
  
)

## Salinity on top axis

par(new = TRUE)

plot(
  
  profile$sal_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "blue",
  
  axes = FALSE,
  xlab = "",
  ylab = "",
  
  ylim = rev(depth_range),
  xlim = sal_range
  
)

axis(3)

mtext(
  
  "Salinity (PSU)",
  
  side = 3,
  
  line = 2
  
)

############################################################
## Legend
############################################################

legend(
  
  "bottomleft",
  
  legend = c(
    "Temperature",
    "Salinity"
  ),
  
  col = c(
    "red",
    "blue"
  ),
  
  lwd = 2,
  
  bty = "n"
  
)

dev.off()



#Sice A
png(
  filename = "SiteA_Temperature_Salinity.png",
  width = 1800,
  height = 2600,
  res = 300
)

setwd("/Users/langjiawen/Desktop/DY206_CTD008_align_CTM_Derive_2Hz_Strip")

############################################################
## Function
############################################################

process_ctd_group <- function(files){
  
  all_data <- list()
  
  for(i in seq_along(files)){
    
    ctd <- read.ctd(files[i])
    
    df <- data.frame(
      
      depth = ctd[["depth"]],
      temp  = ctd[["temperature"]],
      sal   = ctd[["salinity"]]
      
    )
    
    df <- df %>%
      filter(complete.cases(.)) %>%
      filter(depth > 2)
    
    df$depth_bin <- floor(df$depth/2)*2
    
    df_bin <- df %>%
      group_by(depth_bin) %>%
      summarise(
        
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
      
      temp_mean = mean(temp),
      sal_mean  = mean(sal),
      
      .groups = "drop"
      
    )
  
  return(summary_profile)
  
}

############################################################

files_pre <- c(
  
  "DY206_CTD068_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD069_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD070_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD071_align_CTM_Derive_2Hz_Strip.cnv",
  "DY206_CTD072_align_CTM_Derive_2Hz_Strip.cnv"
  
)

############################################################

profile <- process_ctd_group(files_pre)

profile %>%
  summarise(
    temp_min = min(temp_mean, na.rm = TRUE),
    temp_max = max(temp_mean, na.rm = TRUE),
    sal_min  = min(sal_mean, na.rm = TRUE),
    sal_max  = max(sal_mean, na.rm = TRUE)
  )

############################################################
## Axis ranges
############################################################

depth_range <- range(profile$depth_bin)

temp_range <- range(profile$temp_mean)

sal_range <- range(profile$sal_mean)

############################################################
## Plot
############################################################

par(mar = c(5,5,7,5))


plot(
  
  profile$temp_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "red",
  
  ylim = rev(depth_range),
  xlim = temp_range,
  
  xlab = "Temperature (°C)",
  ylab = "Depth (m)",
  
  main = "Site A"
  
)

## Salinity on top axis

par(new = TRUE)

plot(
  
  profile$sal_mean,
  profile$depth_bin,
  
  type = "l",
  lwd = 2,
  col = "blue",
  
  axes = FALSE,
  xlab = "",
  ylab = "",
  
  ylim = rev(depth_range),
  xlim = sal_range
  
)

axis(3)

mtext(
  
  "Salinity (PSU)",
  
  side = 3,
  
  line = 2
  
)

############################################################
## Legend
############################################################

legend(
  
  "bottomleft",
  
  legend = c(
    "Temperature",
    "Salinity"
  ),
  
  col = c(
    "red",
    "blue"
  ),
  
  lwd = 2,
  
  bty = "n"
  
)

dev.off()





ctd <- read.ctd("DY206_CTD087_align_CTM_Derive_2Hz_Strip.cnv")
max(ctd[["depth"]], na.rm = TRUE)
ctd <- read.ctd("DY206_CTD068_align_CTM_Derive_2Hz_Strip.cnv")
max(ctd[["depth"]], na.rm = TRUE)
ctd <- read.ctd("DY206_CTD002_align_CTM_Derive_2Hz_Strip.cnv")
max(ctd[["depth"]], na.rm = TRUE)
ctd <- read.ctd("DY206_CTD032_align_CTM_Derive_2Hz_Strip.cnv")
max(ctd[["depth"]], na.rm = TRUE)
