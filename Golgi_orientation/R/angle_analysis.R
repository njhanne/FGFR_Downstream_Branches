#### 0.0 Load libraries for packages ####
# Uncomment and run this the first time to make sure everything gets installed
# install.packages(c('circular', 'plyr', 'dplyr', 'tidyr', 'stringr', 'RImageJROI', 'sf', 'terra', 'ggplot2', 'RColorBrewer', 'Hmisc', 'svglite', 'magick'))
# if terra takes more than 60s to download it will fail
# run this line: options(timeout = max(1000, getOption("timeout")))

# circle stats stuff
# library(bpnreg) # not using this as it can't handle size of our data :(
library(circular)

# df assistance
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

# needed for masking out unwanted regions
library(RImageJROI)
library(sf)

# needed for 'globe' plots
library(terra)

# graphing
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(svglite)
library(magick)


#### 0.1 Helper functions ####
#for making the dataframe and correcting angle data
compile_angle_df <- function() {
  filenames <- dir("./cellpose/angle_results/", "*._angle_results.csv") # get all output csv from python
  types <- sub('(\\A*)_angle_results.*', '\\1', filenames) # get the sample names
  # combine all the individual csv into one big dataframe
  # this can be quite slow as there are tons of data
  df<- adply(data.frame(f=I(filenames), t=types), 1,
             with, cbind(read.csv(file.path('./cellpose/angle_results', f)), sample_info=t))
  # splits up the filename into relevant sample info and puts them in their own little columns
  split_text <- str_split(df$sample_info, '_')
  df <- df %>% mutate(treatment = sapply(split_text, function(l) l[[1]]))
  df <- df %>% mutate(id = sapply(split_text, function(l) l[[2]]))
  df <- df %>% mutate(side = sapply(split_text, function(l) l[[3]]))
  df <- df %>% mutate(section = sapply(split_text, function(l) l[[4]]))
  df <- df %>% mutate(view = sapply(split_text, function(l) l[[5]]))
  df <- df %>% mutate(side_num = case_when(side == "control" ~ 0,
                                           TRUE ~ 1))
  df <- df %>% mutate(sample_info = tolower(sample_info))
  
  df$id_old <- df$id
  df <- df %>% unite(id, c('treatment','id_old'), remove=FALSE)
  df$id <- as.factor(df$id)
  df$id <- as.numeric(df$id)
  
  df <- df %>% select(!X)  # deletes unneeded column
  return(df)
} 


adjust_base_angles <- function(df, info_file){
  angle_adjustment_df <- info_file %>% filter(!is.na(angle_adjustment)) %>% select(new_filename, angle_adjustment) %>% distinct()
  angle_adjustment_df <- angle_adjustment_df %>% mutate(angle_adjustment = (angle_adjustment * pi)/180) # convert to radians
  df <- left_join(df, angle_adjustment_df, by=c('t' = 'new_filename')) # put the adjustment into our main df
  df$angle_adjustment[is.na(df$angle_adjustment)] <- 0 # change na's to zero
  df$angle_old <- df$angle # just in case we save the og angles
  df <- df %>% mutate(angle = angle - angle_adjustment) # 'add' the new angle adjustment
  df <- df %>% mutate(angle = case_when(angle > 2*pi ~ angle - 2*pi,
                                        angle < 0 ~ angle + 2*pi,
                                        TRUE ~ angle))
  return(df)
}


apply_mask <- function(filtered_df, imagej_mask, type_column, save_img) {
  centroids <- filtered_df %>% select(nuclei_centroidx, nuclei_centroidy)
  mask_linestring <- st_linestring(imagej_mask$coords, dim='XY')
  mask_polygon <- st_cast(mask_linestring, 'POLYGON')
  keep <- st_intersection(st_multipoint(data.matrix(centroids), dim='XY'), mask_polygon)
  
  if (save_img) {
    p <- ggplot() + geom_sf(data = mask_polygon) +
      geom_point(data = centroids, aes(nuclei_centroidx, nuclei_centroidy), color = 'red') + geom_sf(data=keep)
    file_name <- paste(type_column, "_masked_nuclei.svg")
    ggsave(filename= file.path('./analysis_output/mask_images', file_name), p)
  }

  keep <- as.data.frame(as.matrix(keep))
  filtered_df <- inner_join(filtered_df, keep, by=c('nuclei_centroidx' = 'V1', 'nuclei_centroidy' = 'V2'))
  return(filtered_df)
}


mask_out_centroids <- function(df, mask_filenames, save_img=TRUE) {
  df_masked <- data.frame()
  for (i in 1:length(mask_filenames)) {
    # reading in the imagej roi creates a matrix of xy coordinates which can be plugged into sf
    mask <- read.ijroi(file.path("./imagej_rois/region_mask_rois/", mask_filenames[i]), verbose = FALSE)
    type_col <- sub('(\\A*).roi$', '\\1', basename(mask_filenames[i])) # get column name for matching
    temp_df <- df %>% filter(t == type_col)
    if (nrow(temp_df) == 0) {
      # if the python output data isn't there then don't analyze it...
      print(paste0('Skipping ', type_col))
      filtered_df <- temp_df
      df_masked <- rbind(df_masked, filtered_df)
    } else {
      # apply the mask to exclude nuclei outside of it
      filtered_df <- apply_mask(temp_df, mask, type_col, save_img)
      df_masked <- rbind(df_masked, filtered_df) # probably not supposed to do this in a loop but it works
    }
  }
  df_masked <- rbind(df_masked, df %>% filter(df$t %in% setdiff(unique(df$t), unique(df_masked$t))))
  return(df_masked)
}


apply_distance_filter <- function(filtered_df, baseline_roi, type_column, thresh_dist, save_img) {
  centroids <- filtered_df %>% select(nuclei_centroidx, nuclei_centroidy)
  baseline_linestring <- st_linestring(baseline_roi$coords, dim='XY')
  keep <- st_intersection(st_multipoint(data.matrix(centroids), dim='XY'), st_buffer(baseline_linestring, thresh_dist/0.2840910))
  
  if (save_img) {
    p <- ggplot() + geom_point(data = centroids, aes(nuclei_centroidx, nuclei_centroidy), color = 'purple') + geom_sf(data=keep)
    file_name <- paste(type_column, "_masked_nuclei2.svg")
    ggsave(filename= file.path('./analysis_output/mask_images', file_name), p)
  }
  
  keep <- as.data.frame(as.matrix(keep))
  filtered_df <- inner_join(filtered_df, keep, by=c('nuclei_centroidx' = 'V1', 'nuclei_centroidy' = 'V2'))
  return(filtered_df)
}


filter_baseline_distance <- function(baseline_mask_filenames, df, distance=200, save_img=TRUE) {
  df_masked <- data.frame()
  for (i in 1:length(baseline_mask_filenames)) {
    # reading in the imagej roi creates a matrix of xy coordinates which can be plugged into sf
    baseline_roi <- read.ijroi(file.path("./imagej_rois/baseline_rois", baseline_mask_filenames[i]), verbose = FALSE)
    type_col <- sub('(\\A*)_baseline.roi$', '\\1', basename(baseline_mask_filenames[i])) # get column name for matching
    temp_df <- df %>% filter(t == type_col)
    if (nrow(temp_df) == 0) {
      # if the python output data isn't there then don't analyze it...
      print(paste0('Skipping ', type_col))
      filtered_df <- temp_df
      df_masked <- rbind(df_masked, filtered_df)
    } else {
      # apply the mask to exclude nuclei outside of it
      filtered_df <- apply_distance_filter(temp_df, baseline_roi, type_col, distance, save_img)
      df_masked <- rbind(df_masked, filtered_df) # probably not supposed to do this in a loop but it works
    }
  }
  df_masked <- rbind(df_masked, df %>% filter(df$t %in% setdiff(unique(df$t), unique(df_masked$t))))
  return(df_masked)
}


flip_y_angles <- function(df_to_flip) {
  df_flipped_control <- df_to_flip %>% filter((treatment == 'LY' & side_num == 1) |  (treatment != 'LY' & side_num == 0))
  df_flipped_control <- df_flipped_control %>% mutate(angle = case_when(angle <= pi ~ pi-angle, angle > pi ~ 3*pi - angle))
  if ('angle_old' %in% names(df_flipped_control)) {
    df_flipped_control <- df_flipped_control %>% mutate(angle_old = case_when(angle_old <= pi ~ pi-angle_old, angle_old > pi ~ 3*pi - angle_old))
  }
  
  df_original_control <- df_to_flip %>% filter((treatment != 'LY' & side_num == 1) |  (treatment == 'LY' & side_num == 0))
  df_flipped_control <- rbind(df_flipped_control, df_original_control)
  return(df_flipped_control)
}


#### 0.2 Cellularity helpers ####
mask_inter_area <- function(area_mask, baseline_roi, thresh_dist, res) {
  baseline_linestring <- st_linestring(baseline_roi$coords, dim='XY')
  baseline_buffer <- st_buffer(baseline_linestring, thresh_dist/res)
  
  mask_linestring <- st_linestring(area_mask$coords, dim='XY')
  mask_polygon <- st_cast(mask_linestring, 'POLYGON')
  if (st_is_valid(mask_polygon))  {
    area <- st_area(st_intersection(mask_polygon, baseline_buffer)) * (res^2)
  } else {
    area <- st_area(st_intersection(st_buffer(mask_polygon, 0), baseline_buffer)) * (res^2)
  }
  return(area)
}


get_cellularity <- function(df, matches, mask_filenames, baseline_mask_filenames, distance=200, res=0.2840910) {
  cellularity <- data.frame()
  for (i in 1:nrow(matches)) {
    # reading in the imagej roi creates a matrix of xy coordinates which can be plugged into sf
    j <- matches[i,2]
    area_mask <- read.ijroi(file.path("./imagej_rois/region_mask_rois", mask_filenames[j]), verbose = FALSE)
    type_col <- sub('(\\A*).roi$', '\\1', basename(mask_filenames[j])) # get column name for matching
    # cellularity <- mask_cellularity(df %>% filter(t == type_col), mask_area, type_col)
    
    baseline_roi <- read.ijroi(file.path("./imagej_rois/baseline_rois", baseline_mask_filenames[i]), verbose = FALSE)
    # type_col <- sub('(\\A*)_baseline.roi$', '\\1', basename(baseline_mask_filenames[i])) # get column name for matching
    area <- mask_inter_area(area_mask, baseline_roi, distance, res)
    cells <- nrow(df %>% filter(t == type_col))
    if (cells != 0) {
      cellularity <- rbind(cellularity, c(type_col, area, cells, cells/area)) # probably not supposed to do this in a loop but it works
    }
  }
  colnames(cellularity) <- c('t', 'area', 'cells', 'cells_per_sqmicron')
  df <- full_join(df, cellularity)
  return(df)
}


#### 0.3 Circular statistics helpers ####
get_circular_stats <- function(treatment_name, side_name, df) {
  temp_df <- df %>% filter(treatment == treatment_name, side == side_name)
  circ_info <- circularize_group(temp_df$angle)
  means <- unlist(group_sample_circular_stats(circ_info[[1]], temp_df$id, 'mean'))
  sds <- unlist(group_sample_circular_stats(circ_info[[1]], temp_df$id, 'sd'))
  rhos <- unlist(group_sample_circular_stats(circ_info[[1]], temp_df$id, 'rho'))
  return(list(treatment_name, side_name,
              mean.circular(means),
              circ_info[[2]],
              sd.circular(means),
              mean.circular(sds),
              circ_info[[3]],
              rho.circular(means),
              mean.circular(rhos),
              circ_info[[4]]))
}


circularize_group <- function(angles) {
  group.circ <- circular(angles, units = 'radians')
  group.circ.mean <- mean.circular(group.circ)
  group.circ.sd <- sd.circular(group.circ)
  group.circ.cr <- rho.circular(group.circ)
  return(list(group.circ, group.circ.mean, group.circ.sd, group.circ.cr))
}


group_sample_circular_stats <- function(angles, ids, stat) {
  samples <- length(unique(ids))
  statistic <- list()
  for (i in 1:samples) {
    samplen <- unique(ids)[i]
    if (stat == 'mean') {
      statistic[[i]] <- mean.circular(angles[ids == samplen])
    }
    if (stat == 'median') {
      statistic[[i]] <- median.circular(angles[ids == samplen])
    }
    if (stat == 'rho') {
      statistic[[i]] <- rho.circular(angles[ids == samplen])
    }
    if (stat == 'sd') {
      statistic[[i]] <- sd.circular(angles[ids == samplen])
    }
  }
  return(statistic)
}


WalraffTest <- function(cdat, ndat, g, gID) {
  # http://circstatinr.st-andrews.ac.uk/resources/Chap7-RCommands.txt
  # tests for homoscedasticity (whether measures are drawn from same population)
  N <- length(cdat) ; ndatcsum <- cumsum(ndat) ; tbar <- circular(0) ; distdat <- 0
  for (k in 1:g) {
  dist <- 0 ; sample <- circular(0)
  if (k==1) {low <- 0} else
  if (k > 1) {low <- ndatcsum[k-1]}
  for (j in 1:ndat[k]) { sample[j] <- cdat[j+low] }
  tm1 <- trigonometric.moment(sample, p=1) ; tbar[k] <- tm1$mu
  for (j in 1:ndat[k]) { dist[j] <- pi-abs(pi-abs(sample[j]-tbar[k])) }
  distdat <- c(distdat, dist)
  }
  distdat <- distdat[-1]
  # gID <- c(rep(1,n1), rep(2,n2), rep(3,n3))
  TestRes <- kruskal.test(distdat, g=gID)
  return(TestRes)
}


#### 0.4 Windrose plot helpers ####
# https://stackoverflow.com/a/17266781
# https://stackoverflow.com/questions/62343603/r-windrose-percent-label-on-figure
plot.windrose <- function(data, dirres = 10, color, control) {
  dir.breaks <- c(-dirres/2,
                  seq(dirres/2, 360-dirres/2, by = dirres),
                  360+dirres/2)
  
  dir.labels <- c(paste(0),
                  paste(seq(dirres, 360 - dirres, by = dirres)),
                  paste(0))
  
  # assign each wind direction to a bin
  dir.binned <- cut(data$angle_deg, breaks = dir.breaks, ordered_result = TRUE)
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned
  # generate bin sums for each samples as % of total cells
  T_data <- data %>% group_by(sample_side, dir.binned) %>% dplyr::summarise(count= n()) %>% mutate(y = count/sum(count))
  # average them
  T_data <- T_data %>% group_by(dir.binned) %>% dplyr::summarize(z = mean(y))
  
  
  labels <- data.frame(x = pi, y = scales::extended_breaks()(range(T_data$z)))
  print(as.character(data$treatment[1]))
  print(max(T_data$z)*100)
  
  if(missing(control)) {
    p.windrose <- ggplot(data = T_data, aes(x = dir.binned, y = z, fill = color, color = color)) +
      geom_bar(width = 1, linewidth = .5, stat='identity') +
      scale_x_discrete(drop = FALSE, labels = waiver()) +
      scale_y_continuous(limits = c(0, 0.042), expand = c(0, 0),  breaks = c(0,.01,.02,.03,.04)) +
      coord_polar(start = ((270-(dirres/2))) * pi/180, direction = -1) +
      scale_fill_manual(name = "treated", values = color, drop = FALSE) +
      scale_color_manual(name = "treated", values = c('black','black'), drop = FALSE) +
      theme_bw() +
      theme(axis.title.x = element_blank(), legend.position="none")
  } else {
    # assign each wind direction to a bin
    Cdir.binned <- cut(control$angle_deg, breaks = dir.breaks, ordered_result = TRUE)
    levels(Cdir.binned) <- dir.labels
    control$dir.binned <- Cdir.binned
    # generate bin sums for each samples as % of total cells
    C_data <- control %>% group_by(sample_side, dir.binned) %>% dplyr::summarise(count= n()) %>% mutate(y = count/sum(count))
    # average them
    C_data <- C_data %>% group_by(dir.binned) %>% dplyr::summarize(z = mean(y))
    labels <- data.frame(x = pi, y = scales::extended_breaks()(range(C_data$z)))
    print(max(C_data$z)*100)
    C_data <- C_data %>% mutate(treatment = 'zDMSO_combined')
    T_data <- T_data %>% mutate(treatment = data$treatment[1])
    data_new <- rbind(C_data, T_data)
    
    p.windrose <- ggplot(data = data_new, aes(x = dir.binned, y = z, group=treatment, fill = treatment, color = treatment, alpha = treatment)) +
      geom_bar(width = 1, linewidth = .5, stat='identity',position='identity') +
      scale_x_discrete(drop = FALSE, labels = waiver()) +
      scale_y_continuous(limits = c(0, 0.042), expand = c(0, 0),  breaks = c(0,.01,.02,.03,.04)) +
      coord_polar(start = ((270-(dirres/2))) * pi/180, direction = -1) +
      scale_fill_manual(name = "treatment", values = color, drop = FALSE) +
      scale_color_manual(name = "treatment", values = c('black','black'), drop = FALSE) +
      scale_alpha_discrete(range = c(1, .5)) +
      theme_bw() +
      theme(axis.title.x = element_blank(), legend.position="none")
  }
  return(p.windrose)
}


#### 0.5 Mollweide projection helpers ####
plot_flatearth <- function(data) {
  proj_type <- "+proj=moll  +lat_0=0 +lon_0=0"
  proj_pts <- project(vect(data), proj_type) # makes a spacvector, which is just the points projected onto globe
  r <- rast(nrow = 300, ncol = 600, xmin=-180, xmax = 180, ymin = -90, ymax=90) # creates a spatraster which is the squares
  r <- project(r, proj_type, res=18000000/6) # projects our spatraster squares over globe
  grat <- sf::st_graticule(lon=seq(-180,180, 60),
                      lat = seq(-90,90,45),
                      ndiscr = 5000) %>%
                      vect() %>%
                      project(proj_type)
  total <- nrow(data)
  obs <- rasterize(proj_pts, r, fun=function(i){length(i)/total*100})
  # f <- function(x)
  p <- plot(obs, axes=TRUE, clip = TRUE)
  plot(grat, add=TRUE)
  return(p)
}


generate_flatearth <- function(data, proj_type, mapres) {
  if(missing(mapres)) {
    mapres = 8
  }
  proj_pts <- project(vect(data), proj_type) # makes a spacvector, which is just the points projected onto globe
  r <- rast(nrow = 300, ncol = 600, xmin=-180, xmax = 180, ymin = -90, ymax=90) # creates a spatraster which is the squares
  # r <- rast(resolution=60, xmin=-180, xmax = 180, ymin = -90, ymax=90) # creates a spatraster which is the squares
  r <- project(r, proj_type, res=18000000/mapres) # projects our spatraster squares over globe
  total <- nrow(data)
  obs <- rasterize(proj_pts, r, fun=function(i){length(i)/total*100})
  return(obs)
}


plot_diff_flatearth <- function(baseline_data, diff_data, palette) {
  proj_type <- "+proj=moll  +lat_0=0 +lon_0=0"
  obs_bl <- generate_flatearth(baseline_data, proj_type, 16)
  obs_diff <- generate_flatearth(diff_data, proj_type, 16)
  obs_rel <- (obs_diff - obs_bl) / obs_bl * 100
  obs_rel <- focal(obs_rel, w=3, fun='mean')
  temp <- data_frame(val = values(obs_rel))
  cutoff <- quantile(temp, na.rm=TRUE)[3][[1]] + 1.5*IQR(temp$val, na.rm = TRUE)
  print(cutoff)
  if (cutoff < 100) {cutoff <- 100}
  # p1 <- ggplot(temp, aes(x=val)) + geom_histogram()
  # png('')
  # p <- plot(obs_rel, smooth=TRUE, axes=TRUE, col = pal, range=c(-100,cutoff))
  if(missing(palette)) {
    plot(obs_rel, smooth=TRUE, axes=TRUE, range=c(-100,cutoff), clip=FALSE)
  } else {
    plot(obs_rel, smooth=TRUE, axes=TRUE, range=c(-100,cutoff), clip=FALSE, col=palette)
  }

  grat <- sf::st_graticule(lon=seq(-180,180, 30),
                    lat = seq(-90,90,30),
                    ndiscr = 5000) %>%
                    vect() %>%
                    project(proj_type)
  plot(grat, add=TRUE)
}


#### 0.6 BPNME statistics helpers ####
get_posterior_estimates <- function(fitted) {
  a1 <- fitted$beta1[,1]
  a2 <- fitted$beta2[,1]
  b1 <- fitted$beta1[,2]
  b2 <- fitted$beta2[,2]

  # Compute zeta + spread control
  zeta   <- sqrt((a1)^2 + (a2)^2)^2/4
  spread <- 1 - sqrt((pi * zeta)/2) * exp(-zeta) *
    (besselI(zeta, 0) + besselI(zeta, 1)) #compute posterior spread

  # Compute zeta + spread treatment
  zeta.b   <- sqrt((a1 + b1)^2 + (a2 + b2)^2)^2/4
  spread.b <- 1 - sqrt((pi * zeta.b)/2) * exp(-zeta.b) *
    (besselI(zeta.b, 0) + besselI(zeta.b, 1))  #compute posterior spread

  # Get posterior summary estimates for table
  mode_est(spread)
  mean(spread)
  sd(spread)

  mode_est(spread.b)
  mean(spread.b)
  sd(spread.b)

  hpd_est(spread)
  hpd_est(spread.b)
  print(hpd_est(spread))
  print(hpd_est(spread.b))
  return(list(spread, spread.b))
}


#### Main ####
#### 1 Setting up the DF ####
# R doesn't have a good way to get the path of where this file is located, unfortunately
# if you are running this code in Rstudio, try this:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("../../data/GM130")
getwd() #check our working directory
# if you aren't using rstudio, use the command setwd() 
# and point it to the data/GM130 directory


### Prepare the cellpose assisted matching angles
# This step can take a long time! More than a minute
# if you've run it before then just load the combined data here
df <- readRDS("combined_angles.Rda")
# else you generate them again here and then
# save the combined dataset for faster loading next time
df <- compile_angle_df() # this will take a while! be patient
saveRDS(df, file='combined_angles.Rda') 


### Correct angles on any images that are not 'square'
# rotates all angles in a given image by a set angle
# in this data we have a horizontal angle pointing right as 0 degrees
# CCW rotation is positive (right hand rule)
# I draw a line in in imagej from the left to right along the 'flat' of the FNP, 
# then click 'measure' (ctrl+m) and get the angle
# just subtract this angle to all other angles
# see the example image in the github readme
info_df <- read.csv('GM130_image_log.csv')
df <- adjust_base_angles(df, info_df)


### Flip all the control angles across the 'y' axis so they match w/ treatment side
# NOTE: For some reason the LY group ones are reversed? Probably how the section was put on the slide...
# this is very important part of the analysis
df <- flip_y_angles(df)


### filter pairs based on imagej rois ###
# this gets rid of angles in the ectoderm or other problem areas like neurectoderm
# I manually drew polygon regions of interest in ImageJ around the mesenchyme, 
# excluding ectoderm and neurectoderm.
# see example image in the github readme
# This code loads in these .roi files and converts them to usable polygons for filtering.
# The read.ijroi unfortunately imports the 'polygon' as a series of points instead of lines,
# so we need to use the 'simple features' sf library to convert roi points into lines, 
# and then those lines into polygons. sf can then see which centroids are within the polygon
mask_filenames <- list.files(path="./imagej_rois/region_mask_rois/", pattern=".roi$", recursive=TRUE)
df_masked <- mask_out_centroids(df, mask_filenames, TRUE)


### filter based on distance from the 'baseline' roi line
# this is very similar to the last chuck of code above
# Here I made line roi's in ImageJ along the top of the nasalpit, globular process, and lateral FNP
# I'll put another example in the github readme
# We will only look at nuclei-Golgi pairs that are within 200 um of this line
# Basically we don't want to analyze pairs that are more in the mid-FNP or too 
# close to the neural ectoderm
baseline_mask_filenames <- list.files(path="./imagej_rois/baseline_rois/", pattern=".roi$")
# this is another slow one, be patient
df_baseline_masked <- filter_baseline_distance(baseline_mask_filenames, df_masked, 200, TRUE)


#### 2 Cellularity #### 
### Calculate cells per area in the masks above
# This will match the file names of the two different ImageJ roi files we used
# we will use this to calculate area
matches <- data.frame()
for (i in 1:length(baseline_mask_filenames)) {
  match <- str_which(mask_filenames, str_replace(baseline_mask_filenames[i], '_baseline.roi', '.roi'))
  matches <- rbind(matches, c(i, match))
}


# this bit of code here will calculate area and number of nuclei in that area for 'cellularity'
res <- 0.2840910 # this is specific to the images for this project, they all have same xy resolution
df_baseline_masked <- get_cellularity(df_baseline_masked, matches, mask_filenames, baseline_mask_filenames, 200, res)


### Check cellularity & plot it
df_cellularity <- df_baseline_masked %>% group_by(treatment, side, t) %>% summarise(cells_per_sqmicron)
df_cellularity <- df_cellularity %>% distinct() %>% drop_na()
df_cellularity$cells_per_sqmicron <- as.numeric(df_cellularity$cells_per_sqmicron)
df_cellularity <- df_cellularity %>% mutate(cells_per_mm2 = cells_per_sqmicron * 1000*1000)

cell_plot <- ggplot(df_cellularity, aes(fill = side, y=cells_per_sqmicron, x=treatment)) +
    geom_bar(position="dodge", stat="summary", fun=mean) + geom_errorbar(stat='summary', fun.data = mean_sdl, position = position_dodge(0.9))
plot(cell_plot)


#### 3 Watson U2 tests ####
### First do some traditional summary statistics
circular_statistics <- list()
# get actual circle stats for export
for (i in 1:length(levels(as.factor(df_baseline_masked$treatment)))) {
  for (j in 1:length(levels(as.factor(df_baseline_masked$side)))) {
    temp_circ_stats <- get_circular_stats(levels(as.factor(df_baseline_masked$treatment))[i], levels(as.factor(df_baseline_masked$side))[j], df_baseline_masked)
    circular_statistics <- append(circular_statistics, list(temp_circ_stats))
  }
}

circular_statistics <- do.call(rbind.data.frame, circular_statistics)
colnames(circular_statistics) <- c('treatment', 'side', 'ind_mean', 'overall_mean', 'ind_sd', 'mean_sd', 'overall_sd', 'ind_cr', 'mean_cr', 'overall_cr')
circular_statistics <- circular_statistics %>% mutate(across(ind_mean:overall_sd, ~.x*180/pi)) # convert to degrees
circular_statistics <- circular_statistics %>% mutate(across(ind_mean:overall_sd, ~case_when(. < 0 ~ . + 360, TRUE ~ .))) # remove negative

### compare side v side for each treatment
watson_side_results <- list()
for (i in 1:length(levels(as.factor(df_baseline_masked$treatment)))) {
  watson <- watson.two.test(df_baseline_masked %>% filter(treatment == levels(as.factor(df_baseline_masked$treatment))[i], side =='control') %>% select(angle), df_baseline_masked %>% filter(treatment == levels(as.factor(df_baseline_masked$treatment))[i], side == 'treated') %>% select(angle))
  # print(watson)
  watson_result_list <- list(levels(as.factor(df_baseline_masked$treatment))[i], watson[[1]])
  watson_side_results <- append(watson_side_results, list(watson_result_list))
}

### compare each side-treatment combo against DMSO combined
watson_results <- list()
for (i in 1:length(levels(as.factor(df_baseline_masked$treatment)))) {
  for (j in 1:length(levels(as.factor(df_baseline_masked$side)))) {
    watson <- watson.two.test(df_baseline_masked %>% filter(treatment == 'DMSO') %>% select(angle), df_baseline_masked %>% filter(treatment == levels(as.factor(df_baseline_masked$treatment))[i], side == levels(as.factor(df_baseline_masked$side))[j]) %>% select(angle))
    # print(watson)
    watson_result_list <- list(levels(as.factor(df_baseline_masked$treatment))[i], levels(as.factor(df_baseline_masked$side))[j], watson[[1]])
    watson_results <- append(watson_results, list(watson_result_list))
  }
}

# cleanup results into a nice df for export
watson_results <- do.call(rbind.data.frame, watson_results)
watson_side_results <- do.call(rbind.data.frame, watson_side_results)
colnames(watson_results) <- c('treatment', 'side', 'statistic')
colnames(watson_side_results) <- c('treatment', 'contralateral_statistic')
watson_results <- left_join(watson_results, watson_side_results)
# https://github.com/cran/circular/blob/master/R/watson.two.test.R
watson_results <- watson_results %>% mutate(pval = case_when(statistic > 0.385 ~ '< 0.001',
                                                             statistic > 0.268 ~ '< 0.01',
                                                             statistic > 0.187 ~ '< 0.05',
                                                             statistic > .152 ~ '< 0.1',
                                                             TRUE ~ 'greater than 0.1'))
watson_results <- watson_results %>% mutate(vs_contralateral_pval = case_when(contralateral_statistic > 0.385 ~ '< 0.001',
                                                                              contralateral_statistic > 0.268 ~ '< 0.01',
                                                                              contralateral_statistic > 0.187 ~ '< 0.05',
                                                                              contralateral_statistic > .152 ~ '< 0.1',
                                                                              TRUE ~ 'greater than 0.1'))
# combine with stats from above
combined_results <- full_join(watson_results, circular_statistics, by = c("treatment", "side"))
write.csv(combined_results, 'Golgi_analysis_output.csv')

#  old, maybe delete
cdat <- DMSO_control.circ[[1]]
sample_ids <- unique(cellpose_DMSO_control$id)
g <- length(sample_ids)
ndat <- vector("integer", g)
for (i in 1:g) {
  ndat[i] <- nrow(cellpose_DMSO_control[cellpose_DMSO_control$id == sample_ids[i],])
}
gID <- as.integer(as.factor(cellpose_DMSO_control$id))
WalraffTest(cdat,ndat,g,gID)


#### 4 Windrose plotting ####
graphing_df <- df_baseline_masked
graphing_df$angle_deg <- graphing_df$angle * 180 / pi
graphing_df$rel_z <- graphing_df$delta_z / graphing_df$distance
graphing_df <- graphing_df %>% mutate(side_spd = case_when(side == 'control' ~ 0,
                                                           side == 'treated' ~ 1))
graphing_df <- graphing_df %>% unite(sample_side, sample_info, side, sep='_', remove=FALSE)

graphing_df$treatment <- factor(graphing_df$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', '3Mix'))
graphing_df <- graphing_df %>% mutate(treatment = recode(treatment, '3Mix' = 'Triple'))
graphing_df$side <- factor(graphing_df$side, levels = c('control', 'treated'))
graphing_df <- graphing_df %>% mutate(side = recode(side, 'control' = 'contralateral'))
for (i in 1:length(levels(graphing_df$treatment))) {
  filter_data <- graphing_df %>% filter(as.integer(treatment) == i)
  
  # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"), units='in', width=5, height=5, res=300)
  pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_contralateral.pdf"), width=5, height=5)
  control_plot <- plot.windrose(filter_data %>% filter(side == 'contralateral'), 'white', dirres = 10)
  plot(control_plot)
  dev.off()
  # dev.off()
  
  # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"), units='in', width=5, height=5, res=300)
  pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated.pdf"), width=5, height=5)
  treated_plot <- plot.windrose(filter_data %>% filter(side == 'treated'), 'red', dirres = 10)
  plot(treated_plot)
  dev.off()
  # dev.off()
}

for (i in 1:length(levels(graphing_df$treatment))) {
  filter_data <- graphing_df %>% filter(as.integer(treatment) == i)
  
  # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"), units='in', width=5, height=5, res=300)
  pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_contralateral_vsDMSO.pdf"), width=5, height=5)
  control_plot <- plot.windrose(filter_data %>% filter(side == 'contralateral'), c('white', 'grey'), dirres = 10, graphing_df %>% filter(treatment == 'DMSO'))
  plot(control_plot)
  dev.off()
  # dev.off()
  
  # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"), units='in', width=5, height=5, res=300)
  pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.pdf"), width=5, height=5)
  treated_plot <- plot.windrose(filter_data %>% filter(side == 'treated'), c('red', 'grey'), dirres = 10, graphing_df %>% filter(treatment == 'DMSO'))
  plot(treated_plot)
  dev.off()
  # dev.off()
}


#### 5 Mollweide Plots ####
df_3d_proj <- df_baseline_masked
df_3d_proj$lat <- (pi/2 - acos(df_3d_proj$unit_z))*180/pi
df_3d_proj$lon <- df_3d_proj$angle*180/pi + 270
df_3d_proj <- df_3d_proj %>% mutate(lon = case_when(lon > 360 ~ lon - 360,
                                                    TRUE ~ lon))
df_3d_proj_sf <- st_as_sf(df_3d_proj, coords = c("lon", "lat"), crs = '+proj=longlat +datum=WGS84 +no_defs +type=crs')


df_3d_proj_sf$treatment <- factor(df_3d_proj_sf$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', '3Mix'))
df_3d_proj_sf$side <- factor(df_3d_proj_sf$side, levels = c('control', 'treated'))
proj_type <- "+proj=moll  +lat_0=0 +lon_0=0"
palette<-colorRampPalette(brewer.pal(9,'Purples'))(100)
grat <- sf::st_graticule(lon= seq(-180,180, 60), #c(-180,-120,0,120,180),
                         lat = seq(-90,90,45),
                         ndiscr = 5000) %>%
                         vect() %>%
                         project(proj_type)
for (i in 1:length(levels(df_3d_proj_sf$treatment))) {
  filter_data <- df_3d_proj_sf %>% filter(as.integer(treatment) == i)
  fe_plot_ctrl <- generate_flatearth(filter_data %>% filter(side == 'control'), proj_type, mapres=16)
  # fe_plot_ctrl <- focal(fe_plot_ctrl, w=3, fun="mean")
  fe_plot_trt <- generate_flatearth(filter_data %>% filter(side == 'treated'), proj_type, mapres=16)
  temp_ctrl <- data_frame(val = values(fe_plot_ctrl))
  print(paste0(as.character(filter_data$treatment[1]), ' control ', max(temp_ctrl, na.rm=TRUE)))
  temp_trt <- data_frame(val = values(fe_plot_trt))
  print(paste0(as.character(filter_data$treatment[1]), ' treated ', max(temp_trt, na.rm=TRUE)))

  png(filename=paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_control.png"), units="in", width=5, height=3, res=300)
  plot(fe_plot_ctrl, smooth=TRUE, axes=FALSE, legend=FALSE, clip=FALSE, col=palette, range = c(0,1))
  plot(grat, add=TRUE)
  dev.off()

  png(filename=paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_treated.png"), units="in", width=5, height=3, res=300)
  plot(fe_plot_trt, smooth=TRUE, axes=FALSE, clip=FALSE, legend=FALSE, col=palette, range = c(0,1))
  plot(grat, add=TRUE)
  dev.off()
}


palette<-colorRampPalette(brewer.pal(11,'PuOr'))(100)
for (i in 1:length(levels(df_3d_proj_sf$treatment))) {
  filter_data <- df_3d_proj_sf %>% filter(as.integer(treatment) == i)

  png(filename=paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_control_vsDMSO2.png"), units="in", width=5, height=5, res=300)
  plot_diff_flatearth(df_3d_proj_sf %>% filter(treatment == 'DMSO', side == 'control'), filter_data %>% filter(side == 'control'), palette)
  # plot(grat, add=TRUE)
  dev.off()

  png(filename=paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_treated_vsDMSO2.png"), units="in", width=5, height=5, res=300)
  plot_diff_flatearth(df_3d_proj_sf %>% filter(treatment == 'DMSO', side == 'control'), filter_data %>% filter(side == 'treated'), palette)
  # plot(grat, add=TRUE)
  dev.off()
}


# combine flatearth plots
i <- 1
filter_data <- df_3d_proj_sf %>% filter(as.integer(treatment) == i)
fe_plot_ctrl <- image_read(paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_control.png"))
fe_plot_trt <- image_read(paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_treated.png"))

image_col <- c(fe_plot_ctrl,fe_plot_trt)
imgs <- image_append(image_scale(image_col, "x200"), stack = TRUE)

imgs <- image_border(imgs, "white", "85x15")
imgs <- image_annotate(imgs, as.character(filter_data$treatment[1]), font = "Times",
                       location = "+194+0", size = 40)
imgs <- image_annotate(imgs, paste0("contralateral"), font = "Times",
                       location = "+20+200", size = 32, degrees = -90)
imgs <- image_annotate(imgs, paste0("treated"), font = "Times",
                       location = "+20+370", size = 32, degrees = -90)
image_browse(imgs)


for (i in 2:length(levels(df_3d_proj_sf$treatment))) {
  filter_data <- df_3d_proj_sf %>% filter(as.integer(treatment) == i)
  fe_plot_ctrl <- image_read(paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_control.png"))
  fe_plot_trt <- image_read(paste0("./figs/flatearth_", as.character(filter_data$treatment[1]), "_treated.png"))


  image_col <- c(fe_plot_ctrl,fe_plot_trt)
  img2 <-image_append(image_scale(image_col, "x200"), stack = TRUE)
  img2 <-image_border(img2, "white", "0x15")
  img2 <-image_annotate(img2, as.character(filter_data$treatment[1]), font = "Times",
                       location = "+100+0", size = 40)

  imgs <-c(imgs,img2)

  imgs <-image_append(image_scale(imgs))
}
image_browse(imgs)
image_write(imgs, path = paste0("./figs/flatearth_combined.png"), format = "png")


#### 6 BPNME statistics ###
# This is a much better way of doing the stats for these data, but unfortunately
# This is all old code that does not work with my data
# There are too many angle measurements per sample, it just can't run...
df_LY <- df %>% filter(treatment == 'LY')
df_U0 <- df %>% filter(treatment == 'U0')
df_U73 <- df %>% filter(treatment == 'U73')

# fit.LY <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_LY, its = 100)
fit.LY_corner <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_LY %>% filter((side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx > (2048-800)) | (side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx < 800)), its = 1000, burn = 100, n.lag = 3, seed = 101)
fit.U0_corner <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_U0 %>% filter((side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx > (2048-800)) | (side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx < 800)), its = 1000, burn = 100, n.lag = 3, seed = 101)
fit.U73_corner <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_U73 %>% filter((side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx > (2048-800)) | (side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx < 800)), its = 1000, burn = 100, n.lag = 3, seed = 101)

# Extract posterior estimates linear coefficients

fit.LY_corner$spreads <- get_posterior_estimates(fit.LY_corner) # variance in treated is (nearly) larger than control!
saveRDS(fit.LY_corner, 'LY_corner_fit.rds')
fit.LY_corner <- readRDS('LY_corner_fit.rds')

fit.U0_corner$spreads <- get_posterior_estimates(fit.U0_corner)
saveRDS(fit.U0_corner, 'U0_corner_fit.rds')
fit.U0_corner <- readRDS('U0_corner_fit.rds')

fit.U73_corner$spreads <- get_posterior_estimates(fit.U73_corner)
saveRDS(fit.U73_corner, 'U73_corner_fit.rds')
fit.U73_corner <- readRDS('U73_corner_fit.rds')

### END BPNME Stats ###