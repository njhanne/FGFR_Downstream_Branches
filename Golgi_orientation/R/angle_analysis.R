### Load libraries for packages
# circle stats stuff
library(circular)
# library(bpnreg)

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
# if this one takes more than 60s to download it will fail
# run this line: options(timeout = max(1000, getOption("timeout")))

# graphing
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(magick)


### Directory
setwd("C:\\Users\\njhan\\Box\\FGF_inhibitor_paper_5-26-2020\\data\\GM130") # home
setwd("C:\\Users\\nhanne\\Box\\FGF_inhibitor_paper_5-26-2020\\data\\GM130") # work
setwd("/Users/nhanne/Library/CloudStorage/Box-Box/FGF_inhibitor_paper_5-26-2020/data/GM130") # laptop
getwd() #check our working directory


### Prepare the cellpose assisted matching angles
# if you've run it before then just load the combined data here
df <- readRDS("combined_angles.Rda")

# else you generate them again here
  filenames <- dir(".", "*._angle_results.csv") # get all output csv from python
  types <- sub('(\\A*)_angle_results.*', '\\1', filenames) # get the sample names
  # combine all the individual csv into one big dataframe
  # this can be quite slow as there are tons of data
  df<- adply(data.frame(f=I(filenames), t=types), 1,
                with, cbind(read.csv(f), sample_info=t))
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

  # save the combined dataset for faster loading next time
  saveRDS(df, file='combined_angles.Rda')


### Correct angles on any images that are not 'square'
# rotates all angles in a given image by a set angle
# in this data we have a horizontal angle pointing towards the right side of the image is 0 degrees
# CCW rotation is positive (Right hand rule)
# need to draw a line in in imagej from the left to right along the 'flat' of the FNP, then click 'measure' (ctrl+m) and get the angle
# just add this angle to all other angles (actually subtract it...)
# the angle is in the 'GM130_image_log.csv' but was added after most of the python analysis so isn't in the df
# re-running the python analysis should get it included but that will take forever so this here's a workaround
info_df <- read.csv('GM130_image_log.csv')
angle_adjustment_df <- info_df %>% filter(!is.na(angle_adjustment)) %>% select(new_filename, angle_adjustment) %>% distinct()
angle_adjustment_df <- angle_adjustment_df %>% mutate(angle_adjustment = (angle_adjustment * pi)/180) # convert to radians
df <- left_join(df, angle_adjustment_df, by=c('t' = 'new_filename')) # put the adjustment into our main df
df$angle_adjustment[is.na(df$angle_adjustment)] <- 0 # change na's to zero
df$angle_old <- df$angle # just in case we save the og angles
df <- df %>% mutate(angle = angle - angle_adjustment) # 'add' the new angle adjustment
df <- df %>% mutate(angle = case_when(angle > 2*pi ~ angle - 2*pi,
                                      angle < 0 ~ angle + 2*pi,
                                      TRUE ~ angle))

### filter pairs based on imagej rois ###
# this helps get rid of angles in the ectoderm or other problem areas
# need to use the 'simple features' sf library to turn imagej roi points into lines, and then those lines into polygons
# sf can then see which centroids are within the polygon

# if you've run it before then just load the combined data here
df_masked <- readRDS("combined_angles_masked.Rda")

# if not then rerun all this
  mask_out_centroids <- function(filtered_df, imagej_mask, type_column) {
    centroids <- filtered_df %>% select(nuclei_centroidx, nuclei_centroidy)
    mask_linestring <- st_linestring(imagej_mask$coords, dim='XY')
    mask_polygon <- st_cast(mask_linestring, 'POLYGON')
    keep <- st_intersection(st_multipoint(data.matrix(centroids), dim='XY'), mask_polygon)

    p <- ggplot() + geom_sf(data = mask_polygon) +
      geom_point(data = centroids, aes(nuclei_centroidx, nuclei_centroidy), color = 'red') + geom_sf(data=keep)
    file_name <- paste(type_column, "_masked_nuclei.svg")
    ggsave(filename=file_name, p)

    keep <- as.data.frame(as.matrix(keep))
    filtered_df <- inner_join(filtered_df, keep, by=c('nuclei_centroidx' = 'V1', 'nuclei_centroidy' = 'V2'))
    return(filtered_df)
  }

  mask_filenames <- list.files(path="./new_images", pattern=".roi$", recursive=TRUE)
  df_masked <- data.frame()
  for (i in 1:length(mask_filenames)) {
    # reading in the imagej roi creates a matrix of xy coordinates which can be plugged into sf
    mask <- read.ijroi(file.path("./new_images", mask_filenames[i]), verbose = FALSE)
    type_col <- sub('(\\A*).roi$', '\\1', basename(mask_filenames[i])) # get column name for matching
    filtered_df <- mask_out_centroids(df %>% filter(t == type_col), mask, type_col)
    df_masked <- rbind(df_masked, filtered_df) # probably not supposed to do this in a loop but it works
  }
  df_masked <- rbind(df_masked, df %>% filter(df$t %in% setdiff(unique(df$t), unique(df_masked$t))))
  # save the combined dataset for faster loading next time
  saveRDS(df_masked, file='combined_angles_masked.Rda')

### filter based on distance from the 'baseline' roi line
filter_distance_centroids <- function(filtered_df, baseline_roi, type_column, thresh_dist) {
  centroids <- filtered_df %>% select(nuclei_centroidx, nuclei_centroidy)
  baseline_linestring <- st_linestring(baseline_roi$coords, dim='XY')
  keep <- st_intersection(st_multipoint(data.matrix(centroids), dim='XY'), st_buffer(baseline_linestring, thresh_dist/0.2840910))

  # p <- ggplot() + geom_sf(data = baseline_linestring) +
  #   geom_point(data = centroids, aes(nuclei_centroidx, nuclei_centroidy), color = 'red') + geom_sf(data=keep)
  # file_name <- paste(type_column, "_masked_baseline_nuclei.svg")
  # ggsave(filename=file_name, p)

  keep <- as.data.frame(as.matrix(keep))
  filtered_df <- inner_join(filtered_df, keep, by=c('nuclei_centroidx' = 'V1', 'nuclei_centroidy' = 'V2'))
  return(filtered_df)
}

baseline_mask_filenames <- list.files(path="./sf_rois", pattern=".roi$")
df_baseline_masked <- data.frame()
for (i in 1:length(baseline_mask_filenames)) {
  # reading in the imagej roi creates a matrix of xy coordinates which can be plugged into sf
  baseline_roi <- read.ijroi(file.path("./sf_rois", baseline_mask_filenames[i]), verbose = FALSE)
  type_col <- sub('(\\A*)_baseline.roi$', '\\1', basename(baseline_mask_filenames[i])) # get column name for matching
  filtered_df <- filter_distance_centroids(df_masked %>% filter(t == type_col), baseline_roi, type_col, 200)
  df_baseline_masked <- rbind(df_baseline_masked, filtered_df) # probably not supposed to do this in a loop but it works
}
df_baseline_masked <- rbind(df_baseline_masked, df_masked %>% filter(df_masked$t %in% setdiff(unique(df_masked$t), unique(df_baseline_masked$t))))


### little aside to calculate cells per area in the masks above
mask_inter_area <- function(area_area, baseline_roi, thresh_dist) {
  res <- 0.2840910
  baseline_linestring <- st_linestring(baseline_roi$coords, dim='XY')
  baseline_buffer <- st_buffer(baseline_linestring, thresh_dist/res)

  mask_linestring <- st_linestring(area_mask$coords, dim='XY')
  mask_polygon <- st_cast(mask_linestring, 'POLYGON')
  if (st_is_valid(mask_polygon))  {
    area <- st_area(st_intersection(mask_polygon, baseline_buffer)) * (res^2)
  } else {
    area <- st_area(st_intersection(st_buffer(mask_polygon, 0), baseline_buffer)) * (res^2)
  }

  # p <- ggplot() + geom_sf(data = mask_polygon) +
  #   geom_sf(data = baseline_buffer, color='red') +
  #   geom_sf(data = st_intersection(mask_polygon, baseline_buffer), color='purple') #+
  # plot(p)

  return(area)
}

mask_filenames <- list.files(path="./new_images", pattern=".roi$", recursive=TRUE)
baseline_mask_filenames <- list.files(path="./sf_rois", pattern=".roi$")

# match the two different imagej roi things
matches <- data.frame()
for (i in 1:length(baseline_mask_filenames)) {
  match <- str_which(mask_filenames, str_replace(baseline_mask_filenames[i], '_baseline.roi', '.roi'))
  matches <- rbind(matches, c(i, match))
}
cellularity <- data.frame()
for (i in 1:nrow(matches)) {
  # reading in the imagej roi creates a matrix of xy coordinates which can be plugged into sf
  j <- matches[i,2]
  area_mask <- read.ijroi(file.path("./new_images", mask_filenames[j]), verbose = FALSE)
  type_col <- sub('(\\A*).roi$', '\\1', basename(mask_filenames[j])) # get column name for matching
  # cellularity <- mask_cellularity(df %>% filter(t == type_col), mask_area, type_col)

  baseline_roi <- read.ijroi(file.path("./sf_rois", baseline_mask_filenames[i]), verbose = FALSE)
  # type_col <- sub('(\\A*)_baseline.roi$', '\\1', basename(baseline_mask_filenames[i])) # get column name for matching
  area <- mask_inter_area(area_mask, baseline_roi, 200)
  cells <- nrow(df_baseline_masked %>% filter(t == type_col))
  cellularity <- rbind(cellularity, c(type_col, area, cells, cells/area)) # probably not supposed to do this in a loop but it works
}

colnames(cellularity) <- c('t', 'area', 'cells', 'cells_per_sqmicron')
df_baseline_masked <- full_join(df_baseline_masked, cellularity)

# flip all the control angles across the 'y' axis so they match w/ treatment side
# don't run this directly, it gets called later
### NOTE: For some reason the LY group ones are reversed? Probably how the section was put on the slide...
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

df <- flip_y_angles(df)
df_masked <- flip_y_angles(df_masked)
df_baseline_masked <- flip_y_angles(df_baseline_masked)


# Check cellularity
df_cellularity <- df_baseline_masked %>% group_by(treatment, side, t) %>% summarise(cells_per_sqmicron)
df_cellularity <- df_cellularity %>% distinct() %>% drop_na()
df_cellularity$cells_per_sqmicron <- as.numeric(df_cellularity$cells_per_sqmicron)
df_cellularity <- df_cellularity %>% mutate(cells_per_mm2 = cells_per_sqmicron * 1000*1000)

cell_plot <- ggplot(df_cellularity, aes(fill = side, y=cells_per_sqmicron, x=treatment)) +
    geom_bar(position="dodge", stat="summary", fun.y=mean) + geom_errorbar(stat='summary', fun.data = mean_sdl, position = position_dodge(0.9))
plot(cell_plot)

# df <- df %>% mutate(plane_angles = case_when(angle < pi/2 ~ angle + pi,
#                                              angle > 3*pi/2 ~ angle + pi,
#                                              TRUE ~ angle))

#OLD
# df_handangle  <- flip_y_angles(df_handangle)
# matching_samples <- intersect(df$sample_info, df_handangle$sample_info)

## OLD
# 0,0 is top left corner of image, all images are 2048x2048
# can choose y > 1024 to get bottom half?
# cellpose_wide <- df %>% pivot_wider(names_from = id, values_from = angle)
# cellpose_wide <- cellpose_wide %>% select(last_col(0:6)) # change based on number of samples
# cellpose_wide <- apply(cellpose_wide,2, na.omit)

# hand <- df_handangle %>% filter(sample_info %in% matching_samples[27:36] & side == 'control')
# cellpose <- df %>% filter(sample_info %in% matching_samples[1:14] & side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx > (2048 - 1200))
# cellpose2 <- df %>% filter(sample_info %in% matching_samples[1:14] & side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx < 1200)#& nuclei_centroidy > 800 & nuclei_centroidx < 1200)
# cellpose <- df %>% filter(sample_info == matching_samples[5] & nuclei_centroidy > 800 & nuclei_centroidx < 1200)
# cellpose2 <- df %>% filter(sample_info == matching_samples[5])

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


grouped.rose.diag <- function(x, pch = 16, cex = 1, axes = TRUE, shrink = 1, bins = NULL, upper = TRUE, ticks = TRUE,
                             tcl = 0.025, tcl.text = 0.125, radii.scale = c("sqrt", "linear"), border = NULL, col = NULL,
                             tol = 0.04, uin = NULL, xlim = c(-1, 1), ylim = c(-1, 1), prop = 1, digits = 2, plot.info = NULL,
                             units = NULL, template = NULL, zero = NULL, rotation = NULL, main = NULL, sub = NULL,
                             xlab = "", ylab = "", add = FALSE, control.circle = circle.control(), grps = NULL, ...) {
  # customize the circular rose.diag to work with multiple samples of unequal sampling
  radii.scale <- match.arg(radii.scale)
  if (is.matrix(x) | is.data.frame(x)) {
    nseries <- ncol(x) # edited
  } else {
    nseries <- 1
  }
  if (!is.null(grps)) {
    nsamples <- length(unique(grps))
  } else {
    nsamples <- 0
  }
  xx <- as.data.frame(x)
  xcircularp <- attr(as.circular(xx[, 1]), "circularp")
  modulo <- xcircularp$modulo
  if (is.null(units)) units <- xcircularp$units
  if (is.null(plot.info)) {
    if (is.null(template)) template <- xcircularp$template
    if (template == "geographics" | template == "clock24") {
      zero <- pi / 2
      rotation <- "clock"
    } else if (template == "clock12") {
      zero <- pi / 2
      rotation <- "clock"
      modulo <- "pi"
    } else {
      if (is.null(zero)) zero <- xcircularp$zero
      if (is.null(rotation)) rotation <- xcircularp$rotation
    }
    next.points <- 0
  } else {
    zero <- plot.info$zero
    rotation <- plot.info$rotation
    next.points <- plot.info$next.points
  }
  if (!add) {
    circular:::CirclePlotRad(xlim = xlim, ylim = ylim, uin = uin, shrink = shrink, tol = tol, main = main, sub = sub, xlab = xlab, ylab = ylab, control.circle = control.circle)
  }
  if (is.null(bins)) {
    bins <- NROW(x)
  } else {
    bins <- round(bins)
    if (bins <= 0) stop("bins must be non negative")
  }
  if (is.null(border)) {
    border <- seq(nseries)
  } else {
    if (length(border) != nseries) {
      border <- rep(border, nseries)[1 : nseries]
    }
  }
  pch <- rep(pch, nseries, length.out = nseries)
  if (axes) {
    axis.circular(units = units, template = template, zero = zero, rotation = rotation, digits = digits, cex = cex, tcl = tcl, tcl.text = tcl.text)
  }
  if (!is.logical(ticks)) stop("ticks must be logical")
  if (ticks) {
    at <- circular((0 : bins) / bins * 2 * pi, zero = zero, rotation = rotation)
    ticks.circular(at, tcl = tcl)
  }
  for (iseries in 1 : nseries) {
    x <- xx[, iseries]
    x <- na.omit(x)
    n <- length(x)
    if (n) {
      x <- conversion.circular(x, units = "radians", modulo = modulo)
      attr(x, "circularp") <- attr(x, "class") <- NULL
      if (template == "clock12") x <- 2 * x
      x <- x %% (2 * pi)
      if (nsamples != 0) {
        groupedRosediagRad(x, zero = zero, rotation, bins, upper, radii.scale, prop, border[iseries], col, grps, ...)
      } else {
        RosediagRad(x, zero = zero, rotation, bins, upper, radii.scale, prop, border[iseries], col, ...)
      }
    }
  }
  return(invisible(list(zero = zero, rotation = rotation, next.points = 0)))
}


groupedRosediagRad <- function (x, zero, rotation, bins, upper, radii.scale, prop, border, col, grps, ...) {
  freq <- as.data.frame(rep(0, bins))
  rel.freq <- as.data.frame(rep(0, bins))
  arc <- (2 * pi)/bins
  if (!is.logical(upper))
    stop("upper must be logical")
  breaks <- seq(0, 2 * pi, length.out = (bins + 1))
  if (upper == TRUE)
    x[x == 0] <- 2 * pi
  x[x >= 2 * pi] <- 2 * pi - 4 * .Machine$double.eps

  samples <- length(unique(grps))
  for (i in 1:samples) {
    samplen <- unique(grps)[i]
    mask <- x
    # mask <- matrix(NA,2,2)
    mask <- grps == samplen
    sample <- x[mask]
    n <- length(sample)

    freq[,i] <- hist.default(sample, breaks = breaks, plot = FALSE, right = upper)$counts
    rel.freq[,i] <- freq[,i]/n
    if (rotation == "clock")
      rel.freq[,i] <- rev(rel.freq[,i])
  }

  if (radii.scale == "sqrt") {
    radius <- sqrt(rel.freq) * prop / samples # not sure if this is correct, but I don't want to use this style anyway
  }
  else {
    grp.rel.freq <- rel.freq * prop / samples
    radius <- apply(grp.rel.freq, 1, sum)
  }

  sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
  mids <- seq(arc/2, 2 * pi - pi/bins, length = bins)
  for (i in 1:bins) {
    if (radius[i] != 0) {
      xx <- c(0, radius[i] * cos(seq(sector[i], sector[i] +
        (2 * pi)/bins, length = 1000/bins) + zero), 0)
      yy <- c(0, radius[i] * sin(seq(sector[i], sector[i] +
        (2 * pi)/bins, length = 1000/bins) + zero), 0)
      polygon(xx, yy, border = border, col = col, ...)
    }
  }
}


cellpose_windrose <- function(df_temp) {
  cellpose_control <- df_temp %>% filter(side == 'control')
  cellpose_treated <- df_temp %>% filter(side == 'treated')
  control.circ <- circularize_group(cellpose_control$angle)
  treated.circ <- circularize_group(cellpose_treated$angle)

  control_means <- group_sample_circular_stats(control.circ[[1]], cellpose_control$id, 'mean')
  treated_means <- group_sample_circular_stats(treated.circ[[1]], cellpose_treated$id, 'mean')
  control_sds <- group_sample_circular_stats(control.circ[[1]], cellpose_control$id, 'sd')
  treated_sds <- group_sample_circular_stats(treated.circ[[1]], cellpose_treated$id, 'sd')
  control_rhos <- group_sample_circular_stats(control.circ[[1]], cellpose_control$id, 'rho')
  treated_rhos <- group_sample_circular_stats(treated.circ[[1]], cellpose_treated$id, 'rho')

  control_mean <- mean.circular(unlist(control_means))
  treated_mean <- mean.circular(unlist(treated_means))
  control_sd <- mean.circular(unlist(control_sds))
  treated_sd <- mean.circular(unlist(treated_sds))
  control_rho <- mean.circular(unlist(control_rhos))
  treated_rho <- mean.circular(unlist(treated_rhos))

  grouped.rose.diag(control.circ[1], prop=15,bins=36, col="white", radii.scale = 'linear', grps = cellpose_control$id, xlab=paste("mean:", control_mean*180/pi, "\nsd:", control_sd*180/pi, "\nsd_mean:", sd.circular(unlist(control_means))*180/pi, "\ncr:", control_rho, "\ncr_mean:", rho.circular(unlist(control_means))), ylab = paste(unique(df_temp$treatment)[1], 'control'))
    arrows.circular(control_mean, y = rho.circular(unlist(control_means)))
    lines(density.circular(unlist(control.circ[1]), bw=10,  n=360, kernel = "vonmises"), lwd=2, lty=1)
  grouped.rose.diag(treated.circ[1], prop=15,bins=36, col="red", radii.scale = 'linear', grps = cellpose_treated$id, xlab=paste("mean:", treated_mean*180/pi, "\nsd:", treated_sd*180/pi, "\nsd_mean:", sd.circular(unlist(treated_means))*180/pi, "\ncr:", treated_rho, "\ncr_mean:", rho.circular(unlist(treated_means))), ylab = paste(unique(df_temp$treatment)[1], 'treated'))
    arrows.circular(treated_mean, y = rho.circular(unlist(treated_means)))
    lines(density.circular(unlist(treated.circ[1]), bw=40), lwd=2, lty=1)
  # return(list(cellpose_control, cellpose_treated, control.circ, treated.circ, control_means, treated_means, control_sds, treated_sds, control_rhos, treated_rhos))
}


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


#filter out 2d results
df_baseline_masked <- df_baseline_masked %>% filter(!is.na(unit_x))

# generate windrose plots w/ some stats on them
plot_df <- df_baseline_masked
for (i in 1:length(levels(as.factor(plot_df$treatment)))) {
  cellpose_windrose(plot_df %>% filter(treatment == levels(as.factor(plot_df$treatment))[i]))
}

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


# Watson U2 tests
# compare side v side for each treatment
watson_side_results <- list()
for (i in 1:length(levels(as.factor(df_baseline_masked$treatment)))) {
  watson <- watson.two.test(df_baseline_masked %>% filter(treatment == levels(as.factor(df_baseline_masked$treatment))[i], side =='control') %>% select(angle), df_baseline_masked %>% filter(treatment == levels(as.factor(df_baseline_masked$treatment))[i], side == 'treated') %>% select(angle))
  # print(watson)
  watson_result_list <- list(levels(as.factor(df_baseline_masked$treatment))[i], watson[[1]])
  watson_side_results <- append(watson_side_results, list(watson_result_list))
}

# compare each side-treatment combo against DMSO combined
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
colnames(watson_side_results) <- c('treatment', 'group_statistic')
watson_results <- left_join(watson_results, watson_side_results)
# https://github.com/cran/circular/blob/master/R/watson.two.test.R
watson_results <- watson_results %>% mutate(pval = case_when(statistic > 0.385 ~ '< 0.001',
                                                             statistic > 0.268 ~ '< 0.01',
                                                             statistic > 0.187 ~ '< 0.05',
                                                             statistic > .152 ~ '< 0.1',
                                                             TRUE ~ 'greater than 0.1'))
watson_results <- watson_results %>% mutate(group_pval = case_when(group_statistic > 0.385 ~ '< 0.001',
                                                             group_statistic > 0.268 ~ '< 0.01',
                                                             group_statistic > 0.187 ~ '< 0.05',
                                                             group_statistic > .152 ~ '< 0.1',
                                                             TRUE ~ 'greater than 0.1'))
# combine with stats from above
combined_results <- full_join(watson_results, circular_statistics, by = c("treatment", "side"))
write.csv(combined_results, 'Golgi_analysis_output.csv')

#OLD
# stats for bimodal distribution
# df_bimodal <- df_masked %>% mutate(old_angle2 = angle)
# df_bimodal <- df_masked %>% mutate(angle = 2*angle) %>% mutate(angle = case_when(angle > 2*pi ~ angle - 2*pi,
#                                                                                          TRUE ~ angle))
# cellpose_bimodal_DMSO_list <- cellpose_windrose(df_bimodal %>% filter(treatment == 'DMSO'))
# cellpose_bimodal_U0126_list <- cellpose_windrose(df_bimodal %>% filter(treatment == 'U0126'))
# cellpose_bimodal_U73122_list <- cellpose_windrose(df_bimodal %>% filter(treatment == 'U73122'))
# cellpose_bimodal_LY294002_list <- cellpose_windrose(df_bimodal %>% filter(treatment == 'LY294002'))
# cellpose_bimodal_mix_list <- cellpose_windrose(df_bimodal %>% filter(treatment == '3Mix'))
# cellpose_bimodal_double_list <- cellpose_windrose(df_bimodal %>% filter(treatment == 'double')) # ***does not have bimodal distribution***



# WindRose.R
# https://stackoverflow.com/a/17266781
# https://stackoverflow.com/questions/62343603/r-windrose-percent-label-on-figure
# require(ggplot2)
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
  T_data <- data %>% group_by(sample_side, dir.binned) %>% summarise(count= n()) %>% mutate(y = count/sum(count))
  T_data <- T_data %>% group_by(dir.binned) %>% summarize(z = mean(y))
  labels <- data.frame(x = pi, y = scales::extended_breaks()(range(T_data$z)))
  print(as.character(data$treatment[1]))
  print(max(T_data$z)*100)

  if(missing(control)) {
    p.windrose <- ggplot(data = T_data, aes(x = dir.binned, y = z, fill = color, color = color)) +
      geom_bar(width = 1, size = .5, stat='identity') +
      scale_x_discrete(drop = FALSE, labels = waiver()) +
      coord_polar(start = ((270-(dirres/2))) * pi/180, direction = -1) +
      scale_fill_manual(name = "treated", values = color, drop = FALSE) +
      scale_color_manual(name = "treated", values = c('black','black'), drop = FALSE) +
      theme(axis.title.x = element_blank())
  } else {
      # assign each wind direction to a bin
    Cdir.binned <- cut(control$angle_deg, breaks = dir.breaks, ordered_result = TRUE)
    levels(Cdir.binned) <- dir.labels
    control$dir.binned <- Cdir.binned
    C_data <- control %>% group_by(sample_side, dir.binned) %>% summarise(count= n()) %>% mutate(y = count/sum(count))
    C_data <- C_data %>% group_by(dir.binned) %>% summarize(z = mean(y))
    labels <- data.frame(x = pi, y = scales::extended_breaks()(range(C_data$z)))
    print(max(C_data$z)*100)
    C_data <- C_data %>% mutate(treatment = 'DMSO')
    T_data <- T_data %>% mutate(treatment = data$treatment[1])
    data_new <- rbind(C_data, T_data)

    p.windrose <- ggplot(data = data_new, aes(x = dir.binned, y = z, group=treatment, fill = treatment, color = treatment, alpha = treatment)) +
      geom_bar(width = 1, size = .5, stat='identity',position='identity') +
      scale_x_discrete(drop = FALSE, labels = waiver()) +
      scale_y_continuous(limits = c(0, 0.0431), expand = c(0, 0),  breaks = c(0,.01,.02,.03,.04)) +
      coord_polar(start = ((270-(dirres/2))) * pi/180, direction = -1) +
      scale_fill_manual(name = "treatment", values = color, drop = FALSE) +
      scale_color_manual(name = "treatment", values = c('black','black'), drop = FALSE) +
      scale_alpha_discrete(range = c(1, 0.5)) +
      theme_bw() +
      theme(axis.title.x = element_blank(), legend.position="none")
  }
return(p.windrose)
}

graphing_df <- df_baseline_masked %>% filter(!is.na(unit_x))
graphing_df$angle_deg <- graphing_df$angle * 180 / pi
graphing_df$rel_z <- graphing_df$delta_z / graphing_df$distance
graphing_df <- graphing_df %>% mutate(side_spd = case_when(side == 'control' ~ 0,
                                                           side == 'treated' ~ 1))
# graphing_df <- graphing_df %>% unite(facet_col, side, treatment, sep='_', remove=FALSE)
graphing_df <- graphing_df %>% unite(sample_side, sample_info, side, sep='_', remove=FALSE)
data <- graphing_df %>% filter(treatment == 'DMSO', side == 'control')
# graphing_df_test <- test_DMSO %>% select(rel_z, angle_deg) %>% dplyr::rename(spd = rel_z, dir = angle_deg)

graphing_df$treatment <- factor(graphing_df$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', '3Mix'))
graphing_df <- graphing_df %>% mutate(treatment = recode(treatment, '3Mix' = 'Triple'))
graphing_df$side <- factor(graphing_df$side, levels = c('control', 'treated'))
for (i in 1:length(levels(graphing_df$treatment))) {
  filter_data <- graphing_df %>% filter(as.integer(treatment) == i)

  png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"), units='in', width=5, height=5, res=300)
  pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.pdf"), width=5, height=5)
  control_plot <- plot.windrose(filter_data %>% filter(side == 'control'), c('white', 'red'), dirres = 10, graphing_df %>% filter(treatment == 'DMSO'))
  plot(control_plot)
  dev.off()
  dev.off()

  png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"), units='in', width=5, height=5, res=300)
  pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.pdf"), width=5, height=5)
  treated_plot <- plot.windrose(filter_data %>% filter(side == 'treated'), c('white', 'red'), dirres = 10, graphing_df %>% filter(treatment == 'DMSO'))
  plot(treated_plot)
  dev.off()
  dev.off()
}


# combine windrose plots
i <- 1
filter_data <- graphing_df %>% filter(as.integer(treatment) == i)
wr_plot_ctrl <- image_read(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"))
wr_plot_trt <- image_read(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"))

image_col <- c(wr_plot_ctrl,wr_plot_trt)
imgs <- image_append(image_scale(image_col, "x200"), stack = TRUE)

imgs <- image_border(imgs, "white", "85x35")
imgs <- image_annotate(imgs, as.character(filter_data$treatment[1]), font = "Times",
                       location = "+130+0", size = 40)
imgs <- image_annotate(imgs, paste0("contralateral"), font = "Times",
                       location = "+50+200", size = 32, degrees = -90)
imgs <- image_annotate(imgs, paste0("treated"), font = "Times",
                       location = "+50+370", size = 32, degrees = -90)


image_browse(imgs)

for (i in 2:length(levels(graphing_df$treatment))) {
  filter_data <- graphing_df %>% filter(as.integer(treatment) == i)
  wr_plot_ctrl <- image_read(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"))
  wr_plot_trt <- image_read(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"))


  image_col <- c(wr_plot_ctrl,wr_plot_trt)
  img2 <-image_append(image_scale(image_col, "x200"), stack = TRUE)
  img2 <-image_border(img2, "white", "0x15")
  img2 <-image_annotate(img2, as.character(filter_data$treatment[1]), font = "Times",
                       location = "+50+0", size = 40)

  imgs <-c(imgs,img2)

  imgs <-image_append(image_scale(imgs))
}
image_browse(imgs)
pdf(paste0("./figs/windrose_combined.pdf"))
plot(imgs)
dev.off()


# ptest <- plot.windrose(graphing_df %>% filter(treatment == 'U73122', side == 'treated'), c('white', 'red'),  dirres = 10, graphing_df %>% filter(treatment == 'DMSO'))
# ptest_facet <- ptest + facet_wrap(~facet_col, ncol = 5)
# ptest_facet <- ptest_facet + theme(axis.text.x = element_blank(),
#           axis.title.x = element_blank())
#
# # test_DMSO <- graphing_df_test %>% filter(treatment == 'DMSO', side == 'control')
# # df_temp <- graphing_df %>% filter (treatment == 'DMSO') %>% mutate(side_spd = 0)
# for (i in 1:length(levels(as.factor(graphing_df$treatment)))) {
#   for (j in 1:length(levels(as.factor(graphing_df$side)))) {
#     file_name <- paste('./figs/',levels(as.factor(graphing_df$treatment))[i], '_', levels(as.factor(graphing_df$side))[j], ".png", sep='')
#     # print(file_name)
#     ptest <- plot.windrose(graphing_df %>% filter(treatment == levels(as.factor(graphing_df$treatment))[i], side == levels(as.factor(graphing_df$side))[j]), 'side_spd', 'angle_deg', dirres = 10,  spdseq = c(-1,.5,1), countmax = 4.5)
#     ggsave(filename=file_name, ptest)
#   }
# }
#
# ptest2 <- plot.windrose(graphing_df %>% filter(treatment == 'DMSO'), 'side_spd', 'angle_deg', dirres = 10,  spdseq = c(-1,.5,1), countmax = 4.8)
# file_name <- paste("./figs/DMSO_combined.png")
# ggsave(filename=file_name, ptest)


plot.windrose_adv <- function(data, spd, dir, spdres = .2, dirres = 10,
                          spdmin = -1, spdmax = 1, spdseq = NULL,
                          palette = "Reds", countmax = NA, debug = 0) {

# Look to see what data was passed in to the function
  if (is.numeric(spd) & is.numeric(dir)){
    # assume that we've been given vectors of the speed and direction vectors
    data <- data.frame(spd = spd,
                       dir = dir)
    spd = "spd"
    dir = "dir"
  } else if (exists("data")){
    # Assume that we've been given a data frame, and the name of the speed
    # and direction columns. This is the format we want for later use.
  }

  # Tidy up input data ----
  n.in <- NROW(data)
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))
  data[[spd]][dnu] <- NA
  data[[dir]][dnu] <- NA

  # figure out the wind speed bins ----
  if (missing(spdseq)){
    spdseq <- seq(spdmin,spdmax,spdres)
  } else {
    if (debug >0){
      cat("Using custom speed bins \n")
    }
  }
  # get some information about the number of bins, etc.
  n.spd.seq <- length(spdseq)
  n.colors.in.range <- n.spd.seq - 1

  # create the color map
  spd.colors <- colorRampPalette(brewer.pal(min(max(3,
                                                    n.colors.in.range),
                                                min(9,
                                                    n.colors.in.range)),
                                            palette))(n.colors.in.range)
  spd.colors <- c('white', 'red')
  if (max(data[[spd]],na.rm = TRUE) > spdmax){
    spd.breaks <- c(spdseq,
                    max(data[[spd]],na.rm = TRUE))
    spd.labels <- c(paste(c(spdseq[1:n.spd.seq-1]),
                          '-',
                          c(spdseq[2:n.spd.seq])),
                    paste(spdmax,
                          "-",
                          max(data[[spd]],na.rm = TRUE)))
    spd.colors <- c(spd.colors, "grey50")
  } else{
    spd.breaks <- spdseq
    spd.labels <- paste(c(spdseq[1:n.spd.seq-1]),
                        '-',
                        c(spdseq[2:n.spd.seq]))
  }
  data$spd.binned <- cut(x = data[[spd]],
                         breaks = spd.breaks,
                         labels = spd.labels,
                         ordered_result = TRUE)
  # clean up the data
  data. <- na.omit(data)

  # figure out the wind direction bins
  dir.breaks <- c(-dirres/2,
                  seq(dirres/2, 360-dirres/2, by = dirres),
                  360+dirres/2)
  # # dir.labels <- c(paste(360-dirres/2,"-",dirres/2),
  #                 paste(seq(dirres/2, 360-3*dirres/2, by = dirres),
  #                       "-",
  #                       seq(3*dirres/2, 360-dirres/2, by = dirres)),
  #                 paste(360-dirres/2,"-",dirres/2))
  dir.labels <- c(paste(0),
                paste(seq(dirres, 360 - dirres, by = dirres)),
                paste(0))
  # assign each wind direction to a bin
  dir.binned <- cut(data[[dir]],
                    breaks = dir.breaks,
                    ordered_result = TRUE)
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned
  T_data <- data %>%
    dplyr::group_by(dir.binned) %>%
    dplyr::summarise(count= n()) %>%
    dplyr::mutate(y = count/sum(count))
    labels <- data.frame(x = pi,
                         y = scales::extended_breaks()(range(T_data$y)))
  print(max(T_data$y)*100)

  # Run debug if required ----
  if (debug>0){
    cat(dir.breaks,"\n")
    cat(dir.labels,"\n")
    cat(levels(dir.binned),"\n")
  }

  # deal with change in ordering introduced somewhere around version 2.2
  if(packageVersion("ggplot2") > "2.2"){
    cat("Hadley broke my code\n")
    data$spd.binned = with(data, factor(spd.binned, levels = rev(levels(spd.binned))))
    spd.colors = rev(spd.colors)
  }

  # create the plot ----
  p.windrose <- ggplot(data = data,
                       aes(x = dir.binned, y = (..count..)/sum(..count..)*100,
                           fill = spd.binned, color = spd.binned)) +
    geom_bar(width = 1, size = .5) +
    scale_x_discrete(drop = FALSE,
                     labels = waiver()) +
    coord_polar(start = ((270-(dirres/2))) * pi/180, direction = -1) +
    scale_fill_manual(name = "relative z",
                      values = spd.colors,
                      drop = FALSE) +
    scale_color_manual(name = "relative z",
                  values = c('black','black'),
                  drop = FALSE) +
    theme(axis.title.x = element_blank())

  # adjust axes if required
  if (!is.na(countmax)){
    p.windrose <- p.windrose +
      ylim(c(0,countmax))
  }

  # print the plot
  print(p.windrose)

  # return the handle to the wind rose
  return(p.windrose)
}



# test <- as.data.frame(c(-90, 0, 270, 270, 270, 270))
# colnames(test) <- 'lon'
# test$lat <- c(0, 0, 1, 2, 3, 4)
# test <- st_as_sf(test, coords = c('lon', 'lat'), crs = '+proj=longlat +datum=WGS84')
# test_vec <- project(vect(test), laea)
# plot(test$geometry, axes=TRUE)
# plot(test_vec, axes=TRUE)
# r <- rast(nrow = 9, ncol = 18, xmin=-180, xmax = 180, ymin = -90, ymax=90)
# r <- project(r, laea)
# obs <- rasterize(test_vec, r, fun=function(i){length(i)})
# plot(obs)

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

df_3d_proj <- df_baseline_masked %>% filter(!is.na(unit_x))
df_3d_proj$lat <- (pi/2 - acos(df_3d_proj$unit_z))*180/pi
# df_3d_proj$lon <- (atan2(df_3d_proj$unit_y, df_3d_proj$unit_x))*180/pi # longitude is 'angle'
df_3d_proj$lon <- df_3d_proj$angle*180/pi + 270
df_3d_proj <- df_3d_proj %>% mutate(lon = case_when(lon > 360 ~ lon - 360,
                                                    TRUE ~ lon))
# df_3d_proj <- df_3d_proj %>% mutate(lat = case_when(lat > 180 ~ lat - 180,
#                                                     TRUE ~ lat))

df_3d_proj_sf <- st_as_sf(df_3d_proj, coords = c("lon", "lat"), crs = '+proj=longlat +datum=WGS84 +no_defs +type=crs')
# data <- df_3d_proj_sf %>% filter(treatment == 'U73122' | treatment == 'DMSO', side == 'treated')
# fe_plot <- plot_flatearth(df_3d_proj_sf %>% filter(treatment == 'U73122', side == 'treated'))
#
# p <- plot_diff_flatearth(df_3d_proj_sf %>% filter(treatment == 'DMSO', side == 'control'), df_3d_proj_sf %>% filter(treatment == 'U0126', side == 'treated'))

# sanity check - see how many cells per image
# tally_count <- df_masked %>% filter(!is.na(unit_x)) %>% group_by(sample_info) %>% tally()

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

cdat <- DMSO_control.circ[[1]]
sample_ids <- unique(cellpose_DMSO_control$id)
g <- length(sample_ids)
ndat <- vector("integer", g)
for (i in 1:g) {
  ndat[i] <- nrow(cellpose_DMSO_control[cellpose_DMSO_control$id == sample_ids[i],])
}
gID <- as.integer(as.factor(cellpose_DMSO_control$id))

WalraffTest(cdat,ndat,g,gID)


### BPNME Statistics ###
# This is a much better way of doing the stats for these data, but unfortunately the bpnme takes a VERY long time with
# how many angle measurements per sample we have. It just can't run...
df_LY <- df %>% filter(treatment == 'LY')
df_U0 <- df %>% filter(treatment == 'U0')
df_U73 <- df %>% filter(treatment == 'U73')

# fit.LY <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_LY, its = 100)
fit.LY_corner <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_LY %>% filter((side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx > (2048-800)) | (side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx < 800)), its = 1000, burn = 100, n.lag = 3, seed = 101)
fit.U0_corner <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_U0 %>% filter((side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx > (2048-800)) | (side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx < 800)), its = 1000, burn = 100, n.lag = 3, seed = 101)
fit.U73_corner <- bpnme(pred.I = angle ~ 1 + side_num + (1|id), data=df_U73 %>% filter((side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx > (2048-800)) | (side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx < 800)), its = 1000, burn = 100, n.lag = 3, seed = 101)

# Extract posterior estimates linear coefficients
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

### Old handangle styles ###

cellpose <- df %>% filter(treatment == 'LY' & side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx > (2048 - 800))
cellpose2 <- df %>% filter(treatment == 'LY' & side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx < 800)

cellpose <- df %>% filter(treatment == 'U0' & side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx < (2048-800))
cellpose2 <- df %>% filter(treatment == 'U0' & side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx > 800)

cellpose <- df %>% filter(treatment == 'U73' & side == 'treated' & nuclei_centroidy > 800 & nuclei_centroidx < (2048-800))
cellpose2 <- df %>% filter(treatment == 'U73' & side == 'control' & nuclei_centroidy > 800 & nuclei_centroidx > 800)

# hand.circ <- circular(hand$angle, units = "radians")
cellpose.circ <- circular(cellpose$angle, units = 'radians')
cellpose2.circ <- circular(cellpose2$angle, units = 'radians')

rose.diag(cellpose.circ, prop=15,bins=36, col="red", radii.scale = 'linear', grps = cellpose$id) +
arrows.circular(fit.U0_corner$circ.coef.means[2])

rose.diag(cellpose2.circ, prop=1,bins=36, col="white", radii.scale = 'linear', grps = cellpose2$id) +
arrows.circular(fit.U0_corner$circ.coef.means[1])



rose.diag(hand.circ, prop=8,bins=18, col="red", radii.scale = 'linear') +
arrows.circular(mean(hand.circ), col = 'red')
rose.diag(cellpose.circ, prop=20,bins=36, col="white", radii.scale = 'linear') #+
# arrows.circular(mean(cellpose.circ))
  arrows.circular(fit.U73$circ.coef.means[2])
rose.diag(cellpose2.circ, prop=20,bins=36, col="red", radii.scale = 'linear') +
# arrows.circular(mean(cellpose2.circ), col = 'red')
  arrows.circular(fit.U73$circ.coef.means[1], col = 'red')



treat <- df_LY %>% filter(side == "treated") %>% select(angle)
ctrl <- df_LY %>% filter(side == "control") %>% select(angle)
fitted = fit.LY_adj

# check if it works
treat.mean <- mean(treat[,1])
ctrl.mean <- mean(ctrl[,1])

#circularize
treat.circ <- circular(treat[,1], units = "radians")
ctrl.circ <- circular(ctrl[,1], units = 'radians')

fplot.circular(treat.circ)

treat.circ.mean <- mean(treat.circ)
ctrl.circ.mean <- mean(ctrl.circ)
treat.circ.mean2 <- fitted$circ.coef.means[2]
ctrl.circ.mean2 <- fitted$circ.coef.means[1]

# plot.circular(ctrl.circ, stack=TRUE, pch = 20, sep=0.04)
# arrows.circular(ctrl.circ.mean)
# points(treat.circ,stack=TRUE, col="red", pch= 20, sep = 0.04)
# arrows.circular(treat.circ.mean, col = 'red')

rose.diag(treat.circ, prop=10,bins=18, col="red", radii.scale = 'linear')
arrows.circular(treat.circ.mean2, col = 'red')
rose.diag(ctrl.circ, prop=10,bins=18, col="white", radii.scale = 'linear')
arrows.circular(ctrl.circ.mean2)

for (i in 1:length(levels(dfdbl$image_num))) {
  inum <- levels(dfdbl$image_num)[i]
  treat <- dfdbl %>% filter(side == "tx",ref=="no",image_num ==inum) %>% select(adj.range)
  ctrl <- dfdbl %>% filter(side == "ctrl",ref=="no",image_num ==inum) %>% select(adj.range)

  treat.circ <- circular(treat[,1], units = 'radians')
  ctrl.circ <- circular(ctrl[,1], units = 'radians')

  treat.circ.mean <- mean.circular(treat.circ)
  treat.circ.sd <- sd.circular(treat.circ)
  treat.circ.cr <- 1-angular.variance(treat.circ)/2
  control.circ.mean <- mean.circular(ctrl.circ)
  control.circ.sd <- sd.circular(ctrl.circ)
  control.circ.cr <- 1-angular.variance(ctrl.circ)/2

  png(paste("double bead", inum, "treated.png"))
  rose.diag(treat.circ, prop=3,bins=36, col="red", main=paste(inum, "treated"), xlab=paste("sd:", treat.circ.sd*180/pi, "\ncr:", treat.circ.cr))
  arrows.circular(treat.circ.mean, col = 'red')
  dev.off()

  png(paste("double bead", inum, "control.png"))
  rose.diag(ctrl.circ, prop=3,bins=36,col="white", main = paste(inum, "control"), xlab=paste("sd:", control.circ.sd*180/pi, "\ncr:", control.circ.cr))
  arrows.circular(control.circ.mean)
  dev.off()
}


cellpose_U0126_control <- df %>% filter(treatment == 'U0126' & side == 'control')
cellpose_U0126_treated <- df %>% filter(treatment == 'U0126' & side == 'treated')
U0126_control.circ <- circularize_group(cellpose_U0126_control$angle)
U0126_treated.circ <- circularize_group(cellpose_U0126_treated$angle)

U0126_control_means <- group_sample_circular_stats(U0126_control.circ[[1]], cellpose_U0126_control$id, 'mean')
U0126_treated_means <- group_sample_circular_stats(U0126_treated.circ[[1]], cellpose_U0126_treated$id, 'mean')

U0126_control_medians <- group_sample_circular_stats(U0126_control.circ[[1]], cellpose_U0126_control$id, 'median')
U0126_treated_medians <- group_sample_circular_stats(U0126_treated.circ[[1]], cellpose_U0126_treated$id, 'median')

U0126_control_rho <- group_sample_circular_stats(U0126_control.circ[[1]], cellpose_U0126_control$id, 'rho')
U0126_treated_rho <- group_sample_circular_stats(U0126_treated.circ[[1]], cellpose_U0126_treated$id, 'rho')


U0126_control_mean <- mean.circular(unlist(U0126_control_means))
U0126_control_median <- mean.circular(unlist(U0126_control_medians))
U0126_control_variance <- mean.circular(unlist(U0126_control_rho))
U0126_treated_mean <- mean.circular(unlist(U0126_treated_means))
U0126_treated_median <- mean.circular(unlist(U0126_treated_medians))
U0126_treated_variance <- mean.circular(unlist(U0126_treated_rho))

grouped.rose.diag(U0126_control.circ[1], prop=20,bins=36, col="white", radii.scale = 'linear', grps = cellpose_U0126_control$id, xlab=paste("sd:", U0126_control.circ[[3]]*180/pi, "\ncr:", U0126_control.circ[[4]]))
  arrows.circular(U0126_control_mean, y = rho.circular(unlist(U0126_control_means)))

grouped.rose.diag(U0126_treated.circ[1], prop=20,bins=36, col="red", radii.scale = 'linear', grps = cellpose_U0126_treated$id, xlab=paste("sd:", U0126_treated.circ[[3]]*180/pi, "\ncr:", U0126_treated.circ[[4]]))
  arrows.circular(U0126_treated_mean, y = rho.circular(unlist(U0126_treated_means)))




cellpose_double_control <- df %>% filter(treatment == 'double' & side == 'control')
cellpose_double_treated <- df %>% filter(treatment == 'double' & side == 'treated')
double_control.circ <- circularize_group(cellpose_double_control$angle)
double_treated.circ <- circularize_group(cellpose_double_treated$angle)

grouped.rose.diag(double_control.circ[1], prop=15,bins=36, col="white", radii.scale = 'linear', grps = cellpose_double_control$id, xlab=paste("sd:", double_control.circ[[3]]*180/pi, "\ncr:", double_control.circ[[4]])) +
  arrows.circular(double_control.circ[[2]])
grouped.rose.diag(double_treated.circ[1], prop=15,bins=36, col="red", radii.scale = 'linear', grps = cellpose_double_treated$id, xlab=paste("sd:", double_treated.circ[[3]]*180/pi, "\ncr:", double_treated.circ[[4]])) +
  arrows.circular(double_treated.circ[[2]])


### Prepare the hand-drawn angles
df_handangle <- read.csv('allbeads.csv')
df_handangle <- df_handangle %>% mutate(treatment = case_when(treatment == 'U0126' ~ 'U0', treatment == 'U73122' ~ 'U73',
                                                              treatment == 'LY49002' ~ 'LY', treatment =='doublemix' ~ '2mix'))
df_handangle <- df_handangle %>% mutate(image_num = str_remove(image_num, '^s')) # peel off  the 's' from many of the ids
df_handangle <- dplyr::rename(df_handangle, view = image)
df_handangle <- dplyr::rename(df_handangle, angle_deg = angle)
df_handangle <- dplyr::rename(df_handangle, angle = radians)
df_handangle <- df_handangle %>% mutate(view = case_when(str_detect(view, image_num) ~ '1', TRUE ~ view))
df_handangle <- df_handangle %>% mutate(view = str_remove(view, '\\D+'))
df_handangle$id_old <- substr(df_handangle$image_num, 1,1)
df_handangle$section <- substr(df_handangle$image_num, 2,3)
df_handangle <- df_handangle %>% mutate(section = case_when(treatment == 'LY' ~ '1', TRUE ~ section))
# df_handangle <- df_handangle %>% separate(image_num, c("id_old", "section"), sep = "(?<=\\S)", remove = FALSE)
# df_handangle <- df_handangle %>% mutate(id_old= str_remove(id_old, '\\D+'))
df_handangle <- df_handangle %>% select(!X)  # deletes unneeded column
df_handangle <- df_handangle %>% filter(ref != 'yes')
df_handangle <- df_handangle %>% mutate(side = case_when(side == "ctrl" ~ 'control',
                                         TRUE ~ 'treated'))
df_handangle <- df_handangle %>% mutate(side_num = case_when(side == "control" ~ 0,
                                         TRUE ~ 1))
df_handangle <- df_handangle %>% unite('sample_info', c('treatment', 'id_old', 'side', 'section', 'view'), remove = FALSE)
df_handangle <- df_handangle %>% mutate(sample_info = tolower(sample_info))

df_handangle$id_old <- df_handangle$id
df_handangle$id <- as.factor(df_handangle$id)
df_handangle$id <- as.numeric(df_handangle$id)


df_handangle <- df_handangle %>% mutate(angle = case_when(angle < 0 ~ 2*pi + angle,
                                                          angle > 2*pi ~ 2*pi - angle,
                                                          TRUE ~ angle))