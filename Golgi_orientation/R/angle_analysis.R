#### 0.0 Load libraries for packages ####
# Uncomment and run this the first time to make sure everything gets installed
# install.packages(c('circular', 'plyr', 'dplyr', 'tidyr', 'stringr', 'RImageJROI', 'sf', 'terra', 'ggplot2', 'RColorBrewer', 'Hmisc', 'svglite', 'magick', 'smoothr))
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
library(sfnetworks)
library(lwgeom)
library(Morpho)
library(smoothr)

# needed for 'globe' plots
library(terra)

# graphing
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(svglite)
library(magick)


#### 0.1.1 Data prep helper functions ####
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
  # I don't think we want to flip the 'old' angles since they are now used for positional angle
  # if ('angle_old' %in% names(df_flipped_control)) {
  #   df_flipped_control <- df_flipped_control %>% mutate(angle_old = case_when(angle_old <= pi ~ pi-angle_old, angle_old > pi ~ 3*pi - angle_old))
  # }
  
  df_original_control <- df_to_flip %>% filter((treatment != 'LY' & side_num == 1) |  (treatment == 'LY' & side_num == 0))
  df_flipped_control <- rbind(df_flipped_control, df_original_control)
  return(df_flipped_control)
}


#### 0.1.2 Positional angle and landmark transformation helpers ####
get_transform_matrices_from_rois <- function(info_df, landmark_filenames) {
  transformation_matrix <- NA
  samples <- unique(info_df$old_filename_generic_noside)
  for (sample in 1:length(samples)) {
    matched_lm_files <- landmark_filenames %>% filter(str_starts(landmark_filenames, samples[sample]))
    if (nrow(matched_lm_files) == 3) {
      print(matched_lm_files[1])
      transformation_matrices <- compute_transform_matrices_from_rois(matched_lm_files)
      # this is actually an interesting issue - I used the name 'sample' in my code here for the iterator
      # but it is also a column name in info_df, which makes this code normally not work
      # the fix (if you don't want to rename the variable, which I probably should) is below with the .env$
      info_df <- info_df %>% mutate(transformation_matrix = ifelse(old_filename_generic_noside == samples[.env$sample]
                                                                   & side == 'control', I(list(transformation_matrices[[1]])), 
                                                                   transformation_matrix))
      info_df <- info_df %>% mutate(transformation_matrix = ifelse(old_filename_generic_noside == samples[.env$sample]
                                                                   & side == 'treated', I(list(transformation_matrices[[2]])), 
                                                                   transformation_matrix))
    }
  }
  return(info_df)
}


compute_transform_matrices_from_rois <- function(landmark_filenames) {
  # first load in the three landmark zips
  for (lm_zip in 1:nrow(landmark_filenames)) {
    if (str_detect(landmark_filenames[lm_zip,1], "overview")) {
      overview_lms <- read.ijzip(file.path("./imagej_rois/overview_landmarks/", landmark_filenames[lm_zip,1]), verbose = FALSE)
      # pulls out the coordinates from the zip and flattens them to a vector
      # first two are contralateral, last two are treated side
      B_full <- unlist(lapply(overview_lms, function(l) as.numeric(l$coords)), use.names = FALSE)
    }
    else if (str_detect(landmark_filenames[lm_zip,1], "control")) {
      contra_lms <- read.ijzip(file.path("./imagej_rois/overview_landmarks/", landmark_filenames[lm_zip,1]), verbose = FALSE)
      contra_lms <- unlist(lapply(contra_lms, function(l) as.numeric(l$coords)), use.names = FALSE)
    } 
    else {
      treated_lms <- read.ijzip(file.path("./imagej_rois/overview_landmarks/", landmark_filenames[lm_zip,1]), verbose = FALSE)
      treated_lms <- unlist(lapply(treated_lms, function(l) as.numeric(l$coords)), use.names = FALSE)
    }
  }
  rotation_bool <- FALSE
  if (str_detect(landmark_filenames[lm_zip,1], "U73122_GM130_d7")) {
    rotation_bool <- TRUE
  }
  
  contra_transform_matrix <- solve_transform_matrix_from_lms(B_full[1:4], contra_lms)
  treated_transform_matrix <- solve_transform_matrix_from_lms(B_full[5:8], treated_lms, rotation_bool)
  return(list(contra_transform_matrix, treated_transform_matrix))
}


solve_transform_matrix_from_lms <- function(overview_lms, stack_lms, rotation =FALSE) {
  # need to do some matrix algebra here, but I think this should work
  # there shouldn't be any image rotation so we can solve with two landmarks
  # S is scale
  # T is translation
  # M = [Sx, 0, Tx
  #      0, Sy, Ty
  #      0,  0,  1]
  #
  # [x', y', 1] = [x, y, 1]*M
  # Here the prime xy are landmark data from overview image and non-prime are
  # from z-stack, but we are ignoring z plane.
  #
  # We can solve for M with these equations
  # A = B
  # x1' = x1*Sx + 0*Sy  + 1*Tx + 0*Ty
  # y1' = 0*Sx  + y1*Sy + 0*Tx + 1*Ty
  # x2' = x2*Sx + 0*Sy  + 1*Tx + 0*Ty
  # y2' = 0*Sx  + y2*Sy + 0*Tx + 1*Ty
  # A = [x1, 0, 1, 0
  #      0, y1, 0, 1
  #      x2, 0, 1, 0
  #      0, y2, 0, 1]
  # B = [x1', y1', x2', y2']
  
  # OK so I think morpho just does all this automatically but I don't have the heart
  # to delete all this work

  A <- matrix(c(stack_lms[1], 0, 1, 0,
                       0, stack_lms[2], 0, 1,
                       stack_lms[3], 0, 1, 0,
                       0, stack_lms[4], 0, 1), 4, 4, byrow=TRUE)
  B <- overview_lms

  
  if (rotation) {
    # only one image is rotated. I have to make a different transform matrix
    # luckily the rotation is 90 degrees so the transform is not too complicated
    A <- matrix(c(0, stack_lms[2], 1, 0,
                  -stack_lms[1], 0, 0, 1,
                  0, stack_lms[4], 1, 0,
                    -stack_lms[3], 0, 0, 1), 4, 4, byrow=TRUE)
    M_vals <- solve(A,B)
    
    M <- matrix(c(0, M_vals[2], M_vals[4],
                  -M_vals[1], 0, M_vals[3],
                  0, 0, 1), 3, 3, byrow=TRUE)
  } else {
    M_vals <- solve(A, B)
    M <- matrix(c(M_vals[1], 0, M_vals[3],
                  0, M_vals[2], M_vals[4],
                  0, 0, 1), 3, 3, byrow=TRUE)
  }
  return(M)
}


transform_centroids_xy <- function(temp_df, info_df) {
  # continuing the matrix algebra from above:
  # [x',y',1] =  M * [x, y, 1]
  centroids_df <- data.frame(matrix(NA, nrow = nrow(temp_df), ncol = 2))
  for (image in 1:nrow(info_df)) {
    if (!is.na(info_df[image,]$transformation_matrix)) {
      rows <- which(temp_df$t == info_df[image,]$new_filename)
      M_vals <- info_df[image,]$transformation_matrix
      #pull out the centroids and add the column of 1 to the end for matrix multiplication
      centroids <- temp_df[rows,2:3] %>% mutate(one_col = 1)
      # %*% is matrix multiplication, the t in beginning transposes the results to be in columns instead of rows, the 1 tells it to go row-wise
      transformed_centroids <- t(apply(centroids, 1, function (x) M_vals[[1]] %*% as.numeric(x)))
      centroids_df[rows,1] <- transformed_centroids[,1]
      centroids_df[rows,2] <- transformed_centroids[,2]
    }
  }
  return(centroids_df)
}


transform_reference_centroids_xy <- function(temp_df, info_df) {
  # same idea as above function but this time for the overview to reference cartoon
  # it's easier to have two functions than a bunch of ifelse statements
  centroids_df <- data.frame(matrix(NA, nrow = nrow(temp_df), ncol = 2))
  for (overview_image in 1:length(unique(info_df$old_filename_generic_noside))) {
    info_row <- info_df %>% filter(old_filename_generic_noside == unique(info_df$old_filename_generic_noside)[overview_image]) %>% slice(1)
    M_vals <- info_row$overview_transformation_matrix
    if (!is.na(M_vals)) {
      rows <- which(temp_df$old_filename_generic_noside == info_row$old_filename_generic_noside)
      #pull out the centroids and add the column of 1 to the end for matrix multiplication
      centroids <- temp_df[rows,2:3] %>% mutate(one_col = 1)
      # %*% is matrix multiplication, the t in beginning transposes the results to be in columns instead of rows, the 1 tells it to go row-wise
      transformed_centroids <- t(apply(centroids, 1, function (x) M_vals[[1]] %*% as.numeric(x)))
      centroids_df[rows,1] <- transformed_centroids[,1]
      centroids_df[rows,2] <- transformed_centroids[,2]
    }
  }
  return(centroids_df)
}


get_positional_angle <- function(df_temp, octile_zips) {
  sample_names <- unique(df_temp$old_filename_generic_noside)
  if (is.null(sample_names)) {
    sample_names <- unique(df_temp$old_filename_generic_noside.x)
  }
  for (image in 1:length(sample_names)) {
    octile_zip <- octile_zips %>% filter(str_starts(octile_zips[,1], sample_names[image]))
    if (length(octile_zip[[1]] != 0)) {
      octile_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", octile_zip[[1]]), verbose = FALSE)
      octile_linestrings <- st_sfc(lapply(octile_rois, function(x) st_linestring(x$coords, dim="XY")))
      rows <- which(df_temp$old_filename_generic_noside == sample_names[image] & !is.na(df_temp$nuclei_centroidx_overview))
      if (length(rows) == 0) {
        rows <- which(df_temp$old_filename_generic_noside.x == sample_names[image] & !is.na(df_temp$nuclei_centroidx_overview))
      }
      xmax <- attributes(octile_linestrings)$bbox[['xmax']]
      positional_angles <- data.frame(matrix(ncol=3, nrow=length(rows)))
      print(octile_zip[[1]])
      for (nuc_pair in 1:length(rows)) {
        row <- df_temp[rows[nuc_pair],]
        # print(row)
        extended_line <- extend_line(row, xmax)
        # intersects gives list of number of intersections b/w extended line and each octile, 
        # the which gives the index for the octile that contains intersection
        intersected_segments <- st_intersects(octile_linestrings, extended_line)
        intersected_segment <- which(lapply(intersected_segments, function(x) unlist(x)) != 0)
        if (length(intersected_segment != 0)) {
          intersection_stats <- distance_along_octile(octile_linestrings, extended_line, intersected_segment[1])
          positional_angles[nuc_pair, 1] <- intersection_stats[1][[1]]
          positional_angles[nuc_pair, 2] <- intersection_stats[2][[1]][1]
          positional_angles[nuc_pair, 3] <- intersection_stats[2][[1]][2]
          # print(positional_angle*180/pi)
          # p <- ggplot() + geom_sf(data = octile_linestrings, color = 'black') +
          #   geom_sf(data = octile_linestrings[intersected_segment[1]], color = 'blue') +
          #   geom_point(x = row$nuclei_centroidx_overview, y = row$nuclei_centroidy_overview, aes(color = 'red')) +
          #   geom_sf(data = st_intersection(octile_linestrings, extended_line), color='red') +
          #   geom_sf(data = extended_line)
          # print('pause')
        }
        else {
          print(nuc_pair)
          # p <- ggplot() + geom_sf(data = octile_linestrings, color = 'black') + geom_sf(data = extended_line)
        }
      }
      df_temp[rows,]$positional_angle <- positional_angles[,1]
      df_temp[rows,]$intersectionx <- positional_angles[,2]
      df_temp[rows,]$intersectiony <- positional_angles[,3]
    }
  }
  return(df_temp)
}


extend_line <- function(df_row, xmax) {
  extended_x <- df_row$nuclei_centroidx_overview + df_row$unit_x * 4*xmax 
  extended_y <- df_row$nuclei_centroidy_overview + df_row$unit_y * 4*xmax
  extended_line <- st_linestring(matrix(c(df_row$nuclei_centroidx_overview, df_row$nuclei_centroidy_overview, extended_x, extended_y), 2, 2, byrow=TRUE), dim="XY")
  return(extended_line)
}


distance_along_octile <- function(linestrings, extended_line, linestring_i) {
  # https://stackoverflow.com/a/77688302
  # I'm sure there are other ways of doing this but the sfnetworks library looks 
  # like it will be much simpler!
  
  # putting the first element here in case it hits two octiles. This could maybe happen
  # near boundaries, and they should be close enough together that the angle
  # will be nearly the same anyway.
  intersected_point <- st_intersection(linestrings[linestring_i], extended_line)
  intersected_point <- st_cast(intersected_point, 'POINT')[1]
  octile_network <- as_sfnetwork(linestrings[linestring_i][1])
  subnet <- st_network_blend(octile_network, intersected_point)
  length_table <- subnet %>% activate("edges") %>% st_as_sf() %>% mutate(Length = st_length(x))
  ratio <- pi/4 - (length_table$Length[2] / (length_table$Length[1]+length_table$Length[2]) * pi/4)
  ratio <- pi + pi/4*(linestring_i-1) + ratio
  return(c(ratio,intersected_point))
}


positional_angle_to_xy <- function(linestrings, pos_angle) {
  # this is basically the inverse of the 'distance along octile' function
  # we want to get the xy position in a linestring from the angle
  # https://stackoverflow.com/a/77688302
  
  # would this be better as just a lookup table? IDK
  # check 'approx()' function. This code is slow 
  # aight I changed it and it's literally 1000x faster lol
  
  # this is a bit obfuscated - the positional angle plus pi rotates it cw so '0' is 180
  # but now all the 180-360 will be 360-540, so we modulo with a full circle 
  # so that they will be 0-180 instead. Divide by pi/4 (1/8 of circle) to get the octile
  # and floor it so it's an integer not a float, then add 1 since it is 1-8 not 0-7
  
  linestring_i <- floor(((pos_angle+pi) %% (2*pi)) / (pi/4)) + 1
  
  ratio <- (pos_angle %% (pi/4)) / (pi/4)
  (pt <- st_linesubstring(linestrings[linestring_i], from = 0, to = ratio) %>% st_endpoint())
  return(pt)
}


generate_overview_positional_LUT <- function(overview_octile_rois, slices = 720) {
  overview_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", overview_octile_rois), verbose = FALSE)
  linestrings <- st_sfc(lapply(overview_rois, function(x) st_linestring(x$coords, dim="XY")))
  
  # https://stackoverflow.com/a/72533271
  # https://stackoverflow.com/a/72267454
  bbox <- st_bbox(linestrings) %>% st_as_sfc()
  
  polygon <- bbox %>% lwgeom::st_split(linestrings) %>% st_collection_extract("POLYGON")
  poly_lms <- polygon[2] # hopefully this isn't random!
  smooth_poly <- smoothr::smooth(poly_lms, method = "ksmooth", smoothness = 20)
  smooth_linestring <- st_cast(smooth_poly, 'LINESTRING')
  endpts <- st_line_sample(linestrings, sample=0)
  st_nearest_points(endpts, smooth_linestring) %>% {. ->> connecting_linestrings}
  new_endpts <- st_line_sample(connecting_linestrings, sample = 1)
  net <- as_sfnetwork(smooth_linestring)
  net <- st_network_blend(net, st_cast(new_endpts, 'POINT'))
  net <- convert(net, to_spatial_smooth, protect = new_endpts)
  smooth_linestrings <- net %>% activate(edges) %>% st_as_sf()
  smooth_linestrings <- st_sfc(smooth_linestrings[[4]])
  
  
  LUT <- setNames(data.frame(matrix(ncol = 3, nrow = slices)), c('pos_angle', 'overview_x', 'overview_y'))
  LUT$pos_angle <- seq(0, 2*pi, length.out = slices)
  
  for (slice in 1:slices) {
    pt_temp <- positional_angle_to_xy(smooth_linestrings, LUT[slice,]$pos_angle)
    LUT[slice,]$overview_x <- pt_temp[[1]][1]
    LUT[slice,]$overview_y <- pt_temp[[1]][2]
  }
  return(LUT)
}


convert_directional_angle_overview_LUT <- function(df_temp, overview_pos_LUT) {
  df_temp$overview_intersectionx <- NA
  df_temp$overview_intersectiony <- NA

  df_temp$overview_intersectionx <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_x, xout = df_temp$positional_angle)$y
  df_temp$overview_intersectiony <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_y, xout = df_temp$positional_angle)$y
  return(df_temp)
} 


get_positional_angle_from_intersection <- function(df_temp, octile_zip, extended_linestrings) {
  octile_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", octile_zip[[1]]), verbose = FALSE)
  octile_linestrings <- st_sfc(lapply(octile_rois, function(x) st_linestring(x$coords, dim="XY")))
  intersection_stats_df <- data.frame(nrow=nrow(df_temp), ncol = 3)
  for (row in 1:nrow(df_temp)) {
    intersected_segments <- st_intersects(octile_linestrings, extended_linestrings[row])
    intersected_segment <- which(lapply(intersected_segments, function(x) unlist(x)) != 0)
    if (length(intersected_segment != 0)) {
      intersection_stat <- distance_along_octile(octile_linestrings, extended_linestrings[row], intersected_segment[1])
      intersection_stats_df[row,1] <- intersection_stat[1][[1]]
      intersection_stats_df[row,2] <- intersection_stat[2][[1]][1]
      intersection_stats_df[row,3] <- intersection_stat[2][[1]][2]
      # print(positional_angle*180/pi)
      # p <- ggplot() + geom_sf(data = octile_linestrings, color = 'black') +
      #   geom_sf(data = octile_linestrings[intersected_segment[1]], color = 'blue') +
      #   geom_point(x = row$nuclei_centroidx_overview, y = row$nuclei_centroidy_overview, aes(color = 'red')) +
      #   geom_sf(data = st_intersection(octile_linestrings, extended_line), color='red') +
      #   geom_sf(data = extended_line)
      # print('pause')
    }
  }
  return(intersection_stats_df)
}


compute_transform_matrices_from_pos_angles <- function(temp_df, min_dist) {
  found_matches <- FALSE
  while (!found_matches) {
    pair = sample.int(nrow(temp_df), 2)
    x_dist = temp_df[pair[1],]$overview_intersectionx - temp_df[pair[2],]$overview_intersectionx
    y_dist = temp_df[pair[1],]$overview_intersectiony - temp_df[pair[2],]$overview_intersectiony
    distance = sqrt((x_dist)^2 + (y_dist)^2)
    if (!is.na(distance)) {
      if ((distance >= min_dist) & (x_dist >= 0) & (y_dist >= 80)) {
        found_matches <- TRUE
      }  
    }
  }

  model_overview_landmarks <- c(temp_df[pair[1],]$overview_intersectionx, temp_df[pair[1],]$overview_intersectiony,
                          temp_df[pair[2],]$overview_intersectionx, temp_df[pair[2],]$overview_intersectiony)
  scout_overview_landmarks <- c(temp_df[pair[1],]$intersectionx, temp_df[pair[1],]$intersectiony,
                                temp_df[pair[2],]$intersectionx, temp_df[pair[2],]$intersectiony)
  
  
  overview_transform_matrix <- solve_transform_matrix_from_lms(model_overview_landmarks, scout_overview_landmarks)
  return(overview_transform_matrix)
}


get_overview_transform_matrices_from_pos_angles <- function(info_df, temp_df, min_dist = 100) {
  transformation_matrix <- NA
  info_df$overview_transformation_matrix <- NA 
  samples <- unique(info_df$old_filename_generic_noside)
  for (sample in 1:length(samples)) {
    calc_transformation_matrix <- compute_transform_matrices_from_pos_angles(temp_df %>% filter(old_filename_generic_noside == unique(info_df$old_filename_generic_noside)[sample]), min_dist)
    # matched_lm_files <- landmark_filenames %>% filter(str_starts(landmark_filenames, samples[sample]))
    # this is actually an interesting issue - I used the name 'sample' in my code here for the iterator
    # but it is also a column name in info_df, which makes this code normally not work
    # the fix (if you don't want to rename the variable, which I probably should) is below with the .env$
    info_df <- info_df %>% mutate(overview_transformation_matrix = ifelse(old_filename_generic_noside == samples[.env$sample], I(list(calc_transformation_matrix)), overview_transformation_matrix))
  }
  return(info_df)
}


plot_extended_lines <- function(filtered_df, octile_files) {
  #load octiles
  
  p <- ggplot() + geom_sf(data = octiles, color = 'black') +
    geom_point(x = row$nuclei_centroidx_overview, y = row$nuclei_centroidy_overview, aes(color = 'red')) +
    geom_sf(data = st_intersection(octile_linestrings, extended_line), color='red') +
    geom_sf(data = extended_line)
}


morph_centroids_to_ref_img <- function(temp_df, overview_pos_LUT, overview_octile_rois, info_df) {
  temp_df$morphedx <- NA
  temp_df$morphedy <- NA
  
  # Let's get a few internal lms along the midline b/w nasalpits, 
  # to match from ones we got when making regions
  # we will do outside of loop since it is from the reference not sample
  overview_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", overview_octile_rois), verbose = FALSE)
  linestrings <- st_sfc(lapply(overview_rois, function(x) st_linestring(x$coords, dim="XY")))
  midline <- st_linestring(rbind(st_coordinates(st_line_sample(linestrings[3], sample=0)), st_coordinates(st_line_sample(linestrings[7], sample=0)))[,1:2])
  # same for nasalpit, then get the pt where the 2 lines intersect (midpoint)
  np_line <- st_linestring(rbind(st_coordinates(st_line_sample(linestrings[1], sample=0)), st_coordinates(st_line_sample(linestrings[4], sample=1)))[,1:2])
  np_midpoint <- st_cast(st_intersection(midline, np_line), 'POINT') 
  contra_np_line <- st_linestring(rbind(st_coordinates(st_line_sample(linestrings[1], sample=0))[1:2], st_coordinates(np_midpoint)))
  treat_np_line <- st_linestring(rbind(st_coordinates(st_line_sample(linestrings[4], sample=1))[1:2], st_coordinates(np_midpoint)))
  
  # now make all the 2/5 2/5 1/5 lines for np...
  np_endpts_ref <- data.frame(matrix(nrow = 5, ncol=2))
  np_endpts_ref[1,] <- st_coordinates(st_linesubstring(contra_np_line, from = 0, to = .4) %>% st_endpoint())
  np_endpts_ref[2,] <- st_coordinates(st_linesubstring(contra_np_line, from = .4, to = .8) %>% st_endpoint())
  np_endpts_ref[3,] <- st_coordinates(st_linesubstring(treat_np_line, from = .8, to = 1) %>% st_startpoint())
  np_endpts_ref[4,] <- st_coordinates(st_linesubstring(treat_np_line, from = .8, to = .4) %>% st_endpoint())
  np_endpts_ref[5,] <- np_midpoint
  
  for (sample_i in 1:length(unique(temp_df$old_filename_generic_noside))) {
    sample_df <- temp_df %>% filter(old_filename_generic_noside == unique(temp_df$old_filename_generic_noside)[sample_i])
    print(unique(temp_df$old_filename_generic_noside)[sample_i])
    num_lms <- nrow(sample_df) / 5
    pull_rows <- as.integer(seq(from=1, to=nrow(sample_df), length.out=num_lms))
    landmarks <- sample_df %>% arrange(positional_angle) %>% slice(pull_rows) %>% select(positional_angle, intersectionx, intersectiony)
    landmarks$referencex <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_x, xout = landmarks$positional_angle)$y
    landmarks$referencey <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_y, xout = landmarks$positional_angle)$y
    
    # this is all code for getting interior landmarks, it doesn't work will so I've commented it out, not sure it's worth figuring out
    # # get random cells
    # landmarks2 <- sample_df %>% na.omit(nuclei_centroidx_overview) %>% slice(sample.int(nrow(sample_df %>% na.omit(nuclei_centroidx_overview)), num_lms)) %>% select(nuclei_centroidx_overview, nuclei_centroidy_overview)
    # # make line b/w overview xy and cell xy
    # lm_pair_lines <- st_sfc(mapply(function(x) st_linestring(matrix(c(landmarks[x,2], landmarks[x,3], landmarks2[x,1], landmarks2[x,2]), 2,2, byrow=TRUE), dim='XY'), 1:num_lms, SIMPLIFY = FALSE))
    # 
    # # extend the line and find intersect on other side of line
    # landmarks2$og_length <- st_length(lm_pair_lines) # get lengths
    # # gets unit vector, multiplies it by a big desired length, and add it to original xy coords
    # landmarks2$extendedx <- landmarks[,2] + ((landmarks2[,1] - landmarks[,2]) / landmarks2[,3]) * max(landmarks[,2])
    # landmarks2$extendedy <- landmarks[,3] + ((landmarks2[,2] - landmarks[,3]) / landmarks2[,3]) * max(landmarks[,2])
    # 
    # lm_extended_lines <- st_sfc(mapply(function(x) st_linestring(matrix(c(landmarks2[x,1], landmarks2[x,2], landmarks2[x,4], landmarks2[x,5]), 2,2, byrow=TRUE), dim='XY'), 1:num_lms, SIMPLIFY = FALSE))
    # 
    # octile_zip <- octile_filenames %>% filter(str_starts(octile_filenames, sample_df$old_filename_generic_noside[1]))
    # intersection_stats <- get_positional_angle_from_intersection(landmarks2, octile_zip, lm_extended_lines)
    # landmarks2$extended_pos_angle <- intersection_stats[,1]
    # landmarks2$extended_intersectx <- intersection_stats[,2]
    # landmarks2$extended_intersecty <- intersection_stats[,3]
    # # clear out NA values
    # to_del_rows <- which(is.na(landmarks2$extended_pos_angle))
    # landmarks2 <- landmarks2[-to_del_rows,]
    # landmarks <- landmarks[-to_del_rows,]
    # 
    # # find length ratio
    # lm_to_lm_lines <- st_sfc(mapply(function(x) st_linestring(matrix(c(landmarks[x,]$intersectionx, landmarks[x,]$intersectiony, landmarks2[x,]$extended_intersectx, landmarks2[x,]$extended_intersecty), 2,2, byrow=TRUE), dim='XY'), 1:nrow(landmarks2), SIMPLIFY = FALSE))
    # landmarks2$length_ratio <- landmarks2$og_length / st_length(lm_to_lm_lines)
    # 
    # # draw line on reference b/w reference xys
    # landmarks2$extended_pos_angle <- as.numeric(landmarks2$extended_pos_angle)
    # landmarks2 <- landmarks2 %>% mutate(extended_pos_angle = case_when(extended_pos_angle > 2*pi ~ extended_pos_angle - 2*pi,
    #                                                                    extended_pos_angle < 0 ~ extended_pos_angle + 2*pi,
    #                                                                TRUE ~ extended_pos_angle))
    # 
    # landmarks2$reference_intersectx <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_x, xout = landmarks2$extended_pos_angle)$y
    # landmarks2$reference_intersecty <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_y, xout = landmarks2$extended_pos_angle)$y
    # ref_to_ref_lines <- st_sfc(mapply(function(x) st_linestring(matrix(c(landmarks[x,4], landmarks[x,5], landmarks2[x,]$reference_intersectx, landmarks2[x,]$reference_intersecty), 2,2, byrow=TRUE), dim='XY'), 1:nrow(landmarks2), SIMPLIFY = FALSE))
    # 
    # ref_intercepts <- data.frame(matrix(NA, nrow = nrow(landmarks2), ncol=2))
    # for (line_i in 1:nrow(landmarks2)) {
    #   (pt <- st_linesubstring(ref_to_ref_lines[line_i], from = 0, to = landmarks2[line_i,]$length_ratio) %>% st_endpoint())
    #   ref_intercepts[line_i,1] <- pt[[1]][1]
    #   ref_intercepts[line_i,2] <- pt[[1]][2]
    # }
    # 
    # # cull interior landmarks that fall outside of perimeter lms
    # # https://stackoverflow.com/a/52671103
    # overview_outline <- st_as_sf(landmarks, coords = c("referencex","referencey"), remove = FALSE)
    # overview_outline <- overview_outline %>% summarise(geometry = st_combine(geometry)) %>% st_cast("POLYGON")
    # intersection_pts <- st_intersects(st_as_sf(ref_intercepts, coords=c(1:2)), overview_outline, sparse=FALSE)
    # 
    # # line_i = 8
    # # test <- st_linesubstring(ref_to_ref_lines[line_i], from = 0, to = landmarks2[line_i,]$length_ratio)
    # # test_pt <- st_linesubstring(ref_to_ref_lines[line_i], from = 0, to = landmarks2[line_i,]$length_ratio) %>% st_endpoint()
    # # plot1<- ggplot()  + geom_point(aes(x = landmarks[,4], y = landmarks[,5]), color='red') + geom_sf(data = ref_to_ref_lines[line_i]) + geom_sf(data=test, color = 'red') + geom_sf(data=test_pt, color = 'red')
    # # plot2 <- ggplot() + geom_point(aes(x = landmarks[,2], y = landmarks[,3])) + geom_sf(data = lm_to_lm_lines[line_i]) + geom_point(aes(x = landmarks2[line_i,1], y = landmarks2[line_i,2]), color='red')
    # # plot_grid(plot1, plot2)
    # 
    # 
    # plot1<- ggplot()  + geom_point(aes(x = full_landmarks[,1], y = full_landmarks[,2]), color='red')
    # plot2 <- ggplot() +  geom_point(aes(x = full_landmarks[,3], y = full_landmarks[,4])) 
    # plot_grid(plot1, plot2)
    
    
    # frustrating you  can't rbind if the column names are different
    perimeter_lms <- cbind(landmarks[4:5],landmarks[2:3])
    np_interior_lms <- info_df[info_df$old_filename_generic_noside == unique(temp_df$old_filename_generic_noside)[sample_i],]$internal_lms[[1]]
    np_interior_lms <- cbind(np_endpts_ref, np_interior_lms)
    # colnames(perimeter_lms) <- c('referencex', 'referencey', 'overviewx', 'overviewy')
    # interior_lms <- cbind(ref_intercepts[intersection_pts,], landmarks2[intersection_pts,1:2])
    # colnames(interior_lms) <- colnames(perimeter_lms)
    # full_landmarks <- rbind(perimeter_lms, interior_lms)
    colnames(np_interior_lms) <- colnames(perimeter_lms)
    full_landmarks <- rbind(perimeter_lms, np_interior_lms)

    rotateM <- getTrafo4x4(rotonto(data.matrix(full_landmarks[1:2]), data.matrix(full_landmarks[3:4]), scale=TRUE))
    centered_pts <- applyTransform(data.matrix(full_landmarks[3:4]), rotateM)
    # make TPS
    morphM <- computeTransform(data.matrix(full_landmarks[1:2]), centered_pts, type='tps')
    # apply warp
    rot_pts <- applyTransform(data.matrix(sample_df[26:27]), rotateM)
    morphed_pts <- applyTransform(rot_pts, morphM)
    edit_rows <- which(temp_df$old_filename_generic_noside == unique(temp_df$old_filename_generic_noside)[sample_i])
    temp_df[edit_rows,]$morphedx <- morphed_pts[,1]
    temp_df[edit_rows,]$morphedy <- morphed_pts[,2]
  }
  return(temp_df)
}


calc_unit_extend_line <- function(start_pts, end_pts, mult) {
  # gets unit vector, multiplies it by a big desired length, and add it to original xy coords
  original_length <- sqrt((end_pts[1] - start_pts[1])^2 + (end_pts[2] - start_pts[2])^2)
  extendedx <- start_pts[1] + ((end_pts[1] - start_pts[1]) / original_length) * mult
  extendedy <- start_pts[2] + ((end_pts[2] - start_pts[2]) / original_length) * mult
  return(data.frame(extendedx,extendedy))
}


assign_tissue_region <- function(temp_df, info_df, octile_filenames) {
  # I think setting up 5 'equally' spaced regions will be good
  # then can choose a distance away from baseline and constrain to the region
  # glob, mid, center, mid, glob
  
  # I will draw 2/5, 2/5, 1/5 on each of the 2 baseline octiles and combine center 1/5ths
  # (don't want to do 1/5 across the entire baseline since the sections aren't 
  # perfectly 'flat' and the treated side is often smaller than the other) 
  # then draw line b/w top of the nasalpit lines and a b/w FNP midline and brain midline
  # Break the nasalpit line in half at intersection and make the 2/5, 2/5, 1/5 on
  # it and draw lines b/w baseline and nasalpit lines, extend them to back of
  # the tissue and these will makeup the regions
  
  # this code is gonna be pretty much the same as the previous 'masking' fncs above
  
  #TODO maybe remake 3mix rois... It's very difficult to determine midline
  temp_df$region <- NA
  info_df$internal_lms <- NA
  for (sample_i in 1:length(unique(temp_df$old_filename_generic_noside))) {
    sample_df_rows <- which(temp_df$old_filename_generic_noside == unique(temp_df$old_filename_generic_noside)[sample_i])
    sample_df <- temp_df[sample_df_rows,]
    octile_zip <- octile_filenames[,1][str_starts(octile_filenames[,1], unique(temp_df$old_filename_generic_noside)[sample_i])]
    print(sample_i)
    if (length(octile_zip) != 0) {
      octile_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", octile_zip[[1]]), verbose = FALSE)
      octile_linestrings <- st_sfc(lapply(octile_rois, function(x) st_linestring(x$coords, dim="XY")))
      # this is bigger than the actual sample but it shouldn't matter as cells aren't outside of the the tissue...
      FNP_poly <- st_convex_hull(st_union(octile_linestrings)) 
      # the baseline are linestring 2 and 3, the nasal pit is 1 and 4, brain is 6 and 7
      # make line from baseline mid to brain mid - I doubt this is the 'right' way to do this...
      midline <- st_linestring(rbind(st_coordinates(st_line_sample(octile_linestrings[3], sample=0)), st_coordinates(st_line_sample(octile_linestrings[7], sample=0)))[,1:2])
      # same for nasalpit, then get the pt where the 2 lines intersect (midpoint)
      np_line <- st_linestring(rbind(st_coordinates(st_line_sample(octile_linestrings[1], sample=0)), st_coordinates(st_line_sample(octile_linestrings[4], sample=1)))[,1:2])
      np_midpoint <- tryCatch( 
        {st_cast(st_intersection(midline, np_line), 'POINT') 
        
        }, error=function(e) {
        print('sample midline is too short!')
        midline_ext_pts <- calc_unit_extend_line(st_coordinates(st_line_sample(octile_linestrings[3], sample=0)), st_coordinates(st_line_sample(octile_linestrings[7], sample=0)), st_length(midline)*2)
        midline_ext <- st_linestring(rbind(st_coordinates(st_line_sample(octile_linestrings[3], sample=0))[1:2], unlist(midline_ext_pts)[1:2])[,1:2])
        st_cast(st_intersection(midline_ext, np_line), 'POINT')
        }
      )
      contra_np_line <- st_linestring(rbind(st_coordinates(st_line_sample(octile_linestrings[1], sample=0))[1:2], st_coordinates(np_midpoint)))
      treat_np_line <- st_linestring(rbind(st_coordinates(st_line_sample(octile_linestrings[4], sample=1))[1:2], st_coordinates(np_midpoint)))
      
      # breakup midline linestrings
      bl_net <- as_sfnetwork(octile_linestrings[2:3])
      bl_endpts <- data.frame(matrix(nrow = 4, ncol=2))
      bl_endpts[1,] <- st_coordinates(st_linesubstring(octile_linestrings[2], from = 0, to = .4) %>% st_endpoint())
      bl_endpts[2,] <- st_coordinates(st_linesubstring(octile_linestrings[2], from = .4, to = .8) %>% st_endpoint())
      bl_endpts[3,] <- st_coordinates(st_linesubstring(octile_linestrings[3], from = 0, to = .2) %>% st_endpoint())
      bl_endpts[4,] <- st_coordinates(st_linesubstring(octile_linestrings[3], from = .2, to = .6) %>% st_endpoint())
      bl_net <- st_network_blend(bl_net, st_as_sf(bl_endpts, coords=c(1:2)))
      
      # now make all the 2/5 2/5 1/5 lines for np...
      np_endpts <- data.frame(matrix(nrow = 4, ncol=2))
      np_endpts[1,] <- st_coordinates(st_linesubstring(contra_np_line, from = 0, to = .4) %>% st_endpoint())
      np_endpts[2,] <- st_coordinates(st_linesubstring(contra_np_line, from = .4, to = .8) %>% st_endpoint())
      np_endpts[3,] <- st_coordinates(st_linesubstring(treat_np_line, from = .8, to = 1) %>% st_startpoint())
      np_endpts[4,] <- st_coordinates(st_linesubstring(treat_np_line, from = .8, to = .4) %>% st_endpoint())

      
      # #make first lines
      # contra_glob_mid2np <- st_linestring(matrix(c(unlist(bl_endpts[1,]), unlist(np_endpts[1,])), 2,2, byrow=TRUE))
      # contra_mid_mid2np <- st_linestring(matrix(c(unlist(bl_endpts[2,]), unlist(np_endpts[2,])), 2,2, byrow=TRUE))
      # treat_mid_mid2np <- st_linestring(matrix(c(unlist(bl_endpts[3,]), unlist(np_endpts[3,])), 2,2, byrow=TRUE))
      # treat_glob_mid2np <- st_linestring(matrix(c(unlist(bl_endpts[4,]), unlist(np_endpts[4,])), 2,2, byrow=TRUE))
      
      # extend the lines
      extended_endpts <- apply(cbind(bl_endpts, np_endpts), 1, function(x) unlist(calc_unit_extend_line(x[1:2], x[3:4], 2*max(sample_df$nuclei_centroidy_overview))))
      extended_endpts <- t(extended_endpts)
      extended_startpts <- apply(cbind(np_endpts, bl_endpts), 1, function(x) unlist(calc_unit_extend_line(x[1:2], x[3:4], max(sample_df$nuclei_centroidy_overview))))
      extended_startpts <- t(extended_startpts)
      # lm_pair_lines <- st_sfc(mapply(function(x) st_linestring(matrix(c(landmarks[x,2], landmarks[x,3], landmarks2[x,1], landmarks2[x,2]), 2,2, byrow=TRUE), dim='XY'), 1:num_lms, SIMPLIFY = FALSE))
      
      cut_lines <- st_sfc(lapply(1:nrow(extended_endpts), function(x) st_linestring(matrix(cbind(extended_startpts, extended_endpts)[x,], 2,2, byrow=TRUE))))
      # compared to other things in sf this almost feels too easy...
      # FNP_region_polys <- st_collection_extract(st_split(FNP_poly, cut_lines))
      # lmao it was too easy, the order of the shapes will not always be the same
      # There's not many regions so we'll just do some embarrassing loops
      # This is so hacky and won't work with a lot of scenarios, but for this
      # dataset it works fine
      # regions <- FNP_poly
      areas <- data.frame(matrix(nrow = 5, ncol=1))
      for (region_i in 1:5) {
        if (region_i == 1) {
          # contra glob
          temp_polys <-  st_collection_extract(st_split(FNP_poly, cut_lines[1]))
          temp_areas <- st_area(temp_polys)
          regions <- st_sfc(temp_polys[which.min(temp_areas)])
          areas[1,] <- min(temp_areas)
        }  else if (region_i == 2) {
          # treat glob
          temp_polys <-  st_collection_extract(st_split(FNP_poly, cut_lines[4]))
          temp_areas <- st_area(temp_polys)
          regions <- st_sfc(rbind(regions, temp_polys[which.min(temp_areas)][1]))
          areas[5,] <- min(temp_areas)
        } else if (region_i == 3) {
          # contra mid
          temp_polys <-  st_collection_extract(st_split(FNP_poly, cut_lines[1:2]))
          temp_areas <- st_area(temp_polys)
          areas[2,] <- min(temp_areas[-match(areas[1,], temp_areas)])
          regions <- st_sfc(rbind(regions, temp_polys[which(temp_areas == areas[2,])]))
        } else if (region_i == 4) {
          # treat mid
          temp_polys <-  st_collection_extract(st_split(FNP_poly, cut_lines[3:4]))
          temp_areas <- st_area(temp_polys)
          areas[4,] <- min(temp_areas[-match(areas[5,], temp_areas)])
          regions <- st_sfc(rbind(regions, temp_polys[which(temp_areas == areas[4,])]))
        } else if (region_i == 5) {
          #center
          temp_polys <-  st_collection_extract(st_split(FNP_poly, cut_lines[2:3]))
          temp_areas <- st_area(temp_polys)
          # this will very likely be the smallest area, but if not we can add some code
          # similar to here. It needs protection from float errors
          # temp_areas[-match(round(areas[1,]+areas[2,]), round(temp_areas))]
          areas[3,] <- min(temp_areas)
          regions <- st_sfc(rbind(regions, temp_polys[which(temp_areas == areas[3,])]))
        }
      }
      regions <- regions[order(c(1,5,2,4,3))]
      centroids <- sample_df %>% select(nuclei_centroidx_overview, nuclei_centroidy_overview)
      # now here is the kind of sf pain I would expect. the 'last' command just reduces the returned index to one value
      # as there are plenty of cells on the border b/w polygons
      intersection_indices <- st_intersects(st_cast(st_sfc(st_multipoint(data.matrix(centroids))), 'POINT'), regions) %>% sapply(last)
      temp_df[sample_df_rows,]$region <- intersection_indices
      
      np_endpts[5,] <- np_midpoint
      info_df[info_df$old_filename_generic_noside == sample_df[1,]$old_filename_generic_noside,]$internal_lms <- list(np_endpts)
    }
  }
  temp_df$region_name <- NA
  temp_df <- temp_df %>% mutate(region_name = case_when((region == 1 | region == 5) ~ 'glob',
                                                         (region == 2 | region == 4) ~ 'mid',
                                                         region == 3 ~ 'center', TRUE ~ region_name))
  return(list(temp_df, info_df))
}



filter_distance_from_octile_lines <- function(temp_df, octile_linestrings, info_df, low_thresh, high_thresh, octiles=c(2,3)) {
  for (sample_i in 1:length(unique(temp_df$old_filename_generic_noside))) {
    sample_df_rows <- which(temp_df$old_filename_generic_noside == unique(temp_df$old_filename_generic_noside)[sample_i])
    sample_df <- temp_df[sample_df_rows,]
    octile_zip <- octile_filenames[,1][str_starts(octile_filenames[,1], unique(temp_df$old_filename_generic_noside)[sample_i])]
    if (length(octile_zip) != 0) {
      octile_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", octile_zip[[1]]), verbose = FALSE)
      octile_linestrings <- st_sfc(lapply(octile_rois, function(x) st_linestring(x$coords, dim="XY")))
      base_linestrings <- octile_linestrings[octiles]
      # this is another that's gonna be a bit hacky... technically the first element
      # of transformation matrix is the scale b/w image and 'overview' image. Since the 
      # scale should be square the first scale term should be fine
      # .28 is um per pixel from zstack,
      scale <- 0.2840910 / abs(info_df[info_df$old_filename_generic_noside == unique(temp_df$old_filename_generic_noside)[sample_i],]$transformation_matrix[[1]][1,1])
      distance_low <- low_thresh/scale
      distance_high <- high_thresh/scale
      
      centroids <- sample_df %>% select(nuclei_centroidx_overview, nuclei_centroidy_overview)
      high_rows <- st_intersects(st_cast(st_sfc(st_multipoint(data.matrix(centroids), dim='XY')), 'POINT'), st_union(st_buffer(base_linestrings, distance_high)))
      keep_rows <- data.frame(high_rows)
      if (distance_low != 0) {
        inner_rows <- st_intersects(st_cast(st_sfc(st_multipoint(data.matrix(centroids), dim='XY')), 'POINT'), st_union(st_buffer(base_linestrings, distance_low)))
        inner_rows <- data.frame(inner_rows)
        keep_rows <- keep_rows[!keep_rows$row.id %in% inner_rows$row.id,]
      }
      if (sample_i == 1) {
        filter_df <- sample_df[keep_rows$row.id,]
      } else {
        filter_df <- rbind(filter_df, sample_df[keep_rows$row.id,])
      }
    }
  }
  return(filter_df)
}


# get_octile_endpoint_network <- function(octile_filename) {
#   # load up the overview octile, make the rois into linestrings
#   octile_rois <- read.ijzip(file.path("./imagej_rois/overview_octiles/", octile_filename), verbose = FALSE)
#   octile_linestrings <- st_sfc(lapply(octile_rois, function(x) st_linestring(x$coords, dim="XY")))
#   
#   # pull out the roi endpoing only, we don't want the first point as they are same as endpoint
#   octile_endpoints <- lapply(octile_linestrings, function(x) tail(st_coordinates(x), 1)[1:2])
#   # find all combinations of these endpoint pairs
#   combinations <- combn(1:length(octile_endpoints), 2)
#   # make linestrings from each combination
#   octile_endpoint_mesh <- st_sfc(mapply(function(x) st_linestring(c(st_point(c(octile_endpoints[[combinations[,x][1]]][1], octile_endpoints[[combinations[,x][1]]][2])), 
#                                                                     st_point(c(octile_endpoints[[combinations[,x][2]]][1], octile_endpoints[[combinations[,x][2]]][2]))), dim="XY"), 1:ncol(combinations), SIMPLIFY = FALSE))
#   # make it into a network
#   octile_endpoint_network <- as_sfnetwork(octile_endpoint_mesh)
#   # find all the points where the lines intersect and make them into points
#   intersection_ps <- st_collection_extract(st_intersection(octile_endpoint_mesh), "POINT")
#   # # add these points back into the network, which will remake the lines
#   octile_endpoint_network <- st_network_blend(octile_endpoint_network, intersection_ps)
#   return(octile_endpoint_network)
# }
# 
# sample_net <- get_octile_endpoint_network(octile_filenames[1,])
# reference_net <- get_octile_endpoint_network(octile_filenames[26,])
# 
# # The edges are not connected.
# as_sfnetwork(edges)
# 
# # 2 apply rotation matrix, this will do linear transformation of the data to the reference overview
# # we will need to do a thin plate spline warp to get it really aligned after
# # probably could do the warp without all this other algebra... but it feels like too far to warp...
# # and I already had to figure out how to do the transformation algebra so I'm leaving it in
# # apply transformation matrix
# df_masked$nuclei_centroidx_reference <- NA
# df_masked$nuclei_centroidy_reference <- NA
# df_masked$GM130_centroidx_reference <- NA
# df_masked$GM130_centroidy_reference <- NA
# df_masked$reference_intersectionx <- NA
# df_masked$reference_intersectiony  <- NA
# 
# transformed_centroids <- transform_reference_centroids_xy(df_masked %>% select(old_filename_generic_noside, nuclei_centroidx_overview, nuclei_centroidy_overview), info_df)
# df_masked$nuclei_centroidx_reference <- transformed_centroids[,1]
# df_masked$nuclei_centroidy_reference <- transformed_centroids[,2]
# 
# transformed_centroids <- transform_reference_centroids_xy(df_masked %>% select(old_filename_generic_noside, GM130_centroidx_overview, GM130_centroidy_overview), info_df)
# df_masked$GM130_centroidx_reference <- transformed_centroids[,1]
# df_masked$GM130_centroidy_reference <- transformed_centroids[,2]
# 
# transformed_centroids <- transform_reference_centroids_xy(df_masked %>% select(old_filename_generic_noside, intersectionx, intersectiony), info_df)
# df_masked$reference_intersectionx <- transformed_centroids[,1]
# df_masked$reference_intersectiony <- transformed_centroids[,2]
#   
# # create grid of points for TPS warp
# # make TPS
# computeTransform(data.matrix(fixed_lms), data.matrix(samples_lms), type='tps')
# # apply warp
# applyTransform(data.matrix(transformed_centroids[,1:2]), test_transform)
# 
# 
# 
# # adjust the angles 0-2pi and flip the treated side to match (it's opposite of above since they are already flipped from overview)
# df_masked <- df_masked %>% mutate(positional_angle = case_when(positional_angle > 2*pi ~ positional_angle - 2*pi,
#                                                                                  positional_angle < 0 ~ angle + 2*pi,
#                                                                                  TRUE ~ positional_angle))
# 
# df_masked$positional_angle_old <- df_masked$positional_angle
# df_masked <- df_masked %>% mutate(positional_angle = case_when(side == 'treated' & positional_angle <= pi ~ pi-positional_angle, side == 'treated' & positional_angle > pi ~ 3*pi - positional_angle, TRUE ~ positional_angle_old))
# 
# saveRDS(df_masked, file='combined_angles_positional_masked.Rda') 


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


#### 0.3.1 Circular statistics helpers ####
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

#### 0.3.2 BPNME statistics helpers ####
get_posterior_estimates <- function(fitted) {
  # this code is from Cremers, unfortunately it will not run  in a reasonable 
  # time for this many samples (cells)
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
      scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0)) + #,  breaks = c(0,.01,.02,.03,.04)) +
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
      scale_y_continuous(limits = c(0, 0.097), expand = c(0, 0),  breaks = c(0,.025,.05,.075,.1)) +
      coord_polar(start = ((270-(dirres/2))) * pi/180, direction = -1) +
      scale_fill_manual(name = "treatment", values = color, drop = FALSE) +
      scale_color_manual(name = "treatment", values = c('black','black'), drop = FALSE) +
      scale_alpha_discrete(range = c(1, .5)) +
      theme_bw() +
      theme(axis.title.x = element_blank(), legend.position="none")
  }
  return(p.windrose)
}


positional_graph_sliding_mean <- function(temp_df, perim_n, window) {
  angle_window <- window/2*pi/180 # convert window to radians and divide by 2
  moving_avg_pts <- perim_n / 360 * window #number of perimeter pts in each total window (for calculating mean)
  mean_df <- data.frame(matrix(nrow = perim_n))
  for (angle_i in 1:perim_n) {
    angle <- angle_i * pi/(perim_n/2)
    mean_df[angle_i,1] <- angle
    if (angle - angle_window < 0) {
      # I have no idea why this is necessary, seems useless but it doesn't work otherwise
      a <- angle + 2*pi - angle_window
      b <- angle+angle_window
      mean_df[angle_i,2] <- (temp_df %>% filter(between(positional_angle, a, 2*pi) | between(positional_angle, 0, b)) %>% count()) / moving_avg_pts
    } else if (angle + angle_window > 2*pi) {
      a <- angle - angle_window
      b <- angle + angle_window - 2*pi
      mean_df[angle_i,2] <- (temp_df %>% filter(between(positional_angle, a, 2*pi) | between(positional_angle, 0, b)) %>% count()) / moving_avg_pts
    } else {
      a <- angle - angle_window
      b <- angle + angle_window
      mean_df[angle_i,2] <- (temp_df %>% filter(between(positional_angle, a, b)) %>% count()) / moving_avg_pts
    }
    mean_df[angle_i,3] <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_x, xout = angle)$y*-1
    mean_df[angle_i,4] <- approx(overview_pos_LUT$pos_angle, overview_pos_LUT$overview_y, xout = angle)$y*-1
  }
  mean_df[,2] <- mean_df[,2] * 100 / nrow(temp_df) 
  return(mean_df)
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


#### Main ####
#### 1 Setting up the DF ####
# R doesn't have a good way to get the path of where this file is located, unfortunately
# if you are running this code in Rstudio, try this:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("../../data/GM130")
# if at work...
# setwd("C:/Users/nhanne/Box/FGF_inhibitor_paper_5-26-2020/data/GM130")
getwd() #check our working directory
# if you aren't using rstudio, use the command setwd() 
# and point it to the data/GM130 directory


### Prepare the cellpose assisted matching angles
# This step can take a long time! More than a minute
# if you've run it before then just load the combined data here
df <- readRDS("combined_angles.Rda")
df_masked <- readRDS('combined_angles_positional_masked.Rda') 
info_df <- read.csv('GM130_image_log.csv')
# else you generate them again here and then
# save the combined dataset for faster loading next time
# df <- compile_angle_df() # this will take a while! be patient
# saveRDS(df, file='combined_angles.Rda') 


### Try to align landmarks between overview image and each side stack
# Many of them have been flipped already from the overview, but it seems like the 
# right side of overview is the bead treated side. Let's do landmark 1-2 be contra
# side and landmark 3-4 be treated side
# need to perform this before the following angle adjustment and flipping as we want them to 
# be aligned with the overview image!
landmark_filenames <- list.files(path="./imagej_rois/overview_landmarks/", pattern=".zip$", recursive=TRUE)
landmark_filenames <- data.frame(landmark_filenames)

# this will save all the transform matrices as lists in the info_df
info_df$old_filename_generic <- str_remove_all(info_df$old_filename, "_nuc.tif|_gm130.tif")
info_df$old_filename_generic_noside <- str_remove_all(info_df$old_filename_generic, "_treated|_control")
info_df <- get_transform_matrices_from_rois(info_df, landmark_filenames)


df <- left_join(df, info_df %>% select(new_filename, old_filename_generic_noside) %>% distinct(), by=c('t' = 'new_filename'))

# this will apply them to the actual xy data in the main df
# df$nuclei_centroidx_overview <- NA
# df$nuclei_centroidy_overview <- NA
# df$GM130_centroidx_overview <- NA
# df$GM130_centroidy_overview <- NA
# 
# transformed_centroids <- transform_centroids_xy(df %>% select(t, nuclei_centroidx, nuclei_centroidy), info_df)
# df$nuclei_centroidx_overview <- transformed_centroids[,1]
# df$nuclei_centroidy_overview <- transformed_centroids[,2]
# 
# transformed_centroids <- transform_centroids_xy(df %>% select(t, GM130_centroidx, GM130_centroidy), info_df)
# df$GM130_centroidx_overview <- transformed_centroids[,1]
# df$GM130_centroidy_overview <- transformed_centroids[,2]


### Correct angles on any images that are not 'square'
# rotates all angles in a given image by a set angle
# in this data we have a horizontal angle pointing right as 0 degrees
# CCW rotation is positive (right hand rule)
# I draw a line in in imagej from the left to right along the 'flat' of the FNP, 
# then click 'measure' (ctrl+m) and get the angle
# just subtract this angle to all other angles
# see the example image in the github readme
# df <- adjust_base_angles(df, info_df)
# 
# 
# ### Flip all the control angles across the 'y' axis so they match w/ treatment side
# # NOTE: For some reason the LY group ones are reversed? Probably how the section was put on the slide...
# # this is very important part of the analysis
# df <- flip_y_angles(df)


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
# df_masked <- mask_out_centroids(df, mask_filenames, FALSE)


### filter based on distance from the 'baseline' roi line
# this is very similar to the last chuck of code above
# Here I made line roi's in ImageJ along the top of the nasalpit, globular process, and lateral FNP
# I'll put another example in the github readme
# We will only look at nuclei-Golgi pairs that are within 200 um of this line
# Basically we don't want to analyze pairs that are more in the mid-FNP or too 
# close to the neural ectoderm
# I think I want to get rid of this for now...
# baseline_mask_filenames <- list.files(path="./imagej_rois/baseline_rois/", pattern=".roi$")
# this is another slow one, be patient
# df_baseline_masked <- filter_baseline_distance(baseline_mask_filenames, df_masked, 200, FALSE)


#### 2 Region 'aiming' analysis ####
### Instead of using the angle b/w the nuclei and Golgi we will make a pseudo-angle 
### based on what boundary part of the tissue the angle is pointing at
# This means that a cell in the edge of globular process pointing 'left' at base of
# gp-nasal pit boundary will have same angle as a cell at the top of the nasal pit
# pointing 'down' at tje gp-np boundary. They are pointing at the same feature
# despite having true angles nearly 90 degrees apart.

# To perform this analysis I have drawn eight line rois around the boundary of tissue
# on the overview images that (mostly) capture the entire FNP and nasal pit area.
# Each octile will encompass (ha!) 45 degrees on a windrose plot. The actual angle
# will be calculated as the distance along the octile roi to the intersection 
# of a line radiating (more dumb word play) from the nuclei centroid towards its'
# golgi centroid.


#### 2.1 calculating and transforming data ####
# load in the rois. There are 8, starting at the left (contralateral) nasal pit and moving ccw
# this is counter-intuitive as I had already rotated all the stacks so that treated would be on the left side

# 1 top of left np to gp
# 2 gp to midline
# 3 midline to right gp
# 4 gp to np
# 5 np to neural area lining up in y with gp octile break
# 6 neural area to midline y octile break
# 7 neural midline to left nueral gp break
# 8 neural gp break to top of left np
octile_filenames <- list.files(path="./imagej_rois/overview_octiles/", pattern=".zip$", recursive=TRUE)
octile_filenames <- data.frame(octile_filenames)

# this is quite slow, gets the positional angles
# df_masked$positional_angle <- NA
# df_masked$intersectionx <- NA
# df_masked$intersectiony <- NA
# df_masked <- get_positional_angle(df_masked, octile_filenames)
# df_masked$positional_angle <- as.numeric(df_masked$positional_angle)

# This puts the positional angles onto the main 'overview' cartoon
overview_octile_rois <- octile_filenames %>% filter(octile_filenames == 'overview_octile.zip')
overview_pos_LUT <-  readRDS(file='smooth_LUT.Rda')
# overview_pos_LUT <- generate_overview_positional_LUT(overview_octile_rois)
# saveRDS(overview_pos_LUT, file='smooth_LUT.Rda') 
df_masked <- convert_directional_angle_overview_LUT(df_masked, overview_pos_LUT)
df_masked <- df_masked %>% na.omit(positional_angle)

tmp_list <- assign_tissue_region(df_masked, info_df, octile_filenames)
df_masked <- tmp_list[[1]]
info_df <- tmp_list[[2]]

# need to transform the centroids again here to the overview octile image
# 1 pick two positional angles and use them to re-calculate rotation matrices
# this isn't great, but it gets us closer at least
# the warp will get the rest of the way
# info_df <- get_overview_transform_matrices_from_pos_angles(info_df, df_masked, 100) # not sure this is necessary anymore, may remove
df_masked <- morph_centroids_to_ref_img(df_masked, overview_pos_LUT, overview_octile_rois, info_df)

df_masked$treatment <- as.factor(df_masked$treatment)
df_masked$side <- as.factor(df_masked$side)
df_masked$region_name <- as.factor(df_masked$region_name)
df_baseline_masked <- filter_distance_from_octile_lines(df_masked, octile_filenames, info_df, 0,200)


#### 2.2 overview heatmap plot ####
for (treatment_i in 1:length(levels(df_baseline_masked$treatment))) {
  for (side_i in 1:2) {
    for (region_i in 1:length(levels(df_baseline_masked$region_name))) {
      group_mask <- df_baseline_masked %>% filter(as.integer(treatment) == treatment_i)
      
      # combine centers, there aren't many points here
      if (levels(df_masked$region_name)[region_i] == 'center') {
        filter_mask <- group_mask %>% filter(as.integer(region_name) == region_i)
      } else {
        filter_mask <- group_mask %>% filter(as.integer(region_name) == region_i & as.integer(side) == side_i)
      }

      n_samples <- length(unique(filter_mask$old_filename_generic_noside))
      for (sample_i in 1:n_samples) {
        sample <- unique(filter_mask$old_filename_generic_noside)[sample_i]
        mean_df <- positional_graph_sliding_mean(filter_mask %>% filter(old_filename_generic_noside == sample), 2000, 10)
        mean_df <- mean_df %>% mutate(sample =  sample)
        if (sample_i == 1) {
          collected_means_df <- mean_df
        } else {
          collected_means_df <- rbind(collected_means_df, mean_df)
        }
      }
      
      sum_collected_means_df <- collected_means_df %>% group_by_at(c(1,3,4)) %>% dplyr::summarize(group_mean = mean(n))
      test_pts <- st_as_sf(sum_collected_means_df, coords = c("V3","V4"), remove = FALSE)
      # cat(levels(df_masked$treatment)[treatment_i], ' ', levels(df_masked$side)[side_i], ' ', levels(df_masked$region_name)[region_i])
      # cat('min: ', min(test_pts$group_mean))
      # cat('max: ', max(test_pts$group_mean))
      gradient_max <- 0.3
      if (filter_mask[1,]$region_name == 'center') {
        gradient_max <- 0.5
      }
      
      # this graph is going to be 'upside down' in order to look normal, but it therefore
      # will have reflected the left and right sides. Can be fixed in illustrator
      
      aim_plot <- ggplot(test_pts) + 
                  #geom_point(data=group_mask, aes(x = morphedx*-1, y = morphedy*-1), stroke=0, shape=21, color=alpha('black', .2)) +
                  #geom_point(data=filter_mask, aes(x = morphedx*-1, y = morphedy*-1), color='red') +
                  geom_sf(aes(color = group_mean, size = group_mean)) +
                  scale_colour_viridis_c(values = rescale(c(seq(0,.1,.002), seq(.1,.2,.01), seq(.2,.3,.01))), limits=c(0,gradient_max), option = 'inferno', alpha=.4, breaks = c(0, 0.05, 0.1,0.2,0.3))  +
                  # scale_fill_gradient(limits=c(0,gradient_max)) 
                  scale_size_continuous(limits=c(0,gradient_max), range=c(10,35), breaks = c(0.05, 0.1,0.2,0.3))
                  # scale_fill_viridis_opt(test_pts$group_mean)

      
      file_name <- paste0("./figs/positional_map_0-200_", as.character(levels(df_baseline_masked$treatment)[treatment_i]), '_', as.character(levels(df_baseline_masked$side)[side_i]), '_', as.character(levels(df_baseline_masked$region_name)[region_i]), sep="")
      # ggsave(filename=paste0(file_name, '.tiff'), aim_plot, width = 30, heigh = 30, units='cm')
      
      pdf(paste0(file_name, ".pdf"), width=15, height=15)
      plot(aim_plot)
      dev.off()
        #
    }
  }
}


#### 3 Quadrant comparison and plotting ####
# The windrose plots look good, but it would be nice to get a more quantitative
# analysis of group differences. I will get relative amount in each quadrant
# and compare to contralateral and DMSO groups

df_baseline_masked <- df_baseline_masked  %>% mutate(treatment = factor(treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', '3Mix')))
df_baseline_masked$flipped_positional_angle <- df_baseline_masked$positional_angle
df_baseline_masked <- df_baseline_masked  %>% mutate(flipped_positional_angle = case_when((side == 'control' & positional_angle <= pi) ~ pi-positional_angle, (side == 'control' & positional_angle > pi) ~ 3*pi - positional_angle, TRUE ~ flipped_positional_angle))
df_baseline_masked <- df_baseline_masked  %>% mutate(quad_bin = cut(flipped_positional_angle, breaks = c(0, pi/2, pi, 3*pi/2, 2*pi)))

for (region_i in 1:(length(levels(df_baseline_masked$region_name))+1)) {
  if (region_i == length(levels(df_baseline_masked$region_name))+1) {
  bin_summary <- df_baseline_masked %>% group_by(sample_info, treatment, side) %>% drop_na(quad_bin) %>% count(quad_bin) %>% mutate(freq = n / sum(n))
  region_name <- 'all'
  } else {
  bin_summary <- df_baseline_masked %>% filter(as.integer(region_name) == region_i) %>% group_by(sample_info, treatment, side) %>% drop_na(quad_bin) %>% count(quad_bin) %>% mutate(freq = n / sum(n))
  region_name <- as.character(levels(df_baseline_masked$region_name)[region_i])
  }

  group_bin_summary <- bin_summary %>% group_by(treatment, side, quad_bin) %>% dplyr::summarize(mean_freq = mean(freq))
  
  for (quadrant in 1:length(levels(bin_summary$quad_bin))) {
    group_bin_summary_quad <- group_bin_summary %>% filter(quad_bin == levels(group_bin_summary$quad_bin)[quadrant])
    bin_summary_quad <- bin_summary %>% filter(quad_bin == levels(group_bin_summary$quad_bin)[quadrant])
    
    graph_mean_df <- bin_summary_quad %>% group_by(treatment, side) %>%  summarise(
      y0 = quantile(freq, 0.05),
      y25 = quantile(freq, 0.25),
      y50 = mean(freq),
      y75 = quantile(freq, 0.75),
      y100 = quantile(freq, 0.95),
      ysd = sd(freq))
    
    # p <- ggplot(data=graph_mean_df, aes(treatment, y50, fill = as.factor(side))) +
    #   geom_col(stat = "identity", position = 'dodge') +
    #   geom_errorbar(aes(ymin=y50, ymax=y50+ysd), position = 'dodge', width = 1) +
    #   geom_point(data=bin_summary_quad, aes(treatment, freq, fill=side), size = 1, position=position_jitterdodge(jitter.width = 0.2)) +
    #   theme(legend.position = "none") # + ylim(0,0.6)
    file_name <- paste0("./figs/positional_avg_", region_name, '_', as.character(levels(group_bin_summary$quad_bin)[quadrant]), ".tiff", sep="")
    # ggsave(filename=file_name, p, width = 25, heigh = 15, units='cm')
    # print(p)
    
    
    # pdf(paste0("./figs/positional_avg", as.character(levels(group_bin_summary$quad_bin)[quadrant]), "_contralateral.pdf"), width=10, height=6)
    # p <- ggplot() + geom_bar(data = group_bin_summary_quad, aes(y=mean_freq, x = side), stat="identity") +
    #   geom_jitter(data = bin_summary_quad, aes(x = side, y = freq), shape=16) +
    #   # geom_errorbar(data = cellpose_summarise, aes(y=mean,x=side,ymin=mean-sd,ymax=mean+sd)) +
    #   facet_wrap(~ treatment)
    # print(p)
    # dev.off()
    
    test<- aov(bin_summary_quad$freq ~  bin_summary_quad$treatment * bin_summary_quad$side)
    print(file_name)
    print(summary(test))
    print(TukeyHSD(test))
  }
}


#### 4 Windrose plotting ####
graphing_df <- df_baseline_masked #%>% filter(region_name != 'center')

# it's possible the angles weren't flipped, can be done here just in case
graphing_df <- flip_y_angles(graphing_df)

# graphing_df$angle_deg <- graphing_df$angle * 180 / pi
graphing_df$angle_deg <- as.numeric(graphing_df$angle) * (180 / pi)
graphing_df <- graphing_df %>% drop_na(positional_angle)
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
  
  # # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"), units='in', width=5, height=5, res=300)
  # pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_contralateral.pdf"), width=5, height=5)
  # control_plot <- plot.windrose(filter_data %>% filter(side == 'contralateral'), 'white', dirres = 10)
  # plot(control_plot)
  # dev.off()
  # # dev.off()
  # 
  # # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"), units='in', width=5, height=5, res=300)
  # pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated.pdf"), width=5, height=5)
  # treated_plot <- plot.windrose(filter_data %>% filter(side == 'treated'), 'red', dirres = 10)
  # plot(treated_plot)
  # dev.off()
  # # dev.off()
  
  if (i == 1) {
    pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "combined.pdf"), width=5, height=5)
    combined_plot <- plot.windrose(filter_data, 'white', dirres = 10)
    plot(combined_plot)
    dev.off()
  }
}

## temp code
for (i in 1:length(unique(filter_data$old_filename_generic_noside.x))) {
  i = 6
  further_filt <- filter_data %>% filter(old_filename_generic_noside.x == unique(filter_data$old_filename_generic_noside.x)[i])
  plot(plot.windrose(further_filt %>% filter(side == 'contralateral'), 'white', dirres = 10))
  plot(plot.windrose(further_filt %>% filter(side == 'treated'), 'red', dirres = 10))
}



graphing_positional_df <- df_baseline_masked #%>% filter(region_name != 'center')

# it's possible the angles weren't flipped, can be done here just in case
graphing_positional_df$flipped_positional_angle <- graphing_positional_df$positional_angle
graphing_positional_df <- graphing_positional_df  %>% mutate(flipped_positional_angle = case_when((side == 'control' & positional_angle <= pi) ~ pi-positional_angle, (side == 'control' & positional_angle > pi) ~ 3*pi - positional_angle, TRUE ~ flipped_positional_angle))


# graphing_df$angle_deg <- graphing_df$angle * 180 / pi
graphing_positional_df$angle_deg <- as.numeric(graphing_positional_df$positional_angle) * (180 / pi)
graphing_positional_df <- graphing_positional_df %>% drop_na(positional_angle)
graphing_positional_df$rel_z <- graphing_positional_df$delta_z / graphing_positional_df$distance
graphing_positional_df <- graphing_positional_df %>% mutate(side_spd = case_when(side == 'control' ~ 0,
                                                           side == 'treated' ~ 1))
graphing_positional_df <- graphing_positional_df %>% unite(sample_side, sample_info, side, sep='_', remove=FALSE)

graphing_positional_df$treatment <- factor(graphing_positional_df$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', '3Mix'))
graphing_positional_df <- graphing_positional_df %>% mutate(treatment = recode(treatment, '3Mix' = 'Triple'))
graphing_positional_df$side <- factor(graphing_positional_df$side, levels = c('control', 'treated'))
graphing_positional_df <- graphing_positional_df %>% mutate(side = recode(side, 'control' = 'contralateral'))

for (i in 1:length(levels(graphing_positional_df$treatment))) {
  filter_data <- graphing_positional_df %>% filter(as.integer(treatment) == i)
  
  # # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_control_vsDMSO.png"), units='in', width=5, height=5, res=300)
  # pdf(paste0("./figs/positional_windrose_", as.character(filter_data$treatment[1]), "_contralateral.pdf"), width=5, height=5)
  # control_plot <- plot.windrose(filter_data %>% filter(side == 'contralateral'), 'white', dirres = 10)
  # plot(control_plot)
  # dev.off()
  # # dev.off()
  # 
  # # png(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "_treated_vsDMSO.png"), units='in', width=5, height=5, res=300)
  # pdf(paste0("./figs/positional_windrose_", as.character(filter_data$treatment[1]), "_treated.pdf"), width=5, height=5)
  # treated_plot <- plot.windrose(filter_data %>% filter(side == 'treated'), 'red', dirres = 10)
  # plot(treated_plot)
  # dev.off()
  # # dev.off()
  
  if (i == 1) {
    filter_data$angle_deg <- as.numeric(filter_data$flipped_positional_angle) * (180 / pi)
    pdf(paste0("./figs/windrose_", as.character(filter_data$treatment[1]), "combined.pdf"), width=5, height=5)
    combined_plot <- plot.windrose(filter_data, 'white', dirres = 10)
    plot(combined_plot)
    dev.off()
  }
}


#### 5 Watson U2 tests ####
### First do some traditional summary statistics
circular_statistics <- list()
# get actual circle stats for export
df_watson <- df_baseline_masked #%>% filter(region_name == 'mid')

# it's possible the angles weren't flipped, can be done here just in case
df_watson <- flip_y_angles(df_watson)

for (i in 1:length(levels(as.factor(df_watson$treatment)))) {
  for (j in 1:length(levels(as.factor(df_watson$side)))) {
    temp_circ_stats <- get_circular_stats(levels(as.factor(df_watson$treatment))[i], levels(as.factor(df_watson$side))[j], df_watson)
    circular_statistics <- append(circular_statistics, list(temp_circ_stats))
  }
}

circular_statistics <- do.call(rbind.data.frame, circular_statistics)
colnames(circular_statistics) <- c('treatment', 'side', 'ind_mean', 'overall_mean', 'ind_sd', 'mean_sd', 'overall_sd', 'ind_cr', 'mean_cr', 'overall_cr')
circular_statistics <- circular_statistics %>% mutate(across(ind_mean:overall_sd, ~.x*180/pi)) # convert to degrees
circular_statistics <- circular_statistics %>% mutate(across(ind_mean:overall_sd, ~case_when(. < 0 ~ . + 360, TRUE ~ .))) # remove negative

### compare side v side for each treatment
watson_side_results <- list()
for (i in 1:length(levels(as.factor(df_watson$treatment)))) {
  watson <- watson.two.test(df_watson %>% filter(treatment == levels(as.factor(df_watson$treatment))[i], side =='control') %>% select(angle), 
                            df_watson %>% filter(treatment == levels(as.factor(df_watson$treatment))[i], side == 'treated') %>% select(angle))
  # print(watson)
  watson_result_list <- list(levels(as.factor(df_watson$treatment))[i], watson[[1]])
  watson_side_results <- append(watson_side_results, list(watson_result_list))
}

### compare each side-treatment combo against DMSO combined
watson_results <- list()
for (i in 1:length(levels(as.factor(df_watson$treatment)))) {
  for (j in 1:length(levels(as.factor(df_watson$side)))) {
    watson <- watson.two.test(df_watson %>% filter(treatment == 'DMSO') %>% select(angle),
                              df_watson %>% filter(treatment == levels(as.factor(df_watson$treatment))[i], side == levels(as.factor(df_watson$side))[j]) %>% select(angle))
    # print(watson)
    watson_result_list <- list(levels(as.factor(df_watson$treatment))[i], levels(as.factor(df_watson$side))[j], watson[[1]])
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
write.csv(combined_results, 'Golgi_analysis_output_all.csv')

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


#### 6 Mollweide Plots ####
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


#### 7 BPNME statistics ####
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


#### 8 Cellularity #### 
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
