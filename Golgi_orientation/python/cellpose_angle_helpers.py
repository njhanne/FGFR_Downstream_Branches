import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import cv2
from skimage import measure
import tifffile

from scipy.spatial import distance
from scipy.optimize import linear_sum_assignment
import math
from joblib import Parallel, delayed
import time


from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths

# Image Processing and Geometry Helpers
def get_all_centroids(image, dimension, img_path, centroid_paths, rerun, anis=1):
  save_name = img_path[0].stem + '_centroids.npy'
  centroid_path = [cp for cp in centroid_paths if cp.parts[-1].endswith(save_name)]
  if len(centroid_path) != 0 and not rerun: # don't rerun if there is an existing file, unless we manually force it
    if centroid_path[0].lstat().st_mtime > img_path[0].lstat().st_mtime:
      # if the centroid file is newer than the last time the cellpose image was generated then load it
      centroids = np.load(centroid_path[0])
      print("Loading previous centroids...")
    else:
      # if the centroid file is older than the new cellpose image then regenerate it
      centroids = get_centroids(image, dimension, anis)
      np.save(str(img_path[0].parent / save_name), centroids)
  else:
    centroids = get_centroids(image, dimension, anis)
    np.save(str(img_path[0].parent / save_name), centroids)
  return centroids


def get_centroids(image, dimension, anis=1):
  centroids = []
  max = np.amax(image)+1

  if dimension == 2:
    # TODO add bounding box for 2D
    print('2D centroids code needs simple re-write')
    # centroids.append(calculate_centroid(np.where(image == i, 1, 0)))
  if dimension == 3:
    # first get a bounding box around desired cell, this is way faster than using entire image
    t1 = time.time()
    bbox = Parallel(n_jobs=4)(delayed(bbox_3D)(image, i) for i in range(1, max))
    bbox = np.array(bbox)
    t2 = time.time()
    print(str(max) + ' bounding boxes in ' + str(int(t2-t1)) + ' seconds')
    print((t2-t1) / max)

    # calculate the centroids
    centroids = Parallel(n_jobs=4)(delayed(calculate_3D_centroid)(bbox, image, i, anis) for i in range(0, max-1))
    centroids = np.array(centroids)
    centroids2 = np.transpose(np.array([centroids[:,0] + bbox[:,5],  centroids[:,1] + bbox[:,3],  centroids[:,2] + bbox[:,1]*anis, centroids[:,3]]))
  return centroids2


def bbox_3D(image, i):
  # https://stackoverflow.com/a/31402351
  # this is at least 3x faster than not using a bounding box
  # look into this https://forum.image.sc/t/making-np-where-faster/58641/7
  img = np.where(image == i, 1, 0) # binarize image to only select desired cell
  z = np.any(img, axis=(1, 2)) # return only non-zero axes
  y = np.any(img, axis=(0, 2))
  x = np.any(img, axis=(0, 1))

  xmin, xmax = np.argmax(x), x.size - 1 - np.argmax(x[::-1])
  ymin, ymax = np.argmax(y), y.size - 1 - np.argmax(y[::-1])
  zmin, zmax = np.argmax(z), z.size - 1 - np.argmax(z[::-1])
  bbox = [i, zmin,zmax, ymin,ymax, xmin,xmax]
  return bbox


def calculate_centroid(image):
  # pass in the mask from cellpose .tif and it will output the xy coordinates of the centroid
  M = cv2.moments(image.astype('uint8')) # moments needs uint8 for some reason...
  centroid_X = M["m10"] / M["m00"]
  centroid_Y = M["m01"] / M["m00"]
  return centroid_X, centroid_Y


def calculate_3D_centroid(bboxes, image, i, anis = 1):
  # openCV doesn't handle more than 2D so we have to use skimage
  # pass in the mask from cellpose .tif and it will output the xyz coordinates of the centroid
  # skimage does z,y,x
  img = image[bboxes[i,1]:bboxes[i,2], bboxes[i,3]:bboxes[i,4], bboxes[i,5]:bboxes[i,6]]
  M = measure.moments(np.where(img==i+1, 1, 0).astype('uint8'), order=1, spacing=(anis,1,1)) # moments needs uint8 for some reason...
  centroid_X = M[0,0,1] / M[0,0,0]
  centroid_Y = M[0,1,0] / M[0,0,0]
  centroid_Z = M[1,0,0] / M[0,0,0]
  # 0th moment is the volume, volume of sphere is 4/3 * pi * r^3
  radius = ((3*M[0,0,0])/(4*math.pi))**(1/3)
  return [centroid_X, centroid_Y, centroid_Z, radius]


def build_distance_cost_matrix(nuc_centroids, golgi_centroids, nuc_radii, smallest_nuc, matches = None):
  # make a mask same size as the list of centroids
  # nuc_centroid_cull = np.ones(len(nuc_centroids[:, 0]))
  golgi_centroid_cull = np.ones(len(golgi_centroids[:, 0]))

  # if there are already matches, filter them out of the mask so they aren't included in matrix
  # TODO: rewrite or remove - deprecated
  # if matches is not None:
    # matches are in image order... if there are x nuclei, the centroid index would be x-1
    # nuc_centroid_cull[matches[:,0] - 1] = 0
    # golgi_centroid_cull[matches[:, 1] - 1] = 0

  # filter out very small nuclei (probably bad cellpose segmentation)
  if smallest_nuc is not None:
    nuc_centroid_cull_first = np.ones(len(nuc_centroids[:, 0]))
    nuc_centroid_cull_first[nuc_radii < smallest_nuc] = 0
    nuc_centroid_cull_first = np.where(nuc_centroid_cull_first)[0]

  # filter out centroids that didn't calculate correctly
  nuc_centroid_cull_second = np.arange(0,len(nuc_centroids),1)
  nuc_centroid_cull_second = (np.delete(nuc_centroid_cull_second, np.searchsorted(nuc_centroid_cull_second, np.unique(np.where(np.isnan(nuc_centroids))[0])))) # removes nans

  nuc_centroid_cull = np.intersect1d(nuc_centroid_cull_first, nuc_centroid_cull_second)

  golgi_centroid_cull = np.where(golgi_centroid_cull)
  golgi_centroid_cull = (np.delete(golgi_centroid_cull[0], np.searchsorted(golgi_centroid_cull[0], np.unique(np.where(np.isnan(golgi_centroids))[0])))) # removes nans

  # calculate distance between centroids
  cost_matrix = distance.cdist(nuc_centroids[nuc_centroid_cull], golgi_centroids[golgi_centroid_cull])
  return cost_matrix, nuc_centroid_cull, golgi_centroid_cull


def otsu_test(image):
  # This usually is used with images, but this code works on any data
  # Taken from https://learnopencv.com/otsu-thresholding-with-opencv/
  # Set total number of bins in the histogram
  bins_num = 256

  # Get the image histogram
  hist, bin_edges = np.histogram(image, bins=bins_num)

  # Calculate centers of bins
  bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2.

  # Iterate over all thresholds (indices) and get the probabilities w1(t), w2(t)
  weight1 = np.cumsum(hist)
  weight2 = np.cumsum(hist[::-1])[::-1]

  # Get the class means mu0(t)
  mean1 = np.cumsum(hist * bin_mids) / weight1
  # Get the class means mu1(t)
  mean2 = (np.cumsum((hist * bin_mids)[::-1]) / weight2[::-1])[::-1]

  inter_class_variance = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

  # Maximize the inter_class_variance function val
  index_of_max_val = np.argmax(inter_class_variance)

  threshold = bin_mids[:-1][index_of_max_val]
  print("Otsu's algorithm implementation thresholding result: ", threshold)
  return threshold

def get_angles(matches, nuc_centroids, golgi_centroids, dimension):
  # angles are ccw starting at 0 on the positive x-axis
  angles = []
  delta_z = []
  for pair in matches:
    if dimension == 2:
      x0, y0 = nuc_centroids[pair[0]-1]
      x1, y1 = golgi_centroids[pair[1]-1]
      dz = 0
    else:
      x0, y0, z0 = nuc_centroids[pair[0]-1]
      x1, y1, z1 = golgi_centroids[pair[1]-1]
      dz = z1 - z0
    angle = -math.atan2(y1 -y0, x1-x0)
    if angle < 0:
      angle = angle + 2*math.pi
    if angle > 2*math.pi:
      angle = angle - 2*math.pi
    angles.append(angle)
    delta_z.append(dz)
  return np.array(angles), delta_z


def get_unit_vector(matches, nuc_centroids, golgi_centroids, dimension):
  # get a unit vector with 0,0,0 as nucleus centroid and terminating at golgi centroid
  x_hat = []
  y_hat = []
  z_hat = []
  for pair in matches:
    if dimension == 2:
      x0, y0 = nuc_centroids[pair[0]-1]
      x1, y1 = golgi_centroids[pair[1]-1]
      dx = x1 - x0
      dy = y1 - y0
      dz = 0
    else:
      x0, y0, z0 = nuc_centroids[pair[0]-1]
      x1, y1, z1 = golgi_centroids[pair[1]-1]
      dx = x1 - x0
      dy = y1 - y0
      dz = z1 - z0
    vlen = math.sqrt(dx**2 + dy**2 + dz**2)
    x_hat.append(dx/vlen)
    y_hat.append(dy/vlen)
    z_hat.append(dz/vlen)
  return x_hat, y_hat, z_hat


# End Image Process and Geometry Helpers

# Plotting Helpers
def plot_all_thresh(nuc_img, golgi_img, pairs, nuc_centroids, golgi_centroids, distances, otsu, dimension, anis=1, save_path=None):
  if dimension == 2:
    nuc = np.where(nuc_img != 0, 1, 0)
    golgi = np.where(golgi_img != 0, 2, 0)
    image = nuc + golgi
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(image)
    for ind, pair in enumerate(pairs):
      x0, y0 = nuc_centroids[pair[0] - 1]
      x1, y1 = golgi_centroids[pair[1] - 1]
      if distances[ind] > otsu:
        ax.plot((x0, x1), (y0, y1), '-r', linewidth=2)
      else:
        ax.plot((x0, x1), (y0, y1), '-w', linewidth=1)
    ax.axis((0, 2048, 2048, 0))
  else:
    fig, ax = plt.subplots(figsize=(20, 20))
    xpix = nuc_img.shape[1]
    ypix = nuc_img.shape[2]
    ax.axis((0, xpix, ypix, 0))
    nuc_masks = np.zeros((xpix, ypix))
    golgi_masks = np.zeros((xpix, ypix))
    for ind, pair in enumerate(pairs):
      # if (pair[0] != 1 & pair[1] != 1):
      if (pair[0] != nuc_centroids.shape[0]) & (pair[1] != golgi_centroids.shape[0]):
        if not np.isnan(nuc_centroids[pair[0]-1][0]):
          nuc_mask = np.where(nuc_img[int(nuc_centroids[pair[0]-1][2] /anis),:,:] == (pair[0]), 1, 0)
          golgi_mask = np.where(golgi_img[int(golgi_centroids[pair[1]-1][2] / anis),:,:] == (pair[1]), 1, 0)
          nuc_masks = nuc_masks.astype(int) | nuc_mask.astype(int)
          golgi_masks = golgi_masks.astype(int) | golgi_mask.astype(int)
          # test_helper(nuc_mask, golgi_mask, nuc_centroids, golgi_centroids, pair)
          # print('done')
    image = nuc_masks + golgi_masks*2
    ax.imshow(image)
    for ind, pair in enumerate(pairs):
      x0, y0, z0 = nuc_centroids[pair[0] - 1]
      x1, y1, z0 = golgi_centroids[pair[1] - 1]
      # if distances[ind] > otsu:
      #   # ax.plot((x0, x1), (y0, y1), '-r', linewidth=2)
      # else:
      ax.plot((x0, x1), (y0, y1), '-w', linewidth=1)
  if save_path is not None:
    plt.savefig(save_path)
  plt.show()
  plt.close()


### MAIN ###
# First find all needed directories and load them all up
data_dir = (Path.cwd().parent.parent / 'data' / 'GM130').resolve()

info_csv_path = Path(data_dir / 'GM130_image_log.csv')
sample_info = pd.read_csv(info_csv_path)

# get all the stacks we created from rgb_to_stack.py
nuclei_directory = (data_dir / 'cellpose' / 'nuclei').resolve()
golgi_directory = (data_dir / 'cellpose' / 'golgi').resolve()
results_directory = (data_dir / 'cellpose' / 'angle_results').resolve()

# where to save everything
output_directory = (data_dir / 'cellpose' / 'analysis_output').resolve()

nuclei_dirs, nuclei_mask_paths = find_all_filepaths(nuclei_directory, '.tif')
nuclei_mask_paths = [mp for mp in nuclei_mask_paths if mp.parts[-1].endswith('_masks.tif')]
golgi_dirs, golgi_mask_paths = find_all_filepaths(golgi_directory, '.tif')
golgi_mask_paths = [mp for mp in golgi_mask_paths if mp.parts[-1].endswith('_masks.tif')]

# this loads up the centroids if we've already run the analysis before. It's WAY faster than re-running the whole thing
nuc_centroid_dirs, nuclei_centroid_paths = find_all_filepaths(nuclei_directory, '.npy')
nuclei_centroid_paths = [cp for cp in nuclei_centroid_paths if cp.parts[-1].endswith('_centroids.npy')]
golgi_centroid_dirs, golgi_centroid_paths = find_all_filepaths(golgi_directory, '.npy')
golgi_centroid_paths = [cp for cp in golgi_centroid_paths if cp.parts[-1].endswith('_centroids.npy')]

# this loads up previous results csv, with an option to rerun the analysis even if it's done, like if you've made changes or something
rerun_pair_finding = True
rerun_centroid_finding = True
results_dir, results_paths = find_all_filepaths(results_directory, '.csv')
results_paths = [rp for rp in results_paths if rp.parts[-1].endswith('_angle_results.csv')]


cols = ['sample', 'side', 'section', 'view']
sample_info['sample_side'] = sample_info[cols].astype(str).apply('_'.join, axis=1)
for sample in sample_info['sample_side'].unique():
  this_sample_info = sample_info.loc[sample_info['sample_side'] == sample]
  nuc_img_path = [np for np in nuclei_mask_paths if np.parts[-1].startswith(sample)]
  golgi_img_path = [gp for gp in golgi_mask_paths if gp.parts[-1].startswith(sample)]
  result_angle_path = [rp for rp in results_paths if rp.parts[-1].startswith(sample)]

  if len(result_angle_path) != 0 and not rerun_pair_finding:
    print('Using previous results for ', sample)
  else:
    if len(nuc_img_path) == 0:
      nuc_sample_name = Path(this_sample_info.loc[this_sample_info['channel'] == 'nuclei']['old_filename'].values[0]).stem
      nuc_img_path = [mp for mp in nuclei_mask_paths if mp.parts[-1].startswith(nuc_sample_name)]
    if len(golgi_img_path) == 0:
      golgi_sample_name = Path(this_sample_info.loc[this_sample_info['channel'] == 'golgi']['old_filename'].values[0]).stem
      golgi_img_path = [mp for mp in golgi_mask_paths if mp.parts[-1].startswith(golgi_sample_name)]
    if len(nuc_img_path) != 0 and len(golgi_img_path) != 0: # don't run the code if there's no matching image...
      # this can't happen with this data but maybe someone will find all the 2d stuff helpful so I'll leave it in
      # you'll have to add a column in the info.csv file called 'dimensions' where you put if it's 2 or 3
      # if this_sample_info['dimensions'].iloc[0] == 2:
      #   nuc_img = cv2.imread(nuc_img_path[0], -1)
      #   golgi_img = cv2.imread(golgi_img_path[0], -1)
      #   dimensions = 2
      #   # calculate moments on each mask and get the centroids out
      #   # NOTE: the masks start at 1 but all other arrays in python start with 0, so need to be mindful of adding/subtracting 1
      #   nuc_centroids = get_all_centroids(nuc_img, dimensions)
      #   golgi_centroids = get_all_centroids(golgi_img, dimensions)
      #   scaling = np.array([1,1])
      # else:

      # this is the opencv way but it isn't great for what we want here
      # r, nuc_img = cv2.imreadmulti(nuc_img_path[0], [], -1)
      print('Analyzing sample', sample)
      nuc_img = tifffile.imread(nuc_img_path[0])
      golgi_img = tifffile.imread(golgi_img_path[0])
      anis = this_sample_info['z_space'].iloc[0]
      scaling = np.array([1, 1, anis])
      dimensions = np.size(nuc_img.shape)

      # calculate moments on each mask and get the centroids out
      # this is extraordinarily slow, so we should save and load them if possible
      print('calculating nuclei centroids...')
      nuc_centroids = get_all_centroids(nuc_img, 3, nuc_img_path, nuclei_centroid_paths, rerun_centroid_finding, anis)
      nuc_radii = nuc_centroids[:,3]
      nuc_centroids = nuc_centroids[:,0:3]
      print('calculating Golgi centroids...')
      golgi_centroids = get_all_centroids(golgi_img, dimensions, golgi_img_path, golgi_centroid_paths, rerun_centroid_finding, anis )
      golgi_radii = golgi_centroids[:,3]
      golgi_centroids = golgi_centroids[:,0:3]

      count, bins = np.histogram(nuc_radii, 100)
      plt.stairs(count, bins)
      plt.title(sample + ' nuclei radii')
      plt.show()
      plt.close()

      count, bins = np.histogram(golgi_radii, 100)
      plt.stairs(count, bins)
      plt.title(sample + ' Golgi radii')
      plt.show()
      plt.close()

      #End Else that was commented out above...

      ### Finding matching pairs
      ## This is another large section of code that is helpful for the 2D analysis but not really needed for 3D - leaving it
      ## in so that maybe someone could find it helpful for a different project or something

      ## First we'll just match stain that are co-localized to nuclei
      # Not doing this first messes with the linear sum assignment later, even though it's way faster than these loops
      # makes a true false mask where the two images are not zero

      # print('finding overlaps')
      # overlapping_stain = np.logical_and(nuc_img, stain_img)
      # nuc_w_stain = np.unique(nuc_img[overlapping_stain])
      # overlap_stain_image = stain_img[overlapping_stain]
      # matches = []
      # for nuc in nuc_w_stain:
      #   matched_stain = np.unique(stain_img[np.where(nuc_img == nuc)])
      #   if len(matched_stain) > 1: # this prevents an error  where the entire stain overlaps with nuclei
      #     matched_stain = matched_stain[1:]
      #   if len(matched_stain) > 1: # this chooses the better match
      #     dist = distance.cdist([nuc_centroids[nuc-1]]*scaling, stain_centroids[matched_stain - 1]*scaling) # calculate distances
      #     matched_stain = matched_stain[dist.argsort()][0] # pick the match with the shortest distance
      #   matched_stain = matched_stain[0] # its a list of lists, get rid of the layering
      #
      #   # nuc_mask = np.where(nuc_img[int(nuc_centroids[nuc][2]), :, :] == (nuc), 1, 0)
      #   # stain_mask = np.where(stain_img[int(stain_centroids[matched_stain][2]), :, :] == (matched_stain), 1, 0)
      #   # test_helper(nuc_mask, stain_mask, nuc_centroids, stain_centroids, [nuc, matched_stain])
      #
      #   matches.append([nuc, matched_stain])
      # matches = np.array(matches) # matches are in image index, not centroid index
      # # Find all the duplicates and delete them
      # dups, dup_id, dup_count = np.unique(matches[:,1], return_inverse=True, return_counts=True)
      # count_mask = dup_count > 1
      # duplicates = dups[count_mask]
      # to_del = np.argwhere(np.isin(matches[:,1], duplicates))
      # first_matches = np.delete(matches, to_del, axis=0)
      ### END 2D removed code

      ## Linear sum assignment problem
      # makes an array of euclidian distances between nuclei (rows) and stain (column)
      print('building cost matrix')
      first_matches = None
      smallest_nuc = this_sample_info['smallest_nuc'].iloc[0] # gets lower limit of nuc radius from csv

      t1 = time.time()
      cost_function_matrix, nuc_centroid_cull, golgi_centroid_cull = build_distance_cost_matrix(nuc_centroids, golgi_centroids, nuc_radii, smallest_nuc, first_matches)
      t2 = time.time()
      print('cost matrix building time: ', t2 - t1) # basically instant

      # solve the linear sum assignment problem
      print('solving cost matrix')
      nuc_cost_ind, golgi_cost_ind = linear_sum_assignment(cost_function_matrix)
      # now we have to back out to get the right index before removing already matched ones
      nuc_ind = nuc_centroid_cull[nuc_cost_ind]
      golgi_ind = golgi_centroid_cull[golgi_cost_ind]
      second_matches = np.vstack((nuc_ind + 1, golgi_ind+1)).T

      ## Histogram to remove long distance (likely incorrect) matches
      try:
        matches = np.concatenate((first_matches, second_matches), axis=0)
      except:
        matches = second_matches
      pair_distances = distance.cdist(nuc_centroids[matches[:,0] - 1], golgi_centroids[matches[:,1] -1]).diagonal()
      matches = matches[~np.isnan(pair_distances)]
      pair_distances = pair_distances[~np.isnan(pair_distances)]
      if math.isnan(this_sample_info["max_distance"].iloc[0]):
        # this shouldn't ever run with this data, and the code here isn't well tested on my end, but again it may be of interest to someone...
        otsu_thresh = otsu_test(pair_distances)
        count, bins = np.histogram(pair_distances, 200, range=(0,200))
        plt.stairs(count, bins)
        plt.show()
        plt.close()
      else:
        otsu_thresh = this_sample_info["max_distance"].iloc[0]

      fig_save_name = sample + '_matches.tif'
      fig_save_path = Path(output_directory) / fig_save_name
      matches_clean = np.delete(matches, np.argwhere(pair_distances >= otsu_thresh), axis=0)
      plot_all_thresh(nuc_img, golgi_img, matches_clean, nuc_centroids, golgi_centroids, pair_distances, otsu_thresh, dimensions, anis, fig_save_path)

      angles, delta_z = get_angles(matches_clean, nuc_centroids, golgi_centroids, dimensions)
      unit_x, unit_y, unit_z = get_unit_vector(matches_clean, nuc_centroids, golgi_centroids, dimensions)

      # create a table of all the results for this image
      results_df = pd.DataFrame({'nuclei_index': matches_clean[:,0],
                                 'GM130_index': matches_clean[:,1],
                                 'nuclei_centroidx': nuc_centroids[matches_clean[:,0] - 1][:,0],
                                 'nuclei_centroidy': nuc_centroids[matches_clean[:,0] - 1][:,1],
                                 'GM130_centroidx': golgi_centroids[matches_clean[:,1] - 1][:,0],
                                 'GM130_centroidy': golgi_centroids[matches_clean[:,1] - 1][:,1],
                                 'distance': pair_distances[np.argwhere(pair_distances < otsu_thresh)][:,0],
                                 'angle': angles,
                                 'delta_z': delta_z,
                                 'unit_x': unit_x,
                                 'unit_y': unit_y,
                                 'unit_z': unit_z})

      # save results as csv
      results_save_name = sample + '_angle_results.csv'
      results_df.to_csv((results_directory / results_save_name))