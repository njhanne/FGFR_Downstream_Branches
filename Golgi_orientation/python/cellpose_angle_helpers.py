# from cellpose import plot

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from windrose import WindroseAxes, plot_windrose

import cv2
from skimage import measure
import tifffile

from scipy.spatial import distance
from scipy.optimize import linear_sum_assignment
import math

import os
from pathlib import Path
# import time
# import re
import tkinter as tk
# from tkinter import filedialog

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths

# File Loading and DataFrame Results Helpers
def initialize_results_dataframe(input_csv, length):
  df = pd.DataFrame(columns=input_csv.columns)
  df.drop(columns =['invert'])
  df['nuclei_count']=''
  df['stain_count']=''
  return df


def get_all_vals(d):
  # gets all the values in nested dictionary (d)
  vals = []
  for k, v in d.items():
    if isinstance(v, dict):
      vals.extend(get_all_vals(v))
    else:
      vals.append(v)
  return vals


# End DataFrame Helpers

# Image Processing and Geometry Helpers
def get_all_centroids(image, dimension, img_path, centroid_paths):
  save_name = os.path.splitext(img_path[0])[0] + '_centroids.npy'
  centroid_path = [cp for cp in centroid_paths if cp.endswith(save_name)]
  if len(centroid_path) != 0:
    if os.path.getmtime(centroid_path[0]) > os.path.getmtime(img_path[0]):
      # if the centroid file is newer than the last time the cellpose image was generated then load it
      centroids = np.load(centroid_path[0])
      print("Loading previous centroids...")
    else:
      # if the centroid file is older than the new cellpose image then regenerate it
      centroids = get_centroids(image, dimension)
      np.save(save_name, centroids)
  else:
    centroids = get_centroids(image, dimension)
    np.save(save_name, centroids)
  np.save(save_name, centroids)
  return centroids


def get_centroids(image, dimension):
  centroids = []
  max = np.amax(image)+1
  for i in range(max):
    if i != 0:
      if dimension == 2:
        #TODO add bounding box for 2D
        centroids.append(calculate_centroid(np.where(image == i, 1, 0)))
      if dimension == 3:
        print("Finding centroid " + str(i) +  " of " + str(max))
        xbox, ybox, zbox, bbox_img = bbox2_3D(np.where(image == i, 1, 0))
        centroidx, centroidy, centroidz, radius = calculate_3D_centroid(bbox_img)
        centroids.append((centroidx + xbox, centroidy + ybox, centroidz + zbox, radius))
  return np.array(centroids)


def bbox2_3D(img):
  # https: // stackoverflow.com / a / 31402351 / 19260308
  # this is at least 3x faster than not using a bounding box
  z = np.any(img, axis=(1, 2))
  y = np.any(img, axis=(0, 2))
  x = np.any(img, axis=(0, 1))

  xmin, xmax = np.where(x)[0][[0, -1]]
  ymin, ymax = np.where(y)[0][[0, -1]]
  zmin, zmax = np.where(z)[0][[0, -1]]

  cropped_img = img[zmin:zmax, ymin:ymax, xmin:xmax]
  return xmin, ymin, zmin, cropped_img


def old_get_centroids(image, dimension):
  centroids = []
  max = np.amax(image)+1
  for i in range(max):
    if i != 0:
      if dimension == 2:
        centroids.append(calculate_centroid(np.where(image == i, 1, 0)))
      if dimension == 3:
        print("Finding centroid " + str(i) +  " of " + str(max))
        centroids.append(calculate_3D_centroid(np.where(image == i, 1, 0)))
  return np.array(centroids)


def calculate_centroid(image):
  # pass in the mask from cellpose .png and it will output the xy coordinates of the centroid
  M = cv2.moments(image.astype('uint8')) # moments needs uint8 for some reason...
  centroid_X = M["m10"] / M["m00"]
  centroid_Y = M["m01"] / M["m00"]
  return centroid_X, centroid_Y


def calculate_3D_centroid(image):
  # openCV doesn't handle more than 2D so we have to use skimage
  # skimage doesn't handle z-spacing right now, but I dont think we need it for centroid
  # pass in the mask from cellpose or trackmate .tif and it will output the xyz coordinates of the centroid
  # I think skimage does z,y,x
  M = measure.moments(image.astype('uint8'), order=3) # moments needs uint8 for some reason...
  centroid_X = M[0,0,1] / M[0,0,0]
  centroid_Y = M[0,1,0] / M[0,0,0]
  centroid_Z = M[1,0,0] / M[0,0,0]
  # 0th moment is the volume, volume of sphere is 4/3 * pi * r^3
  radius = ((3*M[0,0,0])/(4*math.pi))**(1/3)
  return centroid_X, centroid_Y, centroid_Z, radius


def build_distance_cost_matrix(nuc_centroids, stain_centroids, scaling, nuc_radii, smallest_nuc, matches = None):
  # make a mask same size as the list of centroids
  nuc_centroid_cull = np.ones(len(nuc_centroids[:, 0]))
  stain_centroid_cull = np.ones(len(stain_centroids[:, 0]))

  # if there are already matches, filter them out of the mask so they aren't included in matrix
  if matches is not None:
    # matches are in image order... if there are x nuclei, the centroid index would be x-1
    nuc_centroid_cull[matches[:,0] - 1] = 0
    stain_centroid_cull[matches[:, 1] - 1] = 0
  # filter out very small nuclei (probably bad cellpose segmentation)
  if smallest_nuc is not None:
    nuc_centroid_cull_first = np.ones(len(nuc_centroids[:, 0]))
    nuc_centroid_cull_first[nuc_radii < smallest_nuc] = 0
    nuc_centroid_cull_first = np.where(nuc_centroid_cull_first)[0]


  nuc_centroid_cull_second = np.arange(0,len(nuc_centroids),1)
  nuc_centroid_cull_second = (np.delete(nuc_centroid_cull_second, np.searchsorted(nuc_centroid_cull_second, np.unique(np.where(np.isnan(nuc_centroids))[0]))))
  nuc_centroid_cull = np.intersect1d(nuc_centroid_cull_first, nuc_centroid_cull_second)
  stain_centroid_cull = np.where(stain_centroid_cull)
  stain_centroid_cull = (np.delete(stain_centroid_cull[0], np.searchsorted(stain_centroid_cull[0], np.unique(np.where(np.isnan(stain_centroids))[0]))))

  # calculate distance between centroids
  cost_matrix = distance.cdist(nuc_centroids[nuc_centroid_cull]*scaling, stain_centroids[stain_centroid_cull]*scaling)
  return cost_matrix, nuc_centroid_cull, stain_centroid_cull


def otsu_test(image):
  # This usually only works with images, but this code works on any data
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

def get_angles(matches, nuc_centroids, stain_centroids, dimension):
  # angles are ccw starting at 0 on the positive x-axis
  angles = []
  delta_z = []
  for pair in matches:
    if dimension == 2:
      x0, y0 = nuc_centroids[pair[0]-1]
      x1, y1 = stain_centroids[pair[1]-1]
      dz = 0
    else:
      x0, y0, z0 = nuc_centroids[pair[0]-1]
      x1, y1, z1 = stain_centroids[pair[1]-1]
      dz = z1 - z0
    angle = -math.atan2(y1 -y0, x1-x0)
    if angle < 0:
      angle = angle + 2*math.pi
    if angle > 2*math.pi:
      angle = angle - 2*math.pi
    angles.append(angle)
    delta_z.append(dz)
  return np.array(angles), delta_z


def get_unit_vector(matches, nuc_centroids, stain_centroids, dimension):
  # get a unit vector with 0,0,0 as nucleus centroid and terminating at golgi centroid
  x_hat = []
  y_hat = []
  z_hat = []
  for pair in matches:
    if dimension == 2:
      x0, y0 = nuc_centroids[pair[0]-1]
      x1, y1 = stain_centroids[pair[1]-1]
      dx = x1 - x0
      dy = y1 - y0
      dz = 0
    else:
      x0, y0, z0 = nuc_centroids[pair[0]-1]
      x1, y1, z1 = stain_centroids[pair[1]-1]
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
def plot_2D_mask(image, mask=None):
  if mask is not None:
    image = np.where(image == mask, 1, 0)
  fig, ax = plt.subplots(figsize=(20, 20))
  ax.imshow(image)
  x0, y0 = calculate_centroid(image)
  ax.plot(x0, y0, '.g', markersize=4)
  ax.axis((0, 2048, 2048, 0))
  plt.show()


def plot_pairs(nuc_img, stain_img, pair):
  nuc = np.where(nuc_img == pair[0], 1, 0)
  stain = np.where(stain_img == pair[1], 2, 0)
  image = nuc + stain
  fig, ax = plt.subplots(figsize=(20, 20))
  ax.imshow(image)
  x0, y0 = calculate_centroid(nuc)
  x1,y1 = calculate_centroid(stain)
  ax.plot(x0, y0, '.g', markersize=4)
  ax.plot(x1,y1, '.r', markersize=4)
  ax.plot((x0, x1), (y0, y1), '-r', linewidth=3)
  ax.axis((0, 2048, 2048, 0))
  plt.show()


def plot_all(nuc_img, stain_img, pairs, nuc_centroids, stain_centroids):
  nuc = np.where(nuc_img != 0, 1, 0)
  stain = np.where(stain_img != 0, 2, 0)
  image = nuc + stain
  fig, ax = plt.subplots(figsize=(20, 20))
  ax.imshow(image)
  for pair in pairs:
    x0, y0 = nuc_centroids[pair[0]-1]
    x1, y1 = stain_centroids[pair[1]-1]
    ax.plot((x0, x1), (y0, y1), '-w', linewidth=1)
  ax.axis((0, 2048, 2048, 0))
  plt.show()

def test_helper(nuc_mask, stain_mask, nuc_centroids, stain_centroids, pair):
  image = nuc_mask + stain_mask * 2

  fig2, ax2 = plt.subplots(figsize=(20, 20))
  ax2.axis((0, 2048, 2048, 0))
  ax2.imshow(image)

  x0, y0, z0 = nuc_centroids[pair[0] - 1]
  x1, y1, z0 = stain_centroids[pair[1] - 1]
  ax2.plot((x0, x1), (y0, y1), '-w', linewidth=1)

  plt.show()
  plt.close()



def plot_all_thresh(nuc_img, stain_img, pairs, nuc_centroids, stain_centroids, distances, otsu, dimension, save_path=None):
  if dimension == 2:
    nuc = np.where(nuc_img != 0, 1, 0)
    stain = np.where(stain_img != 0, 2, 0)
    image = nuc + stain
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(image)
    for ind, pair in enumerate(pairs):
      x0, y0 = nuc_centroids[pair[0] - 1]
      x1, y1 = stain_centroids[pair[1] - 1]
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
    stain_masks = np.zeros((xpix, ypix))
    for ind, pair in enumerate(pairs):
      # if (pair[0] != 1 & pair[1] != 1):
      if (pair[0] != nuc_centroids.shape[0]) & (pair[1] != stain_centroids.shape[0]):
        if not np.isnan(nuc_centroids[pair[0]-1][0]):
          nuc_mask = np.where(nuc_img[int(nuc_centroids[pair[0]-1][2]),:,:] == (pair[0]), 1, 0)
          stain_mask = np.where(stain_img[int(stain_centroids[pair[1]-1][2]),:,:] == (pair[1]), 1, 0)
          nuc_masks = nuc_masks.astype(int) | nuc_mask.astype(int)
          stain_masks = stain_masks.astype(int) | stain_mask.astype(int)
          # test_helper(nuc_mask, stain_mask, nuc_centroids, stain_centroids, pair)
          # print('done')
    image = nuc_masks + stain_masks*2
    ax.imshow(image)
    for ind, pair in enumerate(pairs):
      x0, y0, z0 = nuc_centroids[pair[0] - 1]
      x1, y1, z0 = stain_centroids[pair[1] - 1]
      # if distances[ind] > otsu:
      #   # ax.plot((x0, x1), (y0, y1), '-r', linewidth=2)
      # else:
      ax.plot((x0, x1), (y0, y1), '-w', linewidth=1)
  if save_path is not None:
    plt.savefig(save_path)
  plt.show()
  plt.close()


def plot_distance_histogram(count, count2=None):
  n, bins, patches = plt.hist(count, 16, density=True, facecolor='g', alpha=0.75)
  if count2 is not None:
    plt.hist(count2, 16, density=True, facecolor='r', alpha=0.5)
  plt.show()


def plot_windrose(direction, save_path=None):
  ax = WindroseAxes.from_ax()

  # windrose measures cw from positive y direction, need to adjust...
  direction = -direction + np.pi/2
  # have to convert them so they are all 0-360 # for some reason they use degree instead of rad
  direction = np.degrees((direction + (2 * np.pi)) % (2 * np.pi))

  # 'mirror' the angles so that you get a full windrose
  # direction_flip = [angle + 180 if angle < 180 else angle - 180 for angle in direction]
  # direction = np.append(direction, direction_flip)

  # fill in the 'var' with 1s so that it plots one color
  fake_var = np.ones(len(direction))
  ax.bar(direction=direction, var=fake_var, opening=1, bins=np.arange(-.5, .5, 1), nsector=36)
  if save_path is not None:
    plt.savefig(save_path)
  plt.show()
  plt.close()


### MAIN ###
# First find all needed directories and load them all up
root = tk.Tk()
root.withdraw()

# nuclei_directory = filedialog.askdirectory(title='Select masked nuclei directory')
nuclei_directory = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/cellpose_output/nuclei')
# stain_directory = filedialog.askdirectory(title='Select other stain masked directory')
stain_directory = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/cellpose_output/golgi')
# output_directory = filedialog.askdirectory(title='Select output directory')
output_directory = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/cellpose_output/analysis_output')
# info_csv_path = filedialog.askopenfilename(title='Select sample log')
# info_csv_path = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/GM130_image_log.csv')
# info_csv_path = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/GM130_image_log_temp.csv')
info_csv_path = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/GM130_image_log_temp.csv')


sample_info = pd.read_csv(info_csv_path)

nuclei_dirs, nuclei_mask_paths = find_all_filepaths(nuclei_directory, '.tif')
nuclei_mask_paths = [mp for mp in nuclei_mask_paths if mp.endswith('_masks.tif')]
stain_dirs, stain_mask_paths = find_all_filepaths(stain_directory, '.tif')
stain_mask_paths = [mp for mp in stain_mask_paths if mp.endswith('_masks.tif')]

nuc_centroid_dirs, nuclei_centroid_paths = find_all_filepaths(nuclei_directory, '.npy')
nuclei_centroid_paths = [cp for cp in nuclei_centroid_paths if cp.endswith('_centroids.npy')]
stain_centroid_dirs, stain_centroid_paths = find_all_filepaths(stain_directory, '.npy')
stain_centroid_paths = [cp for cp in stain_centroid_paths if cp.endswith('_centroids.npy')]


# results_df = initialize_results_dataframe(sample_info, len(stain_mask_paths.count))
# results = []

cols = ['sample', 'side', 'section', 'view']
sample_info['sample_side'] = sample_info[cols].astype(str).apply('_'.join, axis=1)
for sample in sample_info['sample_side'].unique():
  this_sample_info = sample_info.loc[sample_info['sample_side'] == sample]
  print('Analyzing sample', sample)
  nuc_img_path = [mp for mp in nuclei_mask_paths if os.path.basename(mp).startswith(sample)]
  stain_img_path = [mp for mp in stain_mask_paths if os.path.basename(mp).startswith(sample)]
  if len(nuc_img_path) == 0:
    nuc_sample_name = Path(this_sample_info.loc[this_sample_info['channel'] == 'nuclei']['old_filename'].values[0]).stem
    nuc_img_path = [mp for mp in nuclei_mask_paths if os.path.basename(mp).startswith(nuc_sample_name)]
  if len(stain_img_path) == 0:
    stain_sample_name = Path(this_sample_info.loc[this_sample_info['channel'] == 'golgi']['old_filename'].values[0]).stem
    stain_img_path = [mp for mp in stain_mask_paths if os.path.basename(mp).startswith(stain_sample_name)]
  if this_sample_info['dimensions'].iloc[0] == 2:
    nuc_img = cv2.imread(nuc_img_path[0], -1)
    stain_img = cv2.imread(stain_img_path[0], -1)

    # calculate moments on each mask and get the centroids out
    # NOTE: the masks start at 1 but all other arrays in python start with 0, so need to be mindful of adding/subtracting 1
    nuc_centroids = get_all_centroids(nuc_img, dimension=2)
    stain_centroids = get_all_centroids(stain_img, dimension=2)
    scaling = np.array([1,1])
  else:
    # this is the opencv way but it isn't great for what we want here
    # temp, nuc_img = cv2.imreadmulti(nuc_img_path[0], [], -1)

    nuc_img = tifffile.imread(nuc_img_path[0])
    stain_img = tifffile.imread(stain_img_path[0])
    scaling = np.array([1, 1, this_sample_info['z_space'].iloc[0]])

    # calculate moments on each mask and get the centroids out
    # this is extraordinarily slow, so we should save and load them if possible
    nuc_centroids = get_all_centroids(nuc_img, 3, nuc_img_path, nuclei_centroid_paths)
    nuc_radii = nuc_centroids[:,3]
    nuc_centroids = nuc_centroids[:,0:3]
    stain_centroids = get_all_centroids(stain_img, 3, stain_img_path, stain_centroid_paths)
    stain_radii = stain_centroids[:,3]
    stain_centroids = stain_centroids[:,0:3]

    count, bins = np.histogram(nuc_radii, 100)
    plt.stairs(count, bins)
    plt.show()
    plt.close()

    count, bins = np.histogram(stain_radii, 100)
    plt.stairs(count, bins)
    plt.show()
    plt.close()

  ### Finding matching pairs
  ## First we'll just match stain that are co-localized to nuclei
  # Not doing this first messes with the linear sum assignment later, even though it's way faster than these loops
  # makes a true false mask where the two images are not zero
  #TODO: I don't know if this is needed for the 3D data - maybe test removing it as it's slow af

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

  ## Linear sum assignment problem
  # makes an array of euclidian distances between nuclei (rows) and stain (column) but only the ones that overlap...
  print('building cost matrix')
  first_matches = None
  smallest_nuc = 4
  cost_function_matrix, nuc_centroid_cull, stain_centroid_cull = build_distance_cost_matrix(nuc_centroids, stain_centroids, scaling, nuc_radii, smallest_nuc, first_matches)

  # solve the linear sum assignment problem
  print('solving cost matrix')
  nuc_cost_ind, stain_cost_ind = linear_sum_assignment(cost_function_matrix)
  # now we have to back out to get the right index before removing already matched ones
  nuc_ind = nuc_centroid_cull[nuc_cost_ind]
  stain_ind = stain_centroid_cull[stain_cost_ind]
  second_matches = np.vstack((nuc_ind + 1, stain_ind+1)).T

  ## Histogram to remove long distance (likely incorrect) matches
  try:
    matches = np.concatenate((first_matches, second_matches), axis=0)
  except:
    matches = second_matches
  pair_distances = distance.cdist(nuc_centroids[matches[:,0] - 1]*scaling, stain_centroids[matches[:,1] -1]*scaling).diagonal()
  matches = matches[~np.isnan(pair_distances)]
  pair_distances = pair_distances[~np.isnan(pair_distances)]
  if math.isnan(this_sample_info["max_distance"].iloc[0]):
    # otsu_thresh = otsu_test(pair_distances)
    # TODO remove this
    count, bins = np.histogram(pair_distances, 200, range=(0,200))
    plt.stairs(count, bins)
    plt.show()
    plt.close()
    otsu_thresh = 70
  else:
    otsu_thresh = this_sample_info["max_distance"].iloc[0]
    # otsu_thresh = 70


  fig_save_name = sample + '_matches.tif'
  fig_save_path = Path(output_directory) / fig_save_name
  matches_clean = np.delete(matches, np.argwhere(pair_distances >= otsu_thresh), axis=0)
  plot_all_thresh(nuc_img, stain_img, matches_clean, nuc_centroids, stain_centroids, pair_distances, otsu_thresh, this_sample_info["dimensions"].iloc[0], fig_save_path)

  angles, delta_z = get_angles(matches_clean, nuc_centroids, stain_centroids, this_sample_info["dimensions"].iloc[0])
  unit_x, unit_y, unit_z = get_unit_vector(matches_clean, nuc_centroids, stain_centroids, this_sample_info["dimensions"].iloc[0])

  fig_save_name = sample + '_windrose.png'
  fig_save_path = Path(output_directory) / fig_save_name
  plot_windrose(angles, fig_save_path)

  ## Old way - will keep for ref
  # # gets the mask index in the stain image where there is overlap with nuclei, then subtracts 1 to be same index as the centroids data
  # overlapping_stain_ind = np.unique(stain_img[overlapping_stain]) - 1
  # # makes an array of euclidian distances between nuclei (rows) and stain (column) but only the ones that overlap...
  # cost_function_matrix_first_run = distance.cdist(nuc_centroids, stain_centroids[overlapping_stain_ind])
  # # solve the linear sum assignment problem
  # nuc_ind, stain_cost_ind = linear_sum_assignment(cost_function_matrix_first_run)
  # # now we have to back out to get the right index for the overlaping stains
  # # first convert from distance matrix index (cost function) to
  # stain_ind = overlapping_stain_ind[stain_cost_ind]
  # matches = np.vstack((nuc_ind + 1, stain_ind+1)).T

  # create a table of all the results for this image
  results_df = pd.DataFrame({'nuclei_index': matches_clean[:,0],
                             'GM130_index': matches_clean[:,1],
                             'nuclei_centroidx': nuc_centroids[matches_clean[:,0] - 1][:,0],
                             'nuclei_centroidy': nuc_centroids[matches_clean[:,0] - 1][:,1],
                             'GM130_centroidx': stain_centroids[matches_clean[:,1] - 1][:,0],
                             'GM130_centroidy': stain_centroids[matches_clean[:,1] - 1][:,1],
                             'distance': pair_distances[np.argwhere(pair_distances < otsu_thresh)][:,0],
                             'angle': angles,
                             'delta_z': delta_z,
                             'unit_x': unit_x,
                             'unit_y': unit_y,
                             'unit_z': unit_z})

  # save results as csv
  results_save_name = sample + '_angle_results.csv'
  results_df.to_csv(os.path.normpath(os.path.dirname(info_csv_path) + '\\' + results_save_name))

