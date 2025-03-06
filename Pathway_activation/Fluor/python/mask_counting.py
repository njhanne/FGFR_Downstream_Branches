import numpy as np
import pandas as pd

import skimage.measure
import skimage.io

from PIL import Image, ImageDraw
import shapely
from read_roi import read_roi_file # pip install read-roi

# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt

from pathlib import Path
import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


# ImageJ ROI helpers
def import_roi(roi_file):
  # imports a '.roi' file, not a .zip file
  # needs to be a string and not a path
  roi = read_roi_file(roi_file)
  return roi


def initialize_mask(image=None):
  # Get image size
  # all our images are 2048 square
  # don't need to run this every time

  if image is not None:
    xsize = image.sizes['y']
    ysize = image.sizes['x']
  else:
    xsize = ysize = 2048
  mask_img = Image.new('1', (xsize, ysize), 0)  # makes a blank image w/ same dimensions as image
  return mask_img


def roi_to_mask(roi, image=None):
  mask_img = initialize_mask(image)
  if roi['type'] == 'composite':
    mask = composite_roi_mask(mask_img, roi)
  elif roi['type'] == 'freehand':
    mask = freehand_roi_mask(mask_img, roi)
  return mask


def freehand_roi_mask(mask_img, roi):
  # Get pixel values of freehand shape
  xvals = np.rint(roi['y']).astype(int)
  yvals = np.rint(roi['x']).astype(int)
  freehand_coords = np.stack((xvals,yvals), axis=1).flatten()
  freehand_coords = np.roll(freehand_coords, 1)  # not sure why it leaves these hangers-on?
  test = ImageDraw.Draw(mask_img) # this draws the mask onto the mask_image - it automatically writes
  test.polygon(freehand_coords.tolist(), outline=1, fill=1)  # draws freehand shape and fills with 1's
  return np.array(mask_img)


def composite_roi_mask(mask_img, roi):
  # unfortunately this is kind of roughly done
  # if I draw multiple rois in ImageJ and then combine them, they get saved as a nested list of freehand rois
  # there's no way to know which list of values is the desired roi and which is the hole cutout (in this case the bead)
  # to get around that I will assume all the bead rois are smaller than the overall roi, so I will rearrange them
  # so that the desired roi gets drawn first (fill with 1s), then the hole will be drawn on top of it (fill with 0s)
  # It should work for this, but I wouldn't count on this code being usable outside of this specific use case
  roi['paths'].sort(key=len, reverse=True) # sorts the lists from long to short, so that big roi is (hopefully) first
  for fh_path_num in range(len(roi['paths'])):
    test = ImageDraw.Draw(mask_img) # this draws the mask onto the mask_image - it automatically writes
    if fh_path_num == 0:
      test.polygon(roi['paths'][0], outline=1, fill=1)
    else:
      test.polygon(roi['paths'][fh_path_num], outline=1, fill=0)
  return np.array(mask_img)


# End ImageJ ROI helpers

### MAIN ###
# First find all needed directories and load them all up
data_dir = (Path.cwd().parent.parent.parent / 'data' / 'Branches_activation').resolve()

sample_img_directory = (data_dir / 'analysis' / 'raw' / 'Fluor' / 'activation').resolve()
sample_mask_directory = (data_dir / 'analysis' / 'cellpose' / 'Fluor' / 'expanded_mask').resolve()
sample_roi_directory = (data_dir / 'analysis' / 'raw' / 'Fluor' / 'Nuc').resolve()
output_dir = (data_dir / 'analysis').resolve()


## Get sample info
sample_info = pd.read_csv((data_dir / 'fluor_samples.csv').resolve())


## find all the orignal images, cellpose mask files, and rois
sample_dirs, img_paths = find_all_filepaths(sample_img_directory, '.tif')
mask_dirs, mask_paths = find_all_filepaths(sample_mask_directory, '.tif')
roi_dirs, roi_paths = find_all_filepaths(sample_roi_directory, '.roi')

## we're going to build the results in a loop, then convert to a dataframe for saving
results = []


for sample in sample_info['sample_name'].unique():
  # match the sample name to the correct images and roi file
  row_indices = sample_info.loc[sample_info['sample_name'] == sample]
  nuc_filename = row_indices.loc[row_indices['channel'] == 'DAPI']['filename'].values[0]
  activation_row = row_indices.loc[row_indices['channel'] == 'activation']

  mask_path = [mp for mp in mask_paths if Path(mp).name.startswith(nuc_filename[:-4])]
  img_path = [sp for sp in img_paths if Path(sp).name.startswith(activation_row['filename'].values[0][:-4])]
  roi_path =[rp for rp in roi_paths if Path(rp).name.startswith(nuc_filename[:-4])]

  if len(mask_path) != 0:
    # load the images
    img = skimage.io.imread(str(img_path[0]))
    cell_mask = skimage.io.imread(str(mask_path[0]))

    print(sample)

    roi = read_roi_file(str(roi_path[0]))
    roi = list(roi.values())[0]
    # https://shapely.readthedocs.io/en/stable/manual.html#linestrings
    line = LineString([(i[0], i[1]) for i in zip(roi['x'], roi['y'])])


    # calculate the mean intensity for each cell
    label_props = skimage.measure.regionprops(cell_mask, img)
    means = np.asarray([d['intensity_mean'] for d in label_props])
    means = np.insert(means, 0, 0, axis=0)

    # calculate thresh and filter dim cells
    filter_thresh = activation_row.iloc[:,-6:].mean(axis=1).values[0] # average the 6 good/bad cells
    n_pos_cells = (means > filter_thresh).sum() # positive cells. sum() only counts True value
    n_total_cells = means.size

    # optional code for viewing the results
    # cell_LUT = np.insert(means, 0, 0, axis=0)
    # all_cells_img = cell_LUT[cell_mask] # view the masked cells but with the mean green channel as their pixel value
    # filtered_LUT = cell_LUT
    # filtered_LUT[filtered_LUT < filter_thresh] = 0
    # filtered_cells_img = filtered_LUT[cell_mask] # view the masked cells but with the mean green channel as their pixel value

    # get all the info into the results table
    result = []
    result.append(activation_row.iloc[0]['filename'])
    result.append(activation_row.iloc[0]['treatment'])
    result.append(activation_row.iloc[0]['pathway'])
    result.append(activation_row.iloc[0]['sample_name'])


    # cell counts
    result.append(filter_thresh)
    result.append(n_total_cells)
    result.append(n_pos_cells)

    results.append(result)

results_df = pd.DataFrame(results, columns = ['old_filename', 'treatment', 'pathway', 'sample_name', 'thresh', 'total_cells', 'pos_cells'])

## save
results_df.to_csv((data_dir / 'fluor_results.csv').resolve())