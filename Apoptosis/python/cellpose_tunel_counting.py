import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

from read_roi import  read_roi_file # pip install read-roi

from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


# ImageJ ROI helpers
def import_roi(roi_file):
  # imports a '.roi' file, not a .zip file
  # needs to be a string and not a path
  roi = read_roi_file(roi_file) # this is from read-roi package
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
  freehand_coords = np.stack((xvals, yvals), axis=1).flatten()
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
data_dir = (Path.cwd().parent.parent / 'data' / 'Apoptosis').resolve()
nuclei_directory = (data_dir / 'cellpose' / 'nuclei').resolve()
tunel_directory = (data_dir / 'cellpose' / 'TUNEL').resolve()
roi_directory = (data_dir / 'cellpose' / 'rois').resolve()

## Get sample info
# I use a csv file with info on all the images to match them to the samples and if they are nuclei or tunel images
sample_info = pd.read_csv((data_dir / 'TUNEL_image_log.csv').resolve())
cols = ['sample', 'side', 'section', 'mag']
sample_info['sample_side'] = sample_info[cols].astype(str).apply('_'.join, axis=1)

## find all the cellpose mask files created in 'run_cellpose_TUNEL.py'
nuclei_dirs, nuclei_mask_paths = find_all_filepaths(nuclei_directory, '.png')
nuclei_mask_paths = [mp for mp in nuclei_mask_paths if Path(mp).name.endswith('_masks.png')]
tunel_dirs, tunel_mask_paths = find_all_filepaths(tunel_directory, '.png')
tunel_mask_paths = [mp for mp in tunel_mask_paths if Path(mp).name.endswith('_masks.png')]
roi_dirs, roi_paths = find_all_filepaths(roi_directory, '.roi')

## Unfortunately this analysis will be a bit more complicated than the proliferation analysis.
## There is a lot of apoptosis in the neurectoderm and ectoderm that we don't want to quantify as being an effect of the
## inhibitors. I will draw rois in ImageJ, import them here, and count how many Cellpose labels it overlaps with.

## we're going to build the results in a loop, then convert to a dataframe for saving
results = []

for sample in sample_info['sample_side'].unique():
  # match the sample name to the correct images and roi file
  nuc_img_path = [mp for mp in nuclei_mask_paths if Path(mp).name.startswith(sample)]
  stain_img_path = [mp for mp in tunel_mask_paths if Path(mp).name.startswith(sample)]
  roi_path = [rp for rp in roi_paths if Path(rp).name.startswith(sample)]

  if len(roi_path) != 0:
    # load the images and roi file
    nuc_img = cv2.imread(str(nuc_img_path[0]), -1) # https://docs.opencv.org/3.4/d8/d6a/group__imgcodecs__flags.html
    stain_img = cv2.imread(str(stain_img_path[0]), -1)
    roi = import_roi(str(roi_path[0])) # import the roi as a dictionary
    roi = roi[roi_path[0].stem] # get the roi out of the dictionary - usually there's only one, but this will get it out

    # create a mask from the roi drawing
    mask = roi_to_mask(roi)

    # img * mask multiples all the image pixel values by either 0 or 1 to delete areas outside of roi
    # np.unique counts all the unique values of pixels, in this case cell labels
    # size calculates the total number of unique cell labels
    # we subtract 1 because it always detects the background as a unique value, but this isn't a cell... it's background
    masked_nuc_count = np.unique(nuc_img * mask).size - 1
    masked_tunel_count = np.unique(stain_img * mask).size - 1

    row_indices = sample_info.loc[sample_info['sample_side'] == sample]
    result = []
    result.append('; '.join([str for str in row_indices.old_filename])) # another 'wonderful' one liner
    result.append(Path(nuc_img_path[0]).name)
    result.append(row_indices.iloc[0]['sample'])
    result.append(row_indices.iloc[0]['treatment'])
    result.append(row_indices.iloc[0]['side'])
    result.append(row_indices.iloc[0]['section'])
    result.append(row_indices.iloc[0]['mag'])
    result.append(masked_nuc_count)
    result.append(masked_tunel_count)

    results.append(result)

results_df = pd.DataFrame(results, columns = ['old_filenames', 'new_filename', 'sample', 'treatment', 'side', 'section',
                                              'mag', 'nuclei_count', 'tunel_count'])

## save
results_df.to_csv((data_dir / 'results.csv').resolve())
