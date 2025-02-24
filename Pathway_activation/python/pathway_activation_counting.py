import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import skimage.measure
import skimage.io

from read_roi import  read_roi_file # pip install read-roi

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
data_dir = (Path.cwd().parent.parent / 'data' / 'Branches_activation').resolve()
dab_directory = (data_dir / 'analysis' / 'deconvolved').resolve()
dab_mask_directory = (data_dir / 'analysis' / 'cellpose' / 'DAB').resolve()


## Get sample info
# I use a csv file with info on all the images to match them to the samples and if they are nuclei or tunel images
sample_info = pd.read_csv((data_dir / 'samples.csv').resolve())

# thresh is just the mean of the pos and neg selected cells
sample_info['mean_thresh'] = sample_info[['pos1','pos2', 'pos3', 'neg1', 'neg2', 'neg3']].mean(axis=1)

## find all the orignal images, cellpose mask files, and rois
sample_dirs, sample_img_paths = find_all_filepaths(dab_directory, '.tiff')
sample_img_paths = [sp for sp in sample_img_paths if Path(sp).name.endswith('_deconvolved.tif')] # we don't want the thresh ones
dab_dirs, dab_mask_paths = find_all_filepaths(dab_mask_directory, '.png')
dab_mask_paths = [mp for mp in dab_mask_paths if Path(mp).name.endswith('_masks.png')]
roi_dirs, roi_paths = find_all_filepaths(dab_directory, '.roi')

## Unfortunately this analysis will be a bit more complicated than the proliferation analysis.
## There is a lot of apoptosis in the neurectoderm and ectoderm that we don't want to quantify as being an effect of the
## inhibitors. I will draw rois in ImageJ, import them here, and count how many Cellpose labels it overlaps with.

## we're going to build the results in a loop, then convert to a dataframe for saving
results = []

for sample in sample_info['new_filename'].unique():
  # match the sample name to the correct images and roi file
  dab_img_path = [mp for mp in dab_mask_paths if Path(mp).name.startswith(sample)]
  sample_img_path = [sp for sp in sample_img_paths if Path(sp).name.startswith(sample)]
  roi_path = [rp for rp in roi_paths if Path(rp).name.startswith(sample)]

  if len(roi_path) != 0:
    # load the images and roi file
    # cv2 isn't playing nice with the sample_img files. It keeps trying to open them as color images and if I force it to
    # treat them as greyscale the intensity values are incorrect. Will use skimage instead.
    dab_img = skimage.io.imread(str(dab_img_path[0]))
    sample_img = skimage.io.imread(str(sample_img_path[0]))
    roi = import_roi(str(roi_path[0]))
    roi = roi[roi_path[0].stem]

    # apply the imageJ hand drawn rois here
    # create a mask from the roi drawing
    mask = roi_to_mask(roi)

    dab_img = dab_img * mask

    # calculate the mean intensity for each cell
    label_props = skimage.measure.regionprops(dab_img, sample_img)
    means = np.asarray([d['intensity_mean'] for d in label_props])

    # get all the info into the results table
    row_indices = sample_info.loc[sample_info['new_filename'] == sample]
    result = []
    result.append(row_indices.iloc[0]['new_filename'])
    result.append(row_indices.iloc[0]['treatment'])
    result.append(row_indices.iloc[0]['time'])
    result.append(row_indices.iloc[0]['AB'])
    result.append(row_indices.iloc[0]['new_id'])
    result.append(row_indices.iloc[0]['side'])

    # total cells
    result.append(means.size)
    # positive cells. sum() only counts True value
    thresh = row_indices['mean_thresh'][0]
    result.append( (means < thresh).sum() )

    results.append(result)

results_df = pd.DataFrame(results, columns = ['new_filename', 'treatment', 'time', 'AB', 'new_id', 'side', 'total_cells', 'DAB_pos'])

## save
results_df.to_csv((data_dir / 'results.csv').resolve())
