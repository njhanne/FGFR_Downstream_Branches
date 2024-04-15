# I noticed some of the images in the phh3 stack appear to be from the same sample, just different views
# This code will help detect those similarities to find duplicates
# https://stackoverflow.com/a/67477060
# https://docs.opencv.org/4.x/d4/dc6/tutorial_py_template_matching.html

import pandas as pd
import cv2

from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


# Image Manipulation Helpers
def confidence(img, template):
  res = cv2.matchTemplate(img, template, cv2.TM_CCOEFF_NORMED)
  conf = res.max()
  return conf


def compare_corners(template, match_image, size = 2):
  ysize = int(match_image.shape[0] / size)
  xsize = int(match_image.shape[1] / size)
  conf = 0
  for i in range(0,4): # there is almost certainly a better way to do this, but for now this is fine
    if i == 0:
      match_corner = match_image[0:ysize, 0:xsize]
    if i == 1:
      match_corner = match_image[0:ysize, match_image.shape[1]-xsize:match_image.shape[1]]
    if i == 2:
      match_corner = match_image[match_image.shape[0]-ysize:match_image.shape[0], 0:xsize]
    if i == 3:
      match_corner = match_image[match_image.shape[0]-ysize:match_image.shape[0], match_image.shape[1]-xsize:match_image.shape[1]]
    conf_temp = confidence(match_corner, template)
    if (conf_temp > conf): conf = conf_temp
  return conf


# End Image Manipulation Helpers

### MAIN ###
# start by getting all the raw images
data_dir = (Path.cwd().parent.parent / 'data' / 'proliferation').resolve()
input_directory = (data_dir / 'new_images2').resolve()
image_dirs, image_paths = find_all_filepaths(input_directory, '.tif')


## Get sample info
# I use a csv file with info on all the images to match them to the samples and if they are nuclei or phh3 images
image_info = pd.read_csv((data_dir / 'pHH3_image_log.csv').resolve())
cols = ['sample', 'side']
image_info['sample_side'] = image_info[cols].astype(str).apply('_'.join, axis=1) # combine sample and side info into one column
cols = ['group', 'side']
image_info['group_side'] = image_info[cols].astype(str).apply('_'.join, axis=1)

names_l = []
conf_l = []
for groups in image_info['group_side'].unique():
  group_info = image_info.loc[(image_info['group_side'] == groups) & (image_info['invert'] == 'no')]

  for idx, info in enumerate(group_info.itertuples(), 0):
    new_filename = info.sample_side + '_40x.tif'
    image_path = [image_path for image_path in image_paths if Path(image_path).name == new_filename]
    r, template_image = cv2.imreadmulti(str(image_path[0]), [], cv2.IMREAD_UNCHANGED)

    match_images = group_info[group_info.index > group_info.iloc[idx].name] # this prevents redoing comparisons
    for match_i in range(0, match_images.shape[0]):
      match_filename = match_images.iloc[match_i].sample_side + '_40x.tif'
      match_path = [image_path for image_path in image_paths if Path(image_path).name == match_filename]
      r, match_image = cv2.imreadmulti(str(match_path[0]), [], cv2.IMREAD_UNCHANGED)

      conf = compare_corners(template_image[0], match_image[0])
      names_l.append(str(new_filename + ' vs ' + match_filename))
      conf_l.append(conf)
results_df = pd.DataFrame({'names': names_l, 'confidence': conf_l})
results_df.to_csv((data_dir / 'output' / 'duplicate_comparison.csv').resolve())
