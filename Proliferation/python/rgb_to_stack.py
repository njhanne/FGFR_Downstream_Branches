# Our lab fluorescence microscope takes RGB images, this converts them to 8bit and then stacks them
# It's been edited / shortened to only handle the pHH3 images for Branches project

import pandas as pd
import numpy as np
import cv2
from tifffile import imwrite

import re
from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


# Image Manipulation Helpers
def manipulate_image(image, info):
  # inverts the pHH3 images so stained cells appear bright instead of dark before converting to 8bit
  if info.invert == 'yes':
    image = cv2.bitwise_not(image)
  # convert RGB to 8bit
  if len(image.shape) == 3:
    image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
  return image


# End Image Manipulation Helpers

### MAIN ###
# start by getting all the raw images
data_dir = (Path.cwd().parent.parent / 'data' / 'proliferation').resolve()
input_directory = (data_dir / 'raw_images').resolve()
output_directory = (data_dir / 'new_images').resolve()
image_dirs, image_paths = find_all_filepaths(input_directory, '.tif')


## Get sample info
# I use a csv file with info on all the images to match them to the samples and if they are nuclei or phh3 images
image_info = pd.read_csv((data_dir / 'pHH3_image_log.csv').resolve())
cols = ['sample', 'side', 'section']
image_info['sample_side'] = image_info[cols].astype(str).apply('_'.join, axis=1) # combine sample and side info into one column

# loop through images, create the stacks, and save them
for image_pair_name in image_info['sample_side'].unique():
  skip = False
  # match sample to the manual info from your .csv
  manual_info = image_info.loc[image_info['sample_side'] == image_pair_name] # find all the rgb images that match the sample&side combo
  stack_size = manual_info.shape[0] # counts number of rgb images that will makeup stack. Should always be 2 for this project

  # optional to only get certain group, comment this out if you want all the images to run
  # if manual_info['group'].iloc[0] != 'DMSO':
  #   skip = True

  if not skip:
    # create the stack
    for idx, info in enumerate(manual_info.itertuples(),0):
      image_name = info.old_filename
      image_path = [image_path for image_path in image_paths if Path(image_path).name == image_name]
      image = cv2.imread(str(image_path[0])) # read rgb image
      if len(image.shape) == 3 & image.shape[2] == 3:
        image = image[..., ::-1]  # convert from BGR to RGB... this is a weird opencv convention

      if idx == 0:
        stack = np.empty([image.shape[0], image.shape[1], stack_size]) # create the stack (empty array)
      if bool(re.search('\d{1,}a.tif', image_name)): # look for the image ending in 'a', this is the pHH3 image
        stack[:,:,1] = manipulate_image(image, info) # add it as second image in stack
      else:
        stack[:,:,0] = manipulate_image(image, info) # add nuclei image as first image in stack

    # find image magnification from old filename
    magnification = re.search('\d{1,}x', manual_info['old_filename'].iloc[0]) # gets magnification
    if magnification is None:
      magnification = ''
    else:
      magnification = magnification[0]

    new_filename = manual_info['sample_side'].iloc[0] + '_' + magnification
    savename = (output_directory / (new_filename + '.tif')).resolve()
    # OME-TIFF should be TZCYX (frustrating) so we need to transpose
    stack = np.moveaxis(stack, -1, 0) # moves the last (-1) dimension to the first dimension (0) so instead of yxz it will be cyx
    imwrite(savename, stack.astype('uint8'))