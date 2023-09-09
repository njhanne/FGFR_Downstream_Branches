import cv2
import numpy as np
import pandas as pd

from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### MAIN ###
# First find all needed directories and load them all up
data_dir = (Path.cwd().parent.parent / 'data' / 'proliferation').resolve()
nuclei_directory = (data_dir / 'cellpose' / 'nuclei').resolve()
phh3_directory = (data_dir / 'cellpose' / 'pHH3').resolve()

## Get sample info
# I use a csv file with info on all the images to match them to the samples and if they are nuclei or phh3 images
sample_info = pd.read_csv((data_dir / 'pHH3_image_log.csv').resolve())
cols = ['sample', 'side', 'section']
sample_info['sample_side'] = sample_info[cols].astype(str).apply('_'.join, axis=1)

## find all the cellpose mask files created in 'run_cellpose_proliferation.py'
nuclei_dirs, nuclei_mask_paths = find_all_filepaths(nuclei_directory, '.png')
nuclei_mask_paths = [mp for mp in nuclei_mask_paths if Path(mp).name.endswith('_masks.png')]
phh3_dirs, phh3_mask_paths = find_all_filepaths(phh3_directory, '.png')
phh3_mask_paths = [mp for mp in phh3_mask_paths if Path(mp).name.endswith('_masks.png')]

## Cellpose created masks of the detected elements. For example, in the nuclei stained images, it would detect the nuclei
## then assign an index to each nucleus. The mask image pixel values on a nucleus are its index, so the 26th nucleus in
## an image will all have pixel value = 26 rather than a value corresponding to intensity. This allows us to just look at
## the largest pixel value to get a count of how many nuclei (or pHH3+ cells) were detected in each image.

## we're going to build the results in a loop, then convert to a dataframe for saving
results = []

for sample in sample_info['sample_side'].unique():
  nuc_img_path = [mp for mp in nuclei_mask_paths if Path(mp).name.startswith(sample)]
  stain_img_path = [mp for mp in phh3_mask_paths if Path(mp).name.startswith(sample)]
  nuc_img = cv2.imread(str(nuc_img_path[0]), -1) # https://docs.opencv.org/3.4/d8/d6a/group__imgcodecs__flags.html
  stain_img = cv2.imread(str(stain_img_path[0]), -1)

  row_indices = sample_info.loc[sample_info['sample_side'] == sample]
  result = []
  result.append('; '.join([str for str in row_indices.old_filename])) # another 'wonderful' one liner
  result.append(Path(nuc_img_path[0]).name)
  result.append(row_indices.iloc[0]['sample'])
  result.append(row_indices.iloc[0]['side'])
  result.append(row_indices.iloc[0]['section'])
  result.append(row_indices.iloc[0]['group'])
  result.append(np.amax(nuc_img))
  result.append(np.amax(stain_img))

  results.append(result)

results_df = pd.DataFrame(results, columns = ['old_filenames', 'new_filename', 'sample', 'side', 'section', 'group', 'nuclei_count', 'phh3_count'])

## save
results_df.to_csv((data_dir / 'results.csv').resolve())
