import numpy as np
import pandas as pd
import skimage.measure
import skimage.io

# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt

from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### MAIN ###
# First find all needed directories and load them all up
data_dir = (Path.cwd().parent.parent.parent / 'data' / 'Branches_activation').resolve()

sample_img_directory = (data_dir / 'analysis' / 'raw' / 'Fluor' / 'activation').resolve()
sample_mask_directory = (data_dir / 'analysis' / 'cellpose' / 'Fluor' / 'expanded_mask').resolve()
output_dir = (data_dir / 'analysis').resolve()


## Get sample info
sample_info = pd.read_csv((data_dir / 'fluor_samples.csv').resolve())


## find all the orignal images, cellpose mask files, and rois
sample_dirs, img_paths = find_all_filepaths(sample_img_directory, '.tif')
mask_dirs, mask_paths = find_all_filepaths(sample_mask_directory, '.tif')


## we're going to build the results in a loop, then convert to a dataframe for saving
results = []


for sample in sample_info['sample_name'].unique():
  # match the sample name to the correct images and roi file
  row_indices = sample_info.loc[sample_info['sample_name'] == sample]
  nuc_filename = row_indices.loc[row_indices['channel'] == 'DAPI']['filename'].values[0]
  activation_row = row_indices.loc[row_indices['channel'] == 'activation']

  mask_path = [mp for mp in mask_paths if Path(mp).name.startswith(nuc_filename[:-4])]
  img_path = [sp for sp in img_paths if Path(sp).name.startswith(activation_row['filename'].values[0][:-4])]

  if len(mask_path) != 0:
    # load the images
    img = skimage.io.imread(str(img_path[0]))
    cell_mask = skimage.io.imread(str(mask_path[0]))

    print(sample)

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