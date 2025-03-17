import pandas as pd
import skimage.io
from scipy.ndimage import distance_transform_edt, label
from skimage.segmentation import watershed

import matplotlib.pyplot as plt

from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### MAIN ###
# First find all needed directories and load them all up
data_dir = (Path.cwd().parent.parent.parent / 'data' / 'Branches_activation').resolve()

sample_img_directory = (data_dir / 'analysis' / 'raw' / 'Fluor' / 'Nuc').resolve()
sample_mask_directory = (data_dir / 'analysis' / 'cellpose' / 'Fluor' / 'nuc_mask').resolve()
output_dir = (data_dir / 'analysis' / 'cellpose' / 'Fluor' / 'expanded_mask').resolve()


## Get sample info
# I use a csv file with info on all the images to match them to the samples and if they are nuclei or activation images
sample_info = pd.read_csv((data_dir / 'fluor_samples.csv').resolve())


## find all the orignal images, cellpose mask files, and rois
sample_dirs, img_paths = find_all_filepaths(sample_img_directory, '.tif')
mask_dirs, mask_paths = find_all_filepaths(sample_mask_directory, '.tif')


for sample in sample_info['filename'].unique():
  # match the sample name to the correct images and roi file
  mask_path = [mp for mp in mask_paths if Path(mp).name[:-13] == sample[:-4]]
  img_path = [sp for sp in img_paths if Path(sp).name[:-4] == sample[:-4]]

  if len(mask_path) != 0:
    # load the images
    img = skimage.io.imread(str(img_path[0]))
    nuc_mask = skimage.io.imread(str(mask_path[0]))

    # https://bioimagebook.github.io/chapters/2-processing/6-transforms/transforms.html
    im_dist_inv = distance_transform_edt(nuc_mask == 0)
    bw_dilated = (im_dist_inv < 15) # expand the cells by 15 pixels from nuclei
    cell_mask = watershed(im_dist_inv, nuc_mask, mask=bw_dilated, watershed_line=True) # resegment them based on cellpose nucs

    print('saving ' + sample)
    skimage.io.imsave(str(output_dir) + '\\' + sample[:-4] + '_cell_mask.tif', cell_mask)