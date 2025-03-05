from cellpose import models, io

from pathlib import Path

import re

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
data_dir = (Path.cwd().parent.parent / 'data' / 'Branches_activation' / 'analysis').resolve()
# get all the stacks we created from rgb_to_stack.py
image_directory = (data_dir / 'cellpose' / 'deconvolved_filtered').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# find the cellpose models that Charlie helped train
dab_model = (data_dir / 'cellpose' / 'train' / 'DAB_bg_deconvolved_HL1').resolve()
# where to save the cellpose output
labelimg_directory = (data_dir / 'cellpose' / 'DAB').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True

dab_model = models.CellposeModel(gpu=cuda, pretrained_model=str(dab_model))


# I did some hand refining of the cyto3 model that comes with Cellpose v3. Even with hand training it does a poor job
# with the DAB staining due to how faint some of the staining is. Since we want 100% of the cells to be tracked regardless
# of how faint the staining is, I further processed the images in ImageJ so they are easier for cellpose to detect
# First I did a color deconvolution to separate the dab - this is used for all the analyses downstream of this step, too
# Next I performed a background subtraction and a 'unsharp' filter which increases local contrast for the dim cells.

for image_path in images:
  filename = Path(image_path).name
  print('analyzing image: ', filename)
  image = io.imread(str(image_path))

  dab_masks, dab_flows, dab_styles = dab_model.eval(image, diameter=20)
  io.save_masks(image, dab_masks, dab_flows, filename, png=True, savedir=str(labelimg_directory), save_txt=False)

  # for this analysis we are also saving the roi for ImageJ so we can load them up and pick a few representative
  # 'positive' and 'negative' levels

  io.save_rois(dab_masks, str(labelimg_directory) + '/' + filename)
