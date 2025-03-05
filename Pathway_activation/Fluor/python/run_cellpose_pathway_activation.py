from cellpose import models, io

from pathlib import Path

import re

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
data_dir = (Path.cwd().parent.parent.parent / 'data' / 'Branches_activation' / 'analysis').resolve()
# get all the stacks we created from rgb_to_stack.py
image_directory = (data_dir / 'raw' / 'Fluor' / 'Nuc').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# where to save the cellpose output
labelimg_directory = (data_dir / 'cellpose' / 'Fluor' / 'nuc_mask').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True

nuc_model = models.CellposeModel(gpu=cuda, model_type='cyto3')


for image_path in images:
  filename = Path(image_path).name
  print('analyzing image: ', filename)
  image = io.imread(str(image_path))

  masks,flows, styles = nuc_model.eval(image, diameter=25)
  io.save_masks(image, masks, flows, filename, png=False, tif=True, savedir=str(labelimg_directory), save_txt=False)
