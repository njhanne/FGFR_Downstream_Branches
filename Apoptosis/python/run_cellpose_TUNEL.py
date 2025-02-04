from cellpose import models, io

from pathlib import Path
import re

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
data_dir = (Path.cwd().parent.parent / 'data' / 'Apoptosis').resolve()
# get all the stacks we created from rgb_to_stack.py
image_directory = (data_dir / 'new_images').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# find the cellpose models that Charlie helped train
tunel_model = (data_dir / 'cellpose' / 'train' / 'pHH3_lab40x_handrefined').resolve()
nuclei_model = (data_dir / 'cellpose' / 'train' / 'nuclei_lab40x_handrefined').resolve()
# where to save the cellpose output
nuclei_directory = (data_dir / 'cellpose' / 'nuclei').resolve()
tunel_directory = (data_dir / 'cellpose' / 'TUNEL').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True

nuc_model = models.CellposeModel(gpu=cuda, pretrained_model=str(nuclei_model))
# tunel_model = models.CellposeModel(gpu=cuda, pretrained_model=str(tunel_model)) # this isn't working well compared to cyto3
tunel_model = models.CellposeModel(gpu=cuda, model_type='cyto3')

# we will use the same models trained for the pHH3, the nuclei are still nuclei shaped after all...
# however, not all these images were captured at 40x, so we need to change the scale of the diameter
# 10x 10, 20x 20, 40x 40
diameters = {'10x': 10,
             '20x': 15,
             '40x': 35}

for image_path in images:
  filename = Path(image_path).name
  print('analyzing image: ', filename)
  image = io.imread(str(image_path))

  mag = re.search('\d{1,}x', Path(image_path).stem)[0]

  nuc_masks, nuc_flows, nuc_styles = nuc_model.eval(image, diameter=diameters[mag], channels=[1,0])
  io.save_masks(image, nuc_masks, nuc_flows, filename, png=True, channels=[1,0], savedir=str(nuclei_directory), save_txt=False)

  tunel_masks, tunel_flows, tunel_styles = tunel_model.eval(image, diameter=diameters[mag], channels=[2,0], normalize={'percentile' : [10,100], 'normalize':True})
  io.save_masks(image, tunel_masks, tunel_flows, filename, png=True, channels=[2,0], savedir=str(tunel_directory), save_txt=False)
