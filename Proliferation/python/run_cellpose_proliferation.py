from cellpose import models, io

from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
data_dir = (Path.cwd().parent.parent / 'data' / 'proliferation').resolve()
# get all the stacks we created from rgb_to_stack.py
image_directory = (data_dir / 'new_images2').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# find the cellpose models that Charlie helped train
phh3_model = (data_dir / 'cellpose' / 'train' / 'pHH3_lab40x_handrefined').resolve()
nuclei_model = (data_dir / 'cellpose' / 'train' / 'nuclei_lab40x_handrefined').resolve()
# where to save the cellpose output
nuclei_directory = (data_dir / 'cellpose' / 'nuclei').resolve()
phh3_directory = (data_dir / 'cellpose' / 'pHH3').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True

nuc_model = models.CellposeModel(gpu=cuda, pretrained_model=str(nuclei_model))
phh3_model = models.CellposeModel(gpu=cuda, pretrained_model=str(phh3_model))

for image_path in images:
  filename = Path(image_path).name
  print('analyzing image: ', filename)
  image = io.imread(str(image_path))

  nuc_masks, nuc_flows, nuc_styles = nuc_model.eval(image, diameter=44, channels=[1,0])
  io.save_masks(image, nuc_masks, nuc_flows, filename, png=True, channels=[1,0], savedir=str(nuclei_directory), save_txt=False)

  phh3_masks, phh3_flows, phh3styles = phh3_model.eval(image, diameter=44, channels=[2,0])
  io.save_masks(image, phh3_masks, phh3_flows, filename, png=True, channels=[2,0], savedir=str(phh3_directory), save_txt=False)
