from cellpose import models, io
from pathlib import Path
import pandas as pd

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
data_dir = (Path.cwd().parent.parent / 'data' / 'GM130').resolve()

info_csv_path = Path(data_dir / 'GM130_image_log.csv')
sample_info = pd.read_csv(info_csv_path)

# get all the stacks we created from stack_channel_split.ijm
image_directory = (data_dir / 'images').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# find the cellpose models that Charlie helped train
golgi_model = (data_dir / 'cellpose' / 'models' / 'golgi_lab40X_GM5').resolve() # TODO rename these
nuclei_model = (data_dir / 'cellpose' / 'models' / 'Charlie_2D_HumanLoop10').resolve()
# where to save the cellpose output
nuclei_directory = (data_dir / 'cellpose' / 'nuclei').resolve()
golgi_directory = (data_dir / 'cellpose' / 'golgi').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True

for image in images:
  filename = Path(image).name
  sample = sample_info.loc[sample_info['old_filename'] == filename]
  if sample.size != 0:
    print('analyzing sample', filename)
    type = sample['channel'].values[0]
    if type == 'golgi':
      diam = 15
      model_path = golgi_model
      output_directory = golgi_directory
    elif type == 'nuclei':
      diam = 20
      model_path = nuclei_model
      output_directory = nuclei_directory
    else:
      diam = None
      model_type = 'cyto2'

    # DEFINE CELLPOSE MODEL
    model = models.CellposeModel(gpu=cuda, pretrained_model=str(model_path))
    chan = [0, 0]
    anis = sample['z_space'].values[0]
    img = io.imread(image)
    masks, flows, styles = model.eval(img, diameter=diam, channels=chan, do_3D=True, anisotropy=anis)

    # save results so you can load in gui
    # io.masks_flows_to_seg(img, masks, flows, diam, filename, chan)
    print('saving sample', filename)
    # save results as tiff
    outname = filename + '_cp.tiff'
    io.save_masks(img, masks, flows, filename, png=False, tif=True, savedir=str(output_directory), save_txt=False)
    # io.save_rois(masks, filename) # this doesn't appear to work ## imagej doesn't have 3D rois...
    print('\n')