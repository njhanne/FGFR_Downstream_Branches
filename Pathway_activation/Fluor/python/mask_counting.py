import numpy as np
import pandas as pd

import skimage.measure
import skimage.io

import matplotlib.pyplot as plt
# for mac
# matplotlib.use('TkAgg')

from PIL import Image, ImageDraw
from shapely.geometry import LineString
from shapely.plotting import plot_line, plot_points, plot_polygon
from read_roi import read_roi_file # pip install read-roi

from pathlib import Path
import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


# ImageJ-Shapely ROI helpers
def roi_to_shapely(roi):
  # https://shapely.readthedocs.io/en/stable/manual.html
  roi = list(roi.values())[0]
  # list comp to combine the x and y dictionary into list((x1,y1) (x2,y2) etc)   - then convert to shapely line
  roi_line = LineString([(i[0], i[1]) for i in zip(roi['x'], roi['y'])])
  # add a buffer to one side of the line to restrict roi distance from ectoderm edge
  # positive distance is 'left hand side' when single sided is set to True
  # pixels / um in these images: 3.296 for lab 20x
  # so what this will do is createe a polygon with rounded caps on both sides of my ROI line
  # then the buffer_neg is the 'left hand side' which should be the part outside of the tissue
  # the negative is then removed from the overall, so you are left with a polygon excluding the ectoderm, but
  # encompassing the mesenchyme with the big rounded caps left in
  # making the negative side bigger gets rid of some of the 'bowtie' issues
  buffer_full = roi_line.buffer(200 * 3.296)
  buffer_neg = roi_line.buffer(250 * 3.296, single_sided=True)
  clean_buffer = buffer_full.difference(buffer_neg)

  if clean_buffer.geom_type != 'Polygon':
    print('roi line has bowtie')
    # test generate the image
    fig, ax = plt.subplots()
    ax.imshow(img)
    plot_polygon(clean_buffer, ax)
    plot_line(roi_line, ax)
    plt.show()

  return clean_buffer


def shapely_to_mask(clean_buffer):
  mask_img = Image.new('1', (2048, 2048), 0)  # makes a blank image w/ same dimensions as image

  # Get pixel values of freehand shape
  if clean_buffer.geom_type != 'Polygon':
    print('roi line has bowtie')
  else:
    if clean_buffer.boundary.geom_type != 'LineString':
      print('roi clean buffer not  clean, trying to fix')
      clean_buffer = clean_buffer.boundary.geoms[0]
      xvals = clean_buffer.xy[0]
      yvals = clean_buffer.xy[1]
    else:
      xvals = clean_buffer.boundary.xy[0]
      yvals = clean_buffer.boundary.xy[1]
    polygon_coords = np.stack((xvals,yvals), axis=1).flatten()

    # draw the polygon onto the empty mask image to generate mask
    test = ImageDraw.Draw(mask_img) # this draws the mask onto the mask_image - it automatically writes
    test.polygon(polygon_coords.tolist(), outline=1, fill=1)  # draws freehand shape and fills with 1's

  return np.array(mask_img)

# End ImageJ ROI helpers

### MAIN ###
# First find all needed directories and load them all up
data_dir = (Path.cwd().parent.parent.parent / 'data' / 'Branches_activation').resolve()

sample_img_directory = (data_dir / 'analysis' / 'raw' / 'Fluor' / 'activation').resolve()
sample_mask_directory = (data_dir / 'analysis' / 'cellpose' / 'Fluor' / 'expanded_mask').resolve()
sample_roi_directory = (data_dir / 'analysis' / 'raw' / 'Fluor' / 'Nuc').resolve()
output_dir = (data_dir / 'analysis' / 'output').resolve()


## Get sample info
sample_info = pd.read_csv((data_dir / 'fluor_samples.csv').resolve())


## find all the orignal images, cellpose mask files, and rois
sample_dirs, img_paths = find_all_filepaths(sample_img_directory, '.tif')
mask_dirs, mask_paths = find_all_filepaths(sample_mask_directory, '.tif')
roi_dirs, roi_paths = find_all_filepaths(sample_roi_directory, '.roi')

## we're going to build the results in a loop, then convert to a dataframe for saving
results = []


for sample in sample_info['sample_name'].unique():
  # match the sample name to the correct images and roi file
  row_indices = sample_info.loc[sample_info['sample_name'] == sample]
  nuc_filename = row_indices.loc[row_indices['channel'] == 'DAPI']['filename'].values[0]
  activation_row = row_indices.loc[row_indices['channel'] == 'activation']

  mask_path = [mp for mp in mask_paths if Path(mp).name.startswith(nuc_filename[:-4])]
  img_path = [sp for sp in img_paths if Path(sp).name.startswith(activation_row['filename'].values[0][:-4])]
  roi_path =[rp for rp in roi_paths if Path(rp).name.startswith(nuc_filename[:-4])]

  if len(mask_path) != 0:
    # load the images
    img = skimage.io.imread(str(img_path[0]))
    cell_mask = skimage.io.imread(str(mask_path[0]))

    print(sample)

    roi = read_roi_file(str(roi_path[0]))
    mask_polygon = roi_to_shapely(roi)
    filter_mask = shapely_to_mask(mask_polygon)

    label_props = skimage.measure.regionprops(cell_mask, img)

    # calculate the mean intensity for each cell
    means = np.asarray([d['intensity_mean'] for d in label_props])
    # calculate thresh and filter dim cells
    filter_thresh = activation_row.iloc[:,-6:].mean(axis=1).values[0] # average the 6 good/bad cells
    n_pos_cells = (means > filter_thresh).sum() # positive cells. sum() only counts True value
    n_total_cells = means.size

    # filter out cells outside the masked region

    labels = np.asarray([d['label'] for d in label_props])
    centroids = np.asarray([d['centroid'] for d in label_props])
    # rounds the centroids then converts to int
    centroids = np.round(centroids).astype(int)
    # sees if the centroids are inside the masked area
    masked_centroids = filter_mask[centroids[:,0],centroids[:,1]]
    # matches the masked_centroids to their label_id
    masked_cells = labels[masked_centroids]
    # deletes the labels from the image that are outside the masked region
    filtered_image = np.isin(cell_mask, masked_cells)*cell_mask

    # save image of filtered cells
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.imshow(filtered_image)
    save_filename = sample + '_masked_cells.png'
    fig.savefig(output_dir / save_filename)
    plt.close(fig)

    # calculate number of masked and filtered cells
    masked_means = means[masked_centroids]
    n_pos_masked_cells = (masked_means > filter_thresh).sum() # positive cells. sum() only counts True value
    n_total_masked_cells = masked_means.size

    # get all the info into the results table
    result = []
    result.append(activation_row.iloc[0]['filename'])
    result.append(activation_row.iloc[0]['treatment'])
    result.append(activation_row.iloc[0]['pathway'])
    result.append(activation_row.iloc[0]['sample_id'])
    result.append(activation_row.iloc[0]['section'])
    result.append(activation_row.iloc[0]['sample_name'])

    # cell counts
    result.append(filter_thresh)
    result.append(n_total_cells)
    result.append(n_pos_cells)
    result.append(n_total_masked_cells)
    result.append(n_pos_masked_cells)

    results.append(result)

results_df = pd.DataFrame(results, columns = ['old_filename', 'treatment', 'pathway', 'sample_id', 'section', 'sample_name', 'thresh', 'total_cells', 'pos_cells', 'total_masked_cells', 'pos_masked_cells'])

## save
results_df.to_csv((data_dir / 'fluor_results.csv').resolve())