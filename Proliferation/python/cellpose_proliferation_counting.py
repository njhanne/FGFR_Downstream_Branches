# from cellpose import plot

import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import os
import re

import tkinter as tk
from tkinter import filedialog

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths

# File Loading and DataFrame Results Helpers
def initialize_results_dataframe(input_csv, length):
  df = pd.DataFrame(columns=input_csv.columns)
  df.drop(columns =['invert'])
  df['nuclei_count']=''
  df['stain_count']=''
  return df


def get_all_vals(d):
  # gets all the values in nested dictionary (d)
  vals = []
  for k, v in d.items():
    if isinstance(v, dict):
      vals.extend(get_all_vals(v))
    else:
      vals.append(v)
  return vals


# End DataFrame Helpers

### MAIN ###
# First find all needed directories and load them all up
root = tk.Tk()
root.withdraw()

nuclei_directory = filedialog.askdirectory(title='Select masked nuclei directory')
stain_directory = filedialog.askdirectory(title='Select other stain masked directory')
info_csv_path = filedialog.askopenfilename(title='Select sample log')
sample_info = pd.read_csv(info_csv_path)

nuclei_dirs, nuclei_mask_paths = find_all_filepaths(nuclei_directory, '.png')
nuclei_mask_paths = [mp for mp in nuclei_mask_paths if mp.endswith('_masks.png')]
stain_dirs, stain_mask_paths = find_all_filepaths(stain_directory, '.png')
stain_mask_paths = [mp for mp in stain_mask_paths if mp.endswith('_masks.png')]

# results_df = initialize_results_dataframe(sample_info, len(stain_mask_paths.count))
results = []

cols = ['sample', 'side']
sample_info['sample_side'] = sample_info[cols].astype(str).apply('_'.join, axis=1)
for sample in sample_info['sample_side'].unique():
  nuc_img_path = [mp for mp in nuclei_mask_paths if os.path.basename(mp).startswith(sample)]
  stain_img_path = [mp for mp in stain_mask_paths if os.path.basename(mp).startswith(sample)]
  nuc_img = cv2.imread(nuc_img_path[0], -1)
  stain_img = cv2.imread(stain_img_path[0], -1)

  row_indices = sample_info.loc[sample_info['sample_side'] == sample]
  result = []
  result.append('; '.join([str for str in row_indices.old_filename])) # another 'wonderful' one liner
  result.append(os.path.basename(nuc_img_path[0]))
  result.append(row_indices.iloc[0]['sample'])
  result.append(row_indices.iloc[0]['side'])
  result.append(row_indices.iloc[0]['group'])
  result.append(np.amax(nuc_img))
  result.append(np.amax(stain_img))

  results.append(result)

results_df = pd.DataFrame(results, columns = ['old_filenames', 'new_filename', 'sample', 'side', 'group', 'nuclei_count', 'stain_count'])

results_df.to_csv(os.path.normpath(os.path.dirname(info_csv_path) + '\\' + 'results.csv'))
print('hey')
