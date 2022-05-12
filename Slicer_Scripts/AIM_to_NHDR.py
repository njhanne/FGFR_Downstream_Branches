# This script will write a .nhdr (.nrrd header) file that points to your .aim file (from Scanco).
# This way you don't have to convert all the .aim files to .nrrd and have the volumes take up twice as much room.

# This is not fully automated, there are a number of steps you need to do yourself!!

# Make sure SlicerMorph and all it's dependencies are installed on Slicer

# You will need the .txt file from IPL (?) with all the header information. If it isn't there you may be able to get it
# from ImageJ. ### Save the header info as a .txt file with the same name as the .aim file. ###

import re, os
from pathlib import Path

def find_all_filepaths(directory, extension):
    filepaths = []
    for f in os.scandir(directory):
        if f.is_file():
            if os.path.splitext(f.name)[1].lower() in extension:
                filepaths.append(f.path)
    return filepaths

process_dir = r'C:\Users\njhan\Box\FGF_inhibitor_paper_5-26-2020\Morphology\3D_data\Raw_Scans\batch2'
aim_files = find_all_filepaths(process_dir, '.aim')
hx_files = find_all_filepaths(process_dir, '.hx')
txt_files = find_all_filepaths(process_dir, '.txt')

for aim in aim_files:
    ipl_txt_file = [file for file in txt_files if Path(file).stem == Path(aim).stem]
    ipl_txt = open(ipl_txt_file[0]).readlines()
    filename = os.path.basename(aim)
    savename_path = aim[0:-4] + '.nhdr'
    ### Things I know how to get from the IPL and HX files ###
    if len(hx_files) == len(aim_files):
        hx_txt_file = [file for file in hx_files if Path(file).stem == Path(aim).stem]
        hx_txt = open(hx_txt_file[0]).readlines()[1].split()
        endian = hx_txt[3]
        type = hx_txt[5]
        encoding = re.search(r'[a-z]+', hx_txt[1])[0]
    else:
        endian = 'little' # must be hardcoded / guessed
        type = next((s for s in ipl_txt if 'Type of data' in s), None).split()[3]
        encoding = 'raw' # hardcoded

    # why did I do this? It's CT data, it will always be 3D... I'm a moron
    # I'll leave this in as proof
    dim_line = next((s for s in ipl_txt if 'dim' in s), None).split()
    dim = len(dim_line) - 1
    size = []
    for i in range(dim):
        size.append(dim_line[i+1])

    offset = ipl_txt[0].split()[-1]

    space_origin = next((s for s in ipl_txt if 'off ' in s), None).split()[1:]
    space_origin = [str(float(s)) for s in space_origin]
    space_origin = ", ".join(space_origin)

    ### Things that need to be in nhdr file but I am not sure where to find without using the Slicer plugin ###
    space = "left-posterior-superior"
    space_directions = "(1.0, 0.0, 0.0) (0.0, 1.0, 0.0) (0.0, 0.0, 1.0)"
    kinds = 'domain domain domain'

    lines = ['NRRD0004',
             'type: ' + type,
             'dimension: ' + str(dim),
             'space: ' + space,
             'sizes: ' + size[0] + ' ' + size[1] + ' ' + size[2],
             'space directions: ' + space_directions,
             'kinds: ' + kinds,
             'endian: ' + endian,
             'encoding: ' + encoding,
             'space origin: ' + '(' + space_origin + ')',
             'byte skip: ' + offset,
             'data file: ' + filename]

    with open(savename_path, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
