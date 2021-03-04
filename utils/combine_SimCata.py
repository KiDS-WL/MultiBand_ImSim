# @Author: lshuns
# @Date:   2021-03-03, 18:21:04
# @Last modified by:   lshuns
# @Last modified time: 2021-03-03, 19:54:58

### a simple script to combine catalogues produced by the main pipeline
###     it can be used to combine catalogues from different running_tags
###             and assgin photo-z based on input id of galaxies
###                 in which case at least one running_tag should contain photo-z
###                     and same galaxies are simulated in different running_tags

import os
import glob

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ the only variables you want to modify
# where to save the final catalogue
### supported formats: fits, csv, feather
###             guess from the file suffix
out_path = ''

# where to find the input catalogues
in_dir = ''

# a list of running tags being combined
run_tag_list = []

# assgin photo-z ?
assgin_photoz = True
# the run_tag contains the photo-z
run_tag_withz = ''

# +++++++++++++++++++++++++++++ workhorse
cata_final = []
for run_tag in run_tag_list:

    # basic info for input shear values
    opened_file = open(os.path.join(in_dir, run_tag, 'basic_info.txt'), 'r')
    all_lines = opened_file.readlines()
    opened_file.close()
    useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
    try:
        g1 = float(useful_line.split()[2])
    except ValueError:
        g1 = float(useful_line.split()[2][:-1])
    g2 = float(useful_line.split()[3])

    # main catalogues from simulation
    cata_tag = []
    file_list = glob.glob(os.path.join(in_dir, run_tag, 'catalogues', '*_combined.*'))
    if not file_list:
        raise Exception('Cannot find any catalogues from the main pipeline!\n\
        make sure the main pipeline ended properly to the end (taskID=7 is required)!')
    for file in file_list:
        file_type = file[-3:]
        if file_type == 'csv':
            cata = pd.read_csv(file)
        elif file_type == 'her':
            cata = pd.read_feather(file)
        elif file_type == 'its':
            with fits.open(file) as hdul:
                cata = Table(hdul[1].data).to_pandas()
        else:
            raise Exception(f'Not supported input file type! {file}')
        cata_tag.append(cata)
    cata_tag = pd.concat(cata_tag)
    if assgin_photoz:
        if run_tag == run_tag_withz:
            cata_tag_withz = cata_tag[['id_input', 'Z_B']]
        else:
            cata_tag = cata_tag.merge(cata_tag_withz, on='id_input', how='left')

    cata_final.append(cata_tag)

cata_final = pd.concat(cata_final)
# dummy values for nan
cata_final.fillna(-999, inplace=True)

# save
file_type = out_path[-3:]
if file_type == 'her':
    data_final.to_feather(out_path)
elif file_type == 'csv':
    data_final.to_csv(out_path, index=False)
elif file_type == 'fits':
    Table.from_pandas(data_final).write(out_path, format='fits')
else:
    raise Exception(f'Not supported output file type! {out_path}')

print(f'Combined catalogue saved as {out_path}')
