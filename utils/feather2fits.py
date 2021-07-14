# @Author: lshuns
# @Date:   2021-06-07, 13:13:51
# @Last modified by:   lshuns
# @Last modified time: 2021-06-07, 13:39:35

### a simple script to transfer feather to fits
###     or vice versa

########################
# usage: feather2fits.py [-h] [--ori_cata ORI_CATA] [--out_cata OUT_CATA]
#                        [--transfer_type tansfer_type]
# feather2fits: transfer feather file to fits (or vice versa).
# optional arguments:
#   -h, --help            show this help message and exit
#   --ori_cata ORI_CATA   the original catalogue path.
#   --out_cata OUT_CATA   the output catalogue path.
#   --transfer_type tansfer_type
#                         transfer type:
#                             feather2fits: from feather file to fits file
#                             fits2feather: from fits file to feather file
########################

import pandas as pd

from astropy.io import fits
from astropy.table import Table

import os
import argparse

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"feather2fits: transfer feather file to fits (or vice versa).",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--ori_cata", type=str,
    help="the original catalogue path.")
parser.add_argument(
    "--out_cata", type=str,
    help="the output catalogue path.")
parser.add_argument(
    "--transfer_type", type=str, choices=['feather2fits', 'fits2feather'], metavar='tansfer_type',
    help="transfer type:\n\
    feather2fits: from feather file to fits file\n\
    fits2feather: from fits file to feather file")

## arg parser
args = parser.parse_args()
ori_cata = args.ori_cata
out_cata = args.out_cata
transfer_type = args.transfer_type

# +++++++++++++++++++++++++++++ workhorse

# from feather to fits
if transfer_type == 'feather2fits':
    cata = pd.read_feather(ori_cata)
    print('original catalogue loaded from', ori_cata)
    Table.from_pandas(cata).write(out_cata, format='fits')
    print('catalogue saved as', out_cata)

# from fits to feather
if transfer_type == 'fits2feather':
    with fits.open(ori_cata) as hdul:
        cata = Table(hdul[1].data).to_pandas()
    print('original catalogue loaded from', ori_cata)
    cata.to_feather(out_cata)
    print('catalogue saved as', out_cata)
