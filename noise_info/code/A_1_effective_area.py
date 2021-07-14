# @Author: lshuns
# @Date:   2021-04-26, 20:52:22
# @Last modified by:   lshuns
# @Last modified time: 2021-05-12, 16:30:56

### obtain the effective area for each tile

import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path

# ++++++ I/O
inDir = '/disks/shear15/KiDS/KiDS-1000/ESO-DR4-9band-masks/'
pathlist = Path(inDir).glob('KiDS_*_ugriZYJHKs_msk.fits.gz')

outpath = '../kids_dr4_effective_pixels.csv'

# ++++++ workhorse
outData = pd.DataFrame({'label':[], 'pixels_28668':[], 'pixels_32764':[]})
for path in pathlist:
    # because path is object not string
    path_in_str = str(path)
    label = path_in_str[61:-23]
    print('file', path_in_str)
    print('label', label)
    with fits.open(path_in_str) as hdul:
        data = hdul[0].data

    # 9-band no AW-r-band masks (used to create the WL mosaic catalogues)
    nine_band_mask = (data & 28668)==0
    data_tmp = data[nine_band_mask]
    N_pixel_28668 = len(data_tmp)
    print('Number (28668)', N_pixel_28668)

    # 9-band masks
    nine_band_mask = (data & 32764)==0
    data_tmp = data[nine_band_mask]
    N_pixel_32764 = len(data_tmp)
    print('Number (32764)', N_pixel_32764)

    outData.loc[len(outData)] = [label, N_pixel_28668, N_pixel_32764]

outData = outData.astype({'pixels_28668': 'int', 'pixels_32764': 'int'})
outData.to_csv(outpath, index=False)
print(f'saved as {outpath}')
