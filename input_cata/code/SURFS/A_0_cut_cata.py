# @Author: lshuns
# @Date:   2021-04-15, 14:24:13
# @Last modified by:   lshuns
# @Last modified time: 2021-04-19, 14:46:18

### cut the whole catalogue into several smaller catalogues
###     cut along the RA direction, dec range is not changed

import numpy as np

from astropy.io import fits
from astropy.table import Table

# ++++++++++++++ I/O

# original cata
infile_surfs = '/disks/shear10/ssli/ImSim/input/SURFS_cata/ori/SURFS_SHARK_LC_100sqdeg_testcols.fits'

# output
outdir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'

# ++++++++++++++ general settings
N_pieces = 4

# ++++++++++++++ workhorse
with fits.open(infile_surfs) as hdul:
    surfs_cata = Table(hdul[1].data)
print('surfs ori', len(surfs_cata))

### index
surfs_cata['index'] = np.arange(len(surfs_cata)).astype(int)
outpath = os.path.join(outdir, 'SURFS_SHARK_LC_all.fits')
surfs_cata.write(outpath, format='fits')
print(f'full cata with index saved as {outpath}')

### split to equally pieces
ra_min = np.min(surfs_cata['ra'])
ra_max = np.max(surfs_cata['ra'])
dec_min = np.min(surfs_cata['dec'])
dec_max = np.max(surfs_cata['dec'])
print(f'whole range [ra], [dec]: [{ra_min}, {ra_max}], [{dec_min}, {dec_max}]')
ra_ranges = np.around(np.linspace(ra_min, ra_max, N_pieces+1), decimals=1)
ra_ranges[0] = ra_min - 0.01 # for edges
ra_ranges[-1] = ra_max + 0.01 # for edges
print(f'ra boundaries: {ra_ranges}')
for i_patch in range(N_pieces):
    ra_min = ra_ranges[i_patch]
    ra_max = ra_ranges[i_patch+1]
    mask_tmp = (surfs_cata['ra']>=ra_min) & (surfs_cata['ra']<ra_max)
    surfs_cata_selec = surfs_cata[mask_tmp]
    ra_min = np.min(surfs_cata_selec['ra'])
    ra_max = np.max(surfs_cata_selec['ra'])
    dec_min = np.min(surfs_cata_selec['dec'])
    dec_max = np.max(surfs_cata_selec['dec'])
    print(f'Patch {i_patch}', len(surfs_cata_selec))
    print(f'=======> [{ra_min}, {ra_max}], [{dec_min}, {dec_max}]')

    outpath = os.path.join(outdir, f'SURFS_SHARK_LC_patch{i_patch}.fits')
    surfs_cata_selec.write(outpath, format='fits')
    print(f'cata saved as {outpath}')

# surfs ori per sq2 711976.2129629629
# whole range [ra], [dec]: [211.5, 223.5], [-4.5, 4.5]
# ra boundaries: [211.49 214.5  217.5  220.5  223.51]
# Patch 0 per sq2 717258.1111111111
# =======> [211.5, 214.49998474121094], [-4.5, 4.499999523162842]
# Patch 1 per sq2 708026.1851851852
# =======> [214.5, 217.49998474121094], [-4.5, 4.5]
# Patch 2 per sq2 705157.1111111111
# =======> [217.5, 220.49998474121094], [-4.5, 4.499999046325684]
# Patch 3 per sq2 717463.4444444445
# =======> [220.5, 223.5], [-4.5, 4.5]
