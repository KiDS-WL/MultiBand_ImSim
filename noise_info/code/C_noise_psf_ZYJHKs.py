# @Author: lshuns
# @Date:   2021-03-25, 17:52:23
# @Last modified by:   ssli
# @Last modified time: 2021-05-11, 15:42:36

### noise info for NIR bands

import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

# ++++++++++++++++++ I/O

# optical noise file
noise_file = '../test_noise_ugri.csv'

# output files
outfile = '../test_noise_NIR.csv'

# original NIR noise file
NIR_file = '/disks/shear10/ssli/KiDS_noise_seeing/frameDetails_1350.dat'

# measured moffat parameters
NIR_moffat_file = '/disks/shear10/ssli/KiDS_noise_seeing/moffat.NIR'

# ++++++++++++++++++ workhorse

bands = ['Z', 'Y', 'J', 'H', 'Ks']

# ====== noise info from state file
# RA Dec from optical tiles
data = pd.read_csv(noise_file)
radec = data['label'].str.split(r"_", expand=True)
radec.columns = ['ra', 'dec']
radec = pd.concat([data['label'], radec], axis=1)
radec = radec.astype({'ra': 'float', 'dec': 'float'})
## sort for fast searching
radec = radec.sort_values(by=['ra','dec'], ignore_index=True)

# for saving results
final_data = pd.DataFrame({'label': radec['label'].values})

# load original NIR info
data_NIR = pd.read_csv(NIR_file, sep=" ")
data_NIR = data_NIR.sort_values(by=['RAchip','DECchip'], ignore_index=True)
print("Original NIR total number", len(data_NIR))

# chip info for NIR camera
Npix_chip = int(np.ceil(0.2*3600./0.34))
# artificial tile info
Npix_tile = int(np.ceil(1.0*3600./0.34))

# loop over bands
for band in bands:

    data_tmp = data_NIR.loc[data_NIR['Filter']==band]
    print(f"Total number of chips in {band} band: {len(data_tmp)}")

    # loop over all the optical tiles
    for index, row in radec.iterrows():
        ra_min = row['ra'] - 0.5
        ra_max = row['ra'] + 0.5
        dec_min = row['dec'] - 0.5
        dec_max = row['dec'] + 0.5
        print('tile', row['label'])

        # chips within selected tile region
        data_tmp2 = data_tmp[(data_tmp.RAchip>ra_min)&(data_tmp.RAchip<ra_max)&(data_tmp.DECchip>dec_min)&(data_tmp.DECchip<dec_max)]
        N_chip = len(data_tmp2)
        print('Number of chips', N_chip)

        # minimum required chips
        ## 4 for J band , 2 for others
        if ((N_chip >= 4) and (band == 'J')) or ((N_chip >= 2) and (band != 'J')):

            # seeing from weighted average
            seeing = (data_tmp2['HeaderSeeing'].values/(data_tmp2['HeaderSkyNoise'].values**2)).sum()/(1./(data_tmp2['HeaderSkyNoise'].values**2)).sum()
            final_data.loc[final_data['label']==row['label'], f'HeaderSeeing_{band}'] = seeing

            seeing = (data_tmp2['MeasSeeing'].values/(data_tmp2['HeaderSkyNoise'].values**2)).sum()/(1./(data_tmp2['HeaderSkyNoise'].values**2)).sum()
            final_data.loc[final_data['label']==row['label'], f'MeasSeeing_{band}'] = seeing

            # noise from built "weight map"
            weight_map = np.zeros((Npix_tile, Npix_tile))
            for index2, row2 in data_tmp2.iterrows():
                # position
                index_row = int((row2['DECchip'] - 0.1 - dec_min)*3600./0.34)
                index_col = int((row2['RAchip'] - 0.1 - ra_min)*3600./0.34)
                weight_map[index_row:(index_row+Npix_chip), index_col:(index_col+Npix_chip)] += 1./(row2['HeaderSkyNoise']**2)

            final_data.loc[final_data['label']==row['label'], f'rms_{band}'] = 1./np.sqrt(np.nanmedian(weight_map[weight_map!=0]))
        else:
            print('Passed tile:', row['label'])

## only preserve those with all bands info
final_data.dropna(inplace=True)
print('Number of tiles with all bands noise', len(final_data))

# ====== Moffat beta
#### 1. inspect relation from representative tiles
moffat_data = pd.read_csv(NIR_moffat_file, sep=' ', names=['file', 'number', 'seeing', 'beta'])
moffat_data['seeing'] *= 0.34 # to arcsec
# colours = ['r', 'b', 'g', 'c', 'm']
# labels  =['H', 'J', 'Ks', 'Y', 'Z']
# for i in range(5):
#     CR = colours[i]
#     label = labels[i]
#     moffat_tmp = moffat_data[(7*i):(7*i+7)]
    # plt.plot(moffat_tmp['seeing'].values, moffat_tmp['beta'].values, marker='o', color=CR, linestyle=' ', markersize=2, label=label)
# plt.xlabel('seeing')
# plt.ylabel('beta')
# plt.legend()
# plt.show()

#### 2. fitting
x = moffat_data['seeing'].values
y = np.log(moffat_data['beta'].values)
def func2(x, a, b, c):
    return a * np.exp(-b * x) + c
popt, pcov = curve_fit(func2, x, y)
# xdata = np.linspace(np.min(x), np.max(x), 50)
# plt.plot(xdata, np.exp(func2(xdata, *popt)), 'k-', label='exp2')
# plt.xlabel('seeing')
# plt.ylabel('beta')
# plt.legend()
# plt.show()

#### 3. apply
for band in bands:
    xdata = np.array(final_data[f'MeasSeeing_{band}'])
    final_data.loc[:, f'MeasBeta_{band}'] = np.exp(func2(xdata, *popt))

# === save
final_data.to_csv(outfile, index=False)
print('final results saved as', outfile)
