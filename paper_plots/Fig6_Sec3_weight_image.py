# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-01-05 16:09:15
# @Last Modified by:   lshuns
# @Last Modified time: 2022-01-05 16:55:31

### transfer fits file to pdf image

import matplotlib.pyplot as plt
from astropy.io import fits

import matplotlib.gridspec as gridspec

# >>>>> read fits file

## simulation
image_file = '/disks/shear16/ssli/ImSim/output/skills_v07D4/part0/m283m283/images/THELI/tile137.0_-1.5_bandr_rot90.weight.fits'
image_sim = fits.getdata(image_file, ext=0)
#### get center
cropx = 18300
cropy = 19700 
y,x = image_sim.shape
startx = x//2-(cropx//2)
starty = y//2-(cropy//2)    
image_sim = image_sim[starty:starty+cropy,startx:startx+cropx]

## real data
image_file = '/disks/shear10/ssli/THELI/KiDS_DR4.0_133.0_0.5_r_det_wei.fits'
image_data = fits.getdata(image_file, ext=0)
#### get center
cropx = 18300
cropy = 19700 
y,x = image_data.shape
startx = x//2-(cropx//2)
starty = y//2-(cropy//2)    
image_data = image_data[starty:starty+cropy,startx:startx+cropx]

## collect
images = [image_sim, image_data]
titles = ['SKiLLS', 'KiDS tile 133.0\_0.5']

# save
outpath = './plots/weight_ima_comp.pdf'

# >>>>> plot

# font size
plt.rc('font', size=14)
# tex
plt.rcParams["text.usetex"] = True

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

fig, axs = plt.subplots(1, 2, constrained_layout=True, figsize=(10, 6))

for i_plot in range(2):

    ax = axs[i_plot]
    image = images[i_plot]
    title = titles[i_plot]

    ax.imshow(image, cmap='gray')
    ax.set_title(title)

    ax.axis('off')

# save
if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("plot saved as", outpath)
