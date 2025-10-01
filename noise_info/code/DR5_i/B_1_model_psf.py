# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-11-04 22:09:08
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-19 15:48:06

### model PSF with Moffat profile for each tile
###### need outside scripts save in:
######      /disks/shear10/ssli/KiDS_noise_seeing/code/PSFmodelling

import os
import glob
import shutil
import subprocess

# >>>>>> I/O
# the parent path for outputs
outDir = '/disks/shear15/ssli/KiDS_noise_seeing/'

# where to find the images
inDir = '/disks/shear15/KiDS_DR5/Awe/images/'

# where to find the extra scripts
codeDir = '/disks/shear10/ssli/KiDS_noise_seeing/code/PSFmodelling/'

# >>>>>> workhorse
for i_label in ('i1', 'i2'):

    # where to save the outputs
    outDir_sub = os.path.join(outDir, f'psf2moffat_{i_label}')
    if os.path.exists(outDir_sub):
        shutil.rmtree(outDir_sub)
    os.mkdir(outDir_sub)

    # copy all extra scripts
    for file in glob.glob(os.path.join(codeDir, '*')):
        shutil.copy2(file, os.path.join(outDir_sub, os.path.basename(file)))

    # find all images
    inpath_list = glob.glob(os.path.join(inDir, f'KIDS_*_{i_label}.fits'))
    print(i_label, 'number of images found', len(inpath_list))
    for inpath in inpath_list:

        # find stars
        psfcat = os.path.join(outDir_sub, os.path.basename(inpath) + '.psf.cat')
        if os.path.exists(psfcat):
            print(f'{psfcat} already exists')
        else:
            code_star = os.path.join(outDir_sub, f'mkpsfcat_kids{i_label}.csh')
            proc = subprocess.run(code_star + f' {inpath} {psfcat}',
                    shell=True, cwd=outDir_sub)
            print(f'{psfcat} saved')

        # psf modelling
        psfmoffat = psfcat + '.moffat'
        if os.path.exists(psfmoffat):
            print(f'{psfmoffat} already exists')
        else:
            ## link to the image
            inimage = os.path.join(outDir_sub, 'inimage.fits')
            if os.path.isfile(inimage):
                os.remove(inimage)
            os.symlink(inpath, inimage)
            ## run
            code_psf = os.path.join(outDir_sub, 'psf2moffat')
            proc = subprocess.run(f'(echo -1; cat {psfcat}) | {code_psf} > {psfmoffat}',
                    shell=True, cwd=outDir_sub)
            print(f'{psfmoffat} saved')

# Elapsed:26:46:22.22,User=83336.501,System=10277.366,CPU=97.1%.