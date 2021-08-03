#!/bin/tcsh

# target PSF size (sigma) needs to be specified in ker.in or as argument 3
# shapelet orders need to be specified in an orders.par file

# parameter 1 is input fits file name
# parameter 2 is input star catalogue name
# gaussianized & tweaked output file (ggpsf) is {input}.ggpsf.fits

if (($# != 2) && ($# != 3)) then
 echo Usage: gpsf-pix2gauss.csh \<image\> \<PSF catalogue\> \[GPSFSIG\]
 echo   "   (if no 3rd argument given, GPSFSIG is read from file ker.in)"
 exit
endif

set im = $1
if ! -e $im then
 echo ====  INPUT IMAGE $im NOT FOUND    ==============
 exit
endif

set psf = $2
if ! -e $psf then
 echo ====  PSF CATALOGUE $psf NOT FOUND    ==============
 exit
endif

if  ($# == 3) then
 set gpsfsig = $3
else
 if (`wc -l ker.in` == "1 ker.in") then
  set gpsfsig = `awk '{print $1}' ker.in`
 else
  echo ====  No GPSFSIG argument given and no one-line ker.in file present    ============
  exit
 endif
endif

# Build the double-shapelet kernel map

echo GAUSSIANIZING IMAGE $im USING KERNEL FROM DIRECT DOUBLE-SHAPELET FIT
echo SCALE RADIUS FOR TARGET PSF IS $gpsfsig

set kersh = $im.ker.sh2
set kermap = $im.ker.map2
set outimage = $im.ggpsf.fits

ln -sf $im inimage.fits

(echo $gpsfsig; cat $psf) | psfcat2gauskerwithtweak_no_recentre > $kersh
fitkermaptwk_noplots < $kersh > $kermap

# convolve the input image with the double-shapelet kernel map

imxshmapwithtweak < $kermap
\mv -f convolved.fits $outimage

exit
