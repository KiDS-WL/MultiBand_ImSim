#!/bin/tcsh

# parameter 1 is input fits file name
# parameter 2 is input star catalogue name
# gaussianized output file (gpsf) is {input}.gpsf.fits
# gaussianized & tweaked output file (ggpsf) is {input}.ggpsf.fits

echo ====================
echo PSF-GAUSSIANIZING IMAGE $1

set im  = $1
if ! -e $im then
 echo ====  INPUT IMAGE $im NOT FOUND    ==============
 exit
endif

set psf = $2
if ! -e $psf then
 echo ====  PSF CATALOGUE $psf NOT FOUND    ==============
 exit
endif

set psfsh = $im.psf.sh
set psfmap = $im.psf.map
set kersh = $im.ker.sh
set kermap = $im.ker.map

#echo 10 > orders.par; echo 4 >> orders.par

if (! -e $psfmap || -z $psfmap) then
 echo BUILDING SHAPELET-BASED PSF MAP
 \rm -f inimage.fits; ln -s $im inimage.fits
 sort -n -k 2 $psf | psfcat2sherr_2 > $psfsh
 fitpsfmap2 < $psfsh > $psfmap
else
 echo ====  USING EXISTING PSF MAP $psfmap
endif

if (! -e $kermap || -z $kermap) then
 echo BUILDING SHAPELET-BASED GAUSSIANIZATION KERNEL
 if (! -e ker.in) then
     echo -0.8 > ker.in
 endif
 \rm -f inimage.fits; ln -s $im inimage.fits
 gaussianizepsf < $psfmap > $kersh
 fitkermap < $kersh > $kermap
 \cp -f ker.in gpsfsig.dat
else
 echo ====  USING EXISTING GAUSSIANIZATION KERNEL $kermap
 head -1 $kermap |awk '{print($2)}' > gpsfsig.dat
endif

set gim  = $im.gpsf.fits
if ! -e $gim then
 \rm -f inimage.fits; ln -s $im inimage.fits
 echo BUILDING GAUSSIANIZED IMAGE
 imxshmap < $kermap
 \mv -f convolved.fits $gim
else
 echo ====  ALREADY HAVE GAUSSIANIZED IMAGE $gim
endif

set beta2 = `awk '{print($1*2.0)}' gpsfsig.dat`
echo SCALE RADIUS FOR PSF TWEAK IS $beta2

set gpsfsh = $im.gpsf.sh
set gpsfmap = $im.gpsf.map

if (! -e $gpsfmap || -z $gpsfmap) then
 echo COMPUTING SHAPELET-BASED PSF MAP FOR GAUSSIANIZED IMAGE
 \rm -f inimage.fits; ln -s $gim inimage.fits
 (echo $beta2; sort -n -k 2 $psf) | psfcat2sherr_fixbeta > $gpsfsh
 #(echo 0 $beta2; sort -n -k 2 $psf) | psfcat2sherr_fixbgbeta > $gpsfsh
 fitpsfmap2 < $gpsfsh > $gpsfmap
else
 echo ====  USING EXISTING PSF MAP FROM GAUSSIANIZED IMAGE $gpsfmap
endif

set ggim  = $im.ggpsf.fits
if ! -e $ggim then
 echo TWEAKING GAUSSIANIZATION
 \rm -f inimage.fits; ln -s $gim inimage.fits
 (cat gpsfsig.dat $gpsfmap) | psfmapmingaus > $im.delker.map
 imxshmaptweak < $im.delker.map
 \mv -f convolved2.fits $ggim
else
 echo ====  ALREADY HAVE TWEAKED GAUSSIANIZED IMAGE $ggim
endif

exit
