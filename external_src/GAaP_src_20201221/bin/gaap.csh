#!/bin/tcsh

# calculate the Gaap photometry from an LDAC catalogue and a GPSF image.

# parameters:
#  $1 LDAC catalogue
#  $2 GPSF image
#  $3 pre-gpsf weight image
#  $4 minimum aperture size (arcsec, added in quadrature to source size)
#  $5 maximum aperture size (arcsec)
#  $6 kernel map
#  $7 output GaaP file


if ! ($# == 7) then
 echo Usage: gaap.csh \<LDACcat\> \<GPSFimg\> \<wtimg\> \<minaper\> \<maxaper\> \<kermap\> \<out\>
 exit
endif

#set code = ~/data/shapelets/kk
#set bigim = $code/bigim
#set script = $code/gpsf
#set theli = ~/SW/theli/bin

set cat = $1
set gpsfimage = $2
set wtimage = $3
set minaper = $4
set maxaper = $5
set kermap = $6
set gaapout = $7

# output files: 
#    A shapelet map of the total ACF (orig * kernel ACF) $3.totacf.map
#        and a fits image of this ACF $3.totacf.map.fits

# Check whether the input image exists
 if ! -e $gpsfimage then
  echo ":<(" GPSF image $gpsfimage missing ===
  goto done
 endif

#set catsky = {`basename $cat`.sky}
#echo "CAT:" $cat
#echo "CATSKY:" $catsky
## Check whether there is a theli catalogue, and a filtered .sky catalogue 
# if ! -e $catsky then
#  if ! -e $cat then
#   echo ":<(" No Theli Lensfit catalogue $cat ===
#   goto done
#  endif
#  ldactoasc -i $cat -t OBJECTS -s -k X_WORLD Y_WORLD MAG_AUTO MAGERR_AUTO FLUX_RADIUS A_WORLD B_WORLD Theta SeqNr  |grep -v \#  > $catsky
# else
#  echo +++ ALREADY HAVE FILTERED $catsky ===
# endif
#
## Make input catalogue for Gaap photometry: needs X,Y,A",B",PAworld,ID
# awk '{print($1,$2)}' $catsky > radec.cat
# sky2xy  $gpsfimage @radec.cat | awk '{print($5,$6)}' > xy.cat
# paste xy.cat $catsky  | awk -v aper=$minaper -v maxaper=$maxaper '{a=$8*3600;b=$9*3600;a=sqrt(aper**2+a**2);if (a>maxaper) a=maxaper;b=sqrt(aper**2+b**2);if (b>maxaper) b=maxaper;print($1,$2,a,b,$10,$11)}' > gaap.in

# Make kernel ACF map if needed
 set keracfmap =  `basename $kermap |sed -e 's/ker.map/keracf.map/' `
 if -e $keracfmap then
  echo +++ ALREADY HAVE KERNEL ACF MAP $keracfmap ===
 else
  # check there is a kernel map frame
  if ! -e $kermap then
   echo ":<(" KERNEL MAP FILE $kermap MISSING ===
   goto done
  endif
  echo === MAKING KERNEL ACF MAP $keracfmap ===
  keracfmap < $kermap | fitkermap > $keracfmap
 endif

# Now do gapphot
 ln -sf $gpsfimage inimage.fits
 ln -sf $wtimage weights.fits
 # write GPSFSIG to gpsfsig.dat file in case GPSFSIG keyword not defined (Astrowise writes it to HISTORY kw)
 #dfits $gpsfimage| grep GPSFSIG | sed 's/.*GPSFSIG =//'|awk '{print($1)}' > gpsfsig.dat
 #gethead "HISTORY GPSFSIG" $gpsfimage > gpsfsig.dat
 (echo $keracfmap ; cat gaap.in) | gapphot > $gaapout
 \mv -f keracf.fits $keracfmap.fits
 \rm -f gpsfsig.dat inimage.fits weights.fits radec.cat
 \rm -f xy.cat gaap.in bgnoise.dat noise-est.ps

done:

exit
