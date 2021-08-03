#!/usr/local/bin/tcsh

setenv GAAP /disks/shear15/kuijken/
setenv Gbin $GAAP/bin

if ! ($# == 5) then
 echo Usage: gpsf-gaap-bin.csh \<PSFimage\> \<GALimage\> \<GALcat\> \<MinAper\> \<MaxAper\>
 exit
endif

set psfimage = $1    # fits image containing bright stars only
set image = $2       # fits image containing galaxies to be photometered
set galaxycat = $3   # SExtractor catalogue from detection image
set minaper = $4     # Aperture to add in quadrature to GAaP aperture (arcsec)
set maxaper = $5     # Maximum aperture radius to  (arcsec)

# names for related files: input weight map (if available - otherwise estimated)
set wtimage     = `basename $image .fits`.weight.fits
# Intermediate and output files
set psfcat      = `basename $psfimage .fits`.psf.cat
set gpsfimage   = `basename $image .fits`.ggpsf.fits
set kermap      = `basename $image .fits`.ker.map
set keracfmap   = `basename $image .fits`.keracf.map
set gaapout     = `basename $image .fits`.gaap.refltheta

# 1. make star catalogue in right SExtractor ascii format
echo
echo ------------- 1. STAR CATALOGUE FROM $psfimage -----------
if -z $psfcat rm -f $psfcat
if ! -e $psfcat then
 if  ! -e $psfimage then
  echo ":<(" PSF image $psfimage missing!! ===
  goto done
 endif
 if ! -e gaap.sex cp $GAAP/lib/gaap.sex .
 if ! -e gaap.sexpar cp $GAAP/lib/gaap.sexpar .
 if ! -e default.conv cp $GAAP/lib/default.conv .
 sex  $psfimage -c gaap.sex -PARAMETERS_NAME gaap.sexpar -DETECT_THRESH 10 -CATALOG_NAME $psfcat
 echo === MADE PSF CATALOGUE $psfcat ===
else
 echo +++ ALREADY HAVE PSF CATALOGUE $psfcat ===
 endif

set gpsfsig = -1

# 2. calculate Gaussianization kernel and Gaussianised image
# allow Gaussianised PSF radius to be determined autimatically
# kernel has shapelet order 8, no spatial variation for the PSF
echo
echo ------------- 2. PSF-GAUSSIANISING $image -----------
if -e $gpsfimage then
 echo +++ ALREADY HAVE PSF-GAUSSIANISED GALAXY IMAGE $gpsfimage ===
else
 if -z $kermap rm -f $kermap
 if ! -e $kermap then
  echo 8 > orders.par
  echo 0 >> orders.par
  ln -sf $psfimage inimage.fits
  (echo $gpsfsig; cat $psfcat) |  $Gbin/psfcat2gauskerwithtweak_no_recentre > tmp_ker.sh2
  $Gbin/fitkermaptwk_noplots < tmp_ker.sh2 > $kermap
  if -z $kermap then
   echo ":<(" kernel map $kermap has zero size!! ===
   goto done
  endif
  echo === MADE GAUSSIANISATION KERNEL $kermap ===
 else
  echo +++ ALREADY HAVE GAUSSIANISATION KERNEL $kermap ===
 endif
 if ! -e $image then
  echo ":<(" galaxy image $image is missing!! ===
  goto done
 endif
 if -e convolved.fits rm -f convolved.fits
 ln -sf $image inimage.fits;  $Gbin/imxshmapwithtweak < $kermap
 mv -f convolved.fits $gpsfimage
 echo === MADE PSF-GAUSSIANISED GALAXY IMAGE $gpsfimage ===
endif


# 3. calculate GAaP photometry
# Make input catalogue for Gaap photometry from the galaxy catalogue:
# want text file with X,Y,A",B",PAworld,ID. where PAworld = -THETA_WORLD
echo
echo ------------- 3. GAaP PHOTOMETRY
echo ------------- \ \ \ \ \ IMAGE: $gpsfimage
echo ------------- \ \ \ \ \ CAT:   $galaxycat
if -z $gaapout rm -f $gaapout
if -e $gaapout then
 echo +++ ALREADY HAVE GAaP CATALOGUE $gaapout ===
else
 echo ------- PREPARE GAaP INPUT FILE ---------------
 set catsky = {`basename $galaxycat`.sky}

# make catalogue of X_WORLD Y_WORLD MAG_AUTO MAGERR_AUTO FLUX_RADIUS A_WORLD B_WORLD Theta SeqNr in $catsky
# note that the position angle is -1 * THETA_WORLD
 awk '{print $3,$4,$9,$10,$19,$23,$24,$26,$33}' $galaxycat |grep -v \# > $catsky

# Make input catalogue for Gaap photometry: needs X,Y,A",B",PAworld,ID
 awk '{print($1,$2)}' $catsky > radec.cat
 sky2xy  $gpsfimage @radec.cat | awk '{print($5,$6)}' > xy.cat
 paste xy.cat $catsky  | awk -v aper=$minaper -v maxaper=$maxaper '{a=$8;b=$9;a=sqrt(aper**2+a**2);if (a>maxaper) a=maxaper;b=sqrt(aper**2+b**2);if (b>maxaper) b=maxaper;print($1,$2,a*3600,b*3600,$10,$11)}' > gaap.in
 rm -f xy.cat radec.cat
# Make kernel ACF map of GPSF image if needed
 echo ------- PREPARE NOISE AUTOCORRELATION FUNCTION ----------
 if -z $keracfmap rm -f $keracfmap
 if -e $keracfmap then
  echo +++ ALREADY HAVE KERNEL ACF MAP $keracfmap ===
 else
  # check there is a kernel map frame
  if ! -e $kermap then
   echo ":<(" KERNEL MAP FILE $kermap MISSING ===
   goto done
  endif
  echo === MAKING KERNEL ACF MAP $keracfmap ===
  if ! -e orders.par then
    head -2 $kermap | awk '{print($1)}' > orders.par
  endif
  $Gbin/keracfmap < $kermap | $Gbin/fitkermap > $keracfmap
 endif

# Estimate pre-GPSF pixel correlation and convolve with kernel ACF
 set totacfmap =  `basename $kermap |sed -e 's/ker.map/totacf.map/' `
 if ! -e $totacfmap then
  ln -sf $image inimage.fits
  ln -sf $wtimage weights.fits
  echo NOMAP|$Gbin/gapphot > /dev/null
  set acfsigorig = `head -1 keracffitted.map|awk '{print($2)}'`
  set acfsigconv = `head -1 $keracfmap |awk '{print($2)}'`
  \mv -f cov.fits $image.cov.fits
  \mv -f covsh.fits $image.covsh.fits
  echo === APPROXIMATED NOISE ACF WITH A GAUSSIAN OF WIDTH $acfsigorig
  (echo $acfsigorig $acfsigconv; cat $keracfmap) |$Gbin/kermapxgauss > $totacfmap
  \rm -f keracffitted.map keracffitted_pknorm.map
  echo === MADE TOTAL NOISE ACF MAP $totacfmap ===
 else
  echo +++ ALREADY HAVE TOTAL NOISE ACF MAP $totacfmap ===
 endif
 echo
 echo ------------- GAaP PHOTOMETRY RUN --------------------
# Now do gapphot
 ln -sf $gpsfimage inimage.fits
 ln -sf $wtimage weights.fits
 # write GPSFSIG to gpsfsig.dat file in case GPSFSIG keyword not defined (Astrowise writes it to HISTORY kw)
 $Gbin/dfits $gpsfimage| grep GPSFSIG | sed 's/.*GPSFSIG =//'|awk '{print($1)}' > gpsfsig.dat
 (echo $totacfmap ; cat gaap.in) | $Gbin/gapphot > $gaapout
 \mv -f keracf.fits $totacfmap.fits
 \rm -f gpsfsig.dat inimage.fits weights.fits xy.cat gaap.in radec.cat
 \rm -f bgnoise.dat noise-est.ps
 echo === GAaP CATALOGUE $gaapout READY ===
endif

done:
exit
