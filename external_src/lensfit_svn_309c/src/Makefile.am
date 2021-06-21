AM_CFLAGS = -W -Wall -O3 -I$(top_srcdir)/swarp -I$(top_srcdir)/swarp/fits -I$(top_srcdir)/swarp/wcs -I$(top_srcdir)/nr
LDADD = -lpthread -lswarp -lwcs_c -lfits -lnr 
AM_LDFLAGS = -g -O3 -L../swarp -L../swarp/fits -L../swarp/wcs -L../nr 
bin_PROGRAMS = makeospsf makeglobalpsf globalshifts lensfit flensfit readlensfit createPSFcube remakeospsf

makeospsf_SOURCES = makeospsf.c \
	mediansub.c \
	xcorr_measure.c \
	shiftft.c \
	extractpostagestamp.c \
	tbswarpextract.c \
	tbswarpextract_nodistortion.c \
	weightfilter.c \
	getimagesize.c \
	readpsfcat_simple_wcs.c \
	readpsfcatsize.c \
	writeglobalcoeffs.c \
	writeoscoeffs.c \
	globalsvdfit.c \
	varylpthresholdf2.c \
	writepsf_fits.c \
	globalreconstruct.c \
	svdfit2dsqc.c \
	writepsfstamp.c \
	util.c \
	psfmoments.c \
	updateglobalcoeffs.c \
	updateoscoeffs.c \
	lensutil.h

remakeospsf_SOURCES = remakeospsf.c \
	mediansub.c \
	xcorr_measure.c \
	shiftft.c \
	extractpostagestamp.c \
	tbswarpextract.c \
	tbswarpextract_nodistortion.c \
	weightfilter.c \
	getimagesize.c \
	readpsfcat_simple_wcs.c \
	readpsfcatsize.c \
	writeglobalcoeffs.c \
	writeoscoeffs.c \
	globalsvdfit.c \
	varylpthresholdf2.c \
	writepsf_fits.c \
	globalreconstruct.c \
	svdfit2dsqc.c \
	writepsfstamp.c \
	util.c \
	psfmoments.c \
	updateglobalcoeffs.c \
	getospsfsizes.c \
	updateoscoeffs.c \
	readoscoeffs.c \
	lensutil.h

createPSFcube_SOURCES = createPSFcube.c \
	tbswarpextract.c \
	svdfit2dsqc.c \
	getimagesize.c \
	getglobalpsfsizes.c \
	getospsfsizes.c \
	getglobalcoeffs.c \
	readoscoeffs.c \
	readpsfcat_simple_wcs.c \
	readpsfcatsize.c \
	globalreconstruct.c \
	psfmoments.c \
	util.c \
	lensutil.h

makeglobalpsf_SOURCES = makeglobalpsf.c \
	mediansub.c \
	xcorr_measure.c \
	shiftft.c \
	extractpostagestamp.c \
	tbswarpextract.c \
	tbswarpextract_nodistortion.c \
	weightfilter.c \
	getimagesize.c \
	readpsfcat_simple_wcs.c \
	readpsfcatsize.c \
	writeglobalcoeffs.c \
	globalsvdfit.c \
	varylpthresholdf2.c \
	writepsf_fits.c \
	globalreconstruct.c \
	svdfit2dsqc.c \
	util.c \
	lensutil.h

globalshifts_SOURCES = globalshifts.c \
	getglobalpsfsize.c \
	getcoeffs_fits_expand.c \
	updateglobalcoeffs.c \
	psfmoments.c \
	tbswarpextract_nodistortion.c \
	swarpextract.c \
	weightfilter.c \
	mediansub.c \
	xcorr_measure.c \
	shiftft.c \
	extractpostagestamp.c \
	getimagesize.c \
	readpsfcat_simple_wcs.c \
	readpsfcatsize.c \
	writeglobalcoeffs.c \
	varylpthresholdf2.c \
	writepsf_fits.c \
	svdfit2dsqc.c \
	globalsvdfit.c \
	globalreconstruct.c \
	util.c


lensfit_SOURCES = tlensfit.c \
	mediansub.c \
	extractpostagestamp.c \
	tbswarpextract_noflag.c \
	tbswarpextract_distortion_measure.c \
	getimagesize.c \
	readgalcatsize.c\
	readpsfcat.c \
	readpsfcatsize.c\
	getcoeffs_fits_expand.c \
	getglobalcoeffs.c \
	getpsfsize.c \
	xcorr_measure.c \
	shiftft.c \
	svdfit2dsqc.c \
	reconstruct.c \
	psfmoments.c \
	psfellipticitycentroid.c \
	convolvemodel.c \
	getglobalpsfsizes.c \
	globalreconstruct.c \
	writepsfstamp.c \
	util.c

flensfit_SOURCES = flensfit.c \
	mediansub.c \
	extractpostagestamp.c \
	tbswarpextract_noflag.c \
	tbswarpextract_distortion_measure.c \
	getimagesize.c \
	readgalcatsize.c\
	readpsfcat.c \
	readpsfcatsize.c\
	getcoeffs_fits_expand.c \
	getglobalcoeffs.c \
	getospsfsizes.c \
	readoscoeffs.c \
	getpsfsize.c \
	xcorr_measure.c \
	shiftft.c \
	svdfit2dsqc.c \
	reconstruct.c \
	psfmoments.c \
	psfellipticitycentroid.c \
	convolvemodel.c \
	getglobalpsfsizes.c \
	globalreconstruct.c \
	writepsfstamp.c \
	util.c

readlensfit_SOURCES = readlensfits.c \
	reconstruct.c \
	svdfit2dsqc.c
