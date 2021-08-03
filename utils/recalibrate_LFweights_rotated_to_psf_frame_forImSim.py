#  measure ellipticity dependence of weights in data as a function of log(SNR)
#  and apply a correction to make the weight uniform in ellipticity in SNR bins

#  prototype

#  Lance Miller November 2015
#  CH: upgrade for input fits catalogues 24.11.15
#  CH: gridding by PSF size and PSF e1/e2
#  CH: Added a cut on size at 0.5 pixels (23.12.15)
#  AK & LM: Fixed bug in ellipticity bins (19.07.18)
#  CH: removed ellipticity bins - this is now done at the catalogue level
#  as the K1000 catalogues were too big to read into memory (20.10.19)
#  CH: accuracy of the recalibration depends on the fineness of the PSF ellipticity bins
#  we're going to rotate the galaxy ellipticities into the PSF frame
#  and then simply bin on |emod| and PSF size
#  SSL: change catalogue format for ImSim outputs (01.07.21)
#  SSL: all print to function for Python3 (01.07.21)
#  SSL: mbin needs to be int (lines 149 & 238) (01.07.21)

import pandas as pd
import numpy as np
import scipy.ndimage as sp

import sys
from scipy import stats

np.set_printoptions(precision=6)
np.set_printoptions(suppress=True) 
np.set_printoptions(linewidth=300)

#==============================

# command line input/output files
if len(sys.argv) <3: 
    print("Usage: %s InputCat OutputCat" % sys.argv[0])
    sys.exit(1)
else:
    infile = sys.argv[1]
    outfile = sys.argv[2] 

#==============================

# Setting up fixed quantities

#  set discretisation correction (see flensfit)
tinterval = 0.02
tintervalsq = tinterval*tinterval
# give prior weight quantities
priorsigma = 0.253091
# give prior weight quantities
priormoment = 2.*priorsigma*priorsigma
efunc_emax = 0.804
maxmoment = efunc_emax*efunc_emax/2.
maxweight = 2*(maxmoment-tintervalsq)/(tintervalsq*maxmoment + priormoment*(maxmoment-tintervalsq))

# define (high resolution) log(snr) bins (these will be smoothed)
numsnrbin = 200
minsnr = np.log(5.)
maxsnr = np.log(100.)
binsize = (maxsnr-minsnr)/numsnrbin

# set sampling of ellipticity plane in emod
de = .01
drange = 1 + int(1./de)

# set azimuthal sampling of ellipticity plane
dt = .01
dtrange = 1 + int(1./dt)

#  set smoothing lengths in ellipticity and log(snr)
smoothing = 10    # ellipticity plane smoothing sigma
snrsmoothing = 30
yvalsmoothing = 10
emodmax = 0.4

#=============================

# Extract the columns of use from the catalogue
# Open feather file
feather_table = pd.read_feather(infile)

# Read in the data
snr=feather_table['SNR_LF_r'].values
var=feather_table['LS_variance_LF_r'].values
weight=feather_table['oldweight_LF_r'].values
fitclass=feather_table['class_LF_r'].values
e1_obs=feather_table['e1_LF_r'].values
e2_obs=feather_table['e2_LF_r'].values

e1_uc_obs=e1_obs-feather_table['e1_corr_LF_r'].values
e2_uc_obs=e2_obs-feather_table['e2_corr_LF_r'].values

scalelength=feather_table['scalelength_LF_r'].values - feather_table['scalelength_corr_LF_r'].values

# local PSF info
PSF_e1=feather_table['psf_e1_LF_r'].values
PSF_e2=feather_table['psf_e2_LF_r'].values
PSF_angle = np.arctan2(PSF_e2,PSF_e1)

# rotate the original and uncorrected ellipticities into the PSF frame
ct=np.cos(PSF_angle)
st=np.sin(PSF_angle)
e1_uc= e1_uc_obs*ct +e2_uc_obs*st
e2_uc= -e1_uc_obs*st + e2_uc_obs*ct
emod_uc = np.sqrt(e1_uc*e1_uc+e2_uc*e2_uc)
e1= e1_obs*ct +e2_obs*st
e2= -e1_obs*st + e2_obs*ct
emod = np.sqrt(e1*e1+e2*e2)

#initialise the recalibrated weight before the PSF loop starts
ngals=len(snr)
recalweight=np.copy(weight)  # I can now modify recalweight without changing the table weight
var_corr=var.copy()

# Select "good" objects and filter

good_obj=((var>0) & (snr>0) & np.logical_or(fitclass==0,fitclass==-9) & 
          np.logical_or(weight>0,snr<20))
var_filt=var[good_obj]
lsnr = np.log(snr[good_obj])
              
print(len(var_filt), "good objects selected out of", ngals)

#  loop through data and accumulate mean variance weighted by
#  snr**2 to deduce trend relation
            
yval = np.zeros(numsnrbin)
ny = np.zeros(numsnrbin)
snrbin = (numsnrbin-1)*(lsnr-minsnr)/(maxsnr-minsnr)
snrbin[snrbin<0] = 0
snrbin[snrbin>=numsnrbin] = numsnrbin-1
snrbin = snrbin.astype(int)

mean_g1 = np.average(e1[good_obj], weights=weight[good_obj])
mean_g2 = np.average(e2[good_obj], weights=weight[good_obj])
print("Mean shear", mean_g1, mean_g2)
              
ngals=len(lsnr)
for i in range(ngals):
    yval[snrbin[i]] += np.log(var_filt[i])
    ny[snrbin[i]] += 1

#========================================================================
#  smooth in 1D snr
syval = np.zeros(numsnrbin)
sny = np.zeros(numsnrbin)
for sbin in range(numsnrbin):
    for dm in range(yvalsmoothing):
        mbin = int(sbin + dm - yvalsmoothing/2)
        if mbin>=0 and mbin<numsnrbin:
            syval[sbin] += yval[mbin]
            sny[sbin] += ny[mbin]
    if sny[sbin]>0:
        syval[sbin] /= sny[sbin]
                
#  force the trend to be non-zero and not too small
logthresh = np.log(0.01)
syval[syval<logthresh] = logthresh

#  initialise arrays for next stage
totw = np.zeros(numsnrbin)
swm = np.zeros(numsnrbin)
stotw = np.zeros(numsnrbin)
sswm = np.zeros(numsnrbin)
fdata = np.zeros([numsnrbin,2*drange,dtrange])
fw = np.zeros([numsnrbin,2*drange,dtrange])
fsdata = np.zeros([numsnrbin,2*drange,dtrange])
fsw = np.zeros([numsnrbin,2*drange,dtrange])
fssdata = np.zeros([numsnrbin,2*drange,dtrange])
fssw = np.zeros([numsnrbin,2*drange,dtrange])
maxemod = np.zeros([numsnrbin,dtrange], dtype = np.int)

#  loop through data again and accumulate mean weights etc
tol = 0.01

# Only select galaxies with good measurements
# Note this is a different selection than the selection above
good_obj=((weight>0) & (var>0) & (snr>0) & (abs(emod-efunc_emax)>tol) & 
          np.logical_or(fitclass==0,fitclass==-9))

#  bin in log(SNR)
lsnr = np.log(snr[good_obj])
snrbin = (numsnrbin-1)*(lsnr-minsnr)/(maxsnr-minsnr)
snrbin[snrbin<0] = 0
snrbin[snrbin>=numsnrbin] = numsnrbin-1
snrbin = snrbin.astype(int)

#  average the detrended measurement variance
lvar = np.log(var[good_obj]) 

#  use uncorrected ellipticity for variance correction binning
angle = np.arctan2(e2_uc[good_obj],e1_uc[good_obj])
emod_uc_filt=emod_uc[good_obj]
e1bin = emod_uc_filt/de
e1bin = e1bin.astype(int)
e2bin = (np.pi+angle)/(2.*np.pi*dt)
e2bin = e2bin.astype(int)
            
ngals=len(lsnr)
if not ((len(e1bin)==ngals)&(len(e2bin)==ngals)): raise ValueError('SNR and ellip bins sizes do not match')

for i in range(ngals):
    tw = syval[snrbin[i]]
    fdata[snrbin[i],e1bin[i],e2bin[i]] += lvar[i] - tw
    fw[snrbin[i],e1bin[i],e2bin[i]] += 1.
    if emod_uc_filt[i] < emodmax:
        totw[snrbin[i]] += lvar[i] - tw
        swm[snrbin[i]] += 1

# pad in emod out to the maximum value
for sbin in range(numsnrbin):
    for ie1bin in range(drange):
        for ie2bin in range(dtrange):
            if fw[sbin,ie1bin,ie2bin]>0:
              maxemod[sbin,ie2bin] = ie1bin
    for ie1bin in range(drange):
        for ie2bin in range(dtrange):
            ibin = maxemod[sbin,ie2bin]
            if ie1bin>ibin:
              fw[sbin,ie1bin,ie2bin] = fw[sbin,ibin,ie2bin]
              fdata[sbin,ie1bin,ie2bin] = fdata[sbin,ibin,ie2bin]

# reflect in emod axis
for sbin in range(numsnrbin):
    for ie1bin in range(drange):
        for ie2bin in range(dtrange):
              fdata[sbin,2*drange-ie1bin-1,ie2bin] = fdata[sbin,ie1bin,ie2bin]
              fw[sbin,2*drange-ie1bin-1,ie2bin] = fw[sbin,ie1bin,ie2bin]

# smooth using scipy in ellipticity
for sbin in range(numsnrbin):
        sp.filters.gaussian_filter(fdata[sbin],smoothing,output=fsdata[sbin],mode='wrap')
        sp.filters.gaussian_filter(fw[sbin],smoothing,output=fsw[sbin],mode='wrap')

#  smooth in 1D snr
for sbin in range(numsnrbin):
    for dm in range(snrsmoothing):
        mbin = int(sbin + dm - snrsmoothing/2)
        if mbin>=0 and mbin<numsnrbin:
              fssdata[sbin] += fsdata[mbin]
              fssw[sbin] += fsw[mbin]
              stotw[sbin] += totw[mbin]
              sswm[sbin] += swm[mbin]

# normalise
stotw[sswm>0] /= sswm[sswm>0]

thresh = 0.5/(2*np.pi*smoothing*smoothing)
for sbin in range(numsnrbin):
    for n in range(2*drange):
        for j in range(dtrange):
            if fw[sbin,n,j]>0:
              fdata[sbin,n,j] /= fw[sbin,n,j]
            if fssw[sbin,n,j] > thresh:
              fssdata[sbin,n,j] /= fssw[sbin,n,j]
            else:
              fssdata[sbin,n,j] = 0.

#  work through all the data and apply correction
                        
snr[snr<=0] = 1e-9
lsnr = np.log(snr)
snrbin = (numsnrbin-1)*(lsnr-minsnr)/(maxsnr-minsnr)
snrbin[snrbin<0] = 0
snrbin[snrbin>=numsnrbin] = numsnrbin-1
snrbin = snrbin.astype(int)

ngals=len(snr)
for i in range(ngals):
    if (snrbin[i]<numsnrbin): 
              #  remove the autocalibration correction 
              # except when this is invalid at |e|=efunc_emax
        if abs(emod[i]-efunc_emax) > tol:
              angle = np.arctan2(e2_uc[i],e1_uc[i])
              ie1bin = int(emod_uc[i]/de)
              ie2bin = int((np.pi+angle)/(2.*np.pi*dt))
              emod[i]=emod_uc[i]
        else:
              angle = np.arctan2(e2[i],e1[i])
              ie1bin = int(emod[i]/de)
              ie2bin = int((np.pi+angle)/(2.*np.pi*dt))

        if weight[i]>0:
            factor = np.exp( stotw[snrbin[i]]-fssdata[snrbin[i],ie1bin,ie2bin] )
            if emod[i]>=efunc_emax and factor<1: factor=1
            fvar = var[i]*factor
            if fvar<tintervalsq: fvar=tintervalsq
            recalweight[i] = 2*(maxmoment-fvar)/(fvar*maxmoment + priormoment*(maxmoment-fvar))
            var_corr[i] = fvar
            if recalweight[i]>maxweight: recalweight[i]=maxweight
            elif recalweight[i]<maxweight/1000: recalweight[i]=0.

#We used this cut for KiDS-450.
#I considered moving this selection criteria to the "TIDY_UP" mode in
#create_catalogue_products_KV450.sh but in the end kept it here for posterity
recalweight[scalelength<0.5] = 0.

# write back to the feather table
feather_table.loc[:, 'recal_weight_LF_r']=recalweight
feather_table.loc[:, 'LS_variance_corr_LF_r']=var_corr

print("writing out file")
feather_table.to_feather(outfile)
print("finished writing out")
