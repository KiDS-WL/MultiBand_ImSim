#  measure ellipticity dependence of weights in data as a function of log(SNR)
#  and apply a correction to make the weight uniform in ellipticity in SNR bins

#  this version applies smoothing in radial and azimuthal coordinates on the e plane

#  Lance Miller November 2015

import numpy as np
import scipy.ndimage as sp
import pyfits
import optparse

np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=300)

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


def save_fits(data, path):
    """
    Save FITS data.
    :param data: FITS data to save, np.array format.
    :param path: Path to where data is to be saved.
    """
    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header
    hdulist.writeto(path, clobber=True)

if __name__ == '__main__':

    argopts = optparse.OptionParser("usage: %prog [options] arg1")
    argopts.add_option("-i", "--input", dest="inputfiles",
                       default="flist",
                       type="string",
                       help="Specify file names separted by a comma.")
    (options, args) = argopts.parse_args()

#print options.inputfiles.split(',')
    lfile = options.inputfiles.split(',')

    fpath = lfile[0] + '_weight.fits'
#    fdpath = lfile[0] + '_unsmoothed_weight.fits'
    ofile = lfile[0] + '.corr'

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

# assume only one input file 
    for i in range(1):

        with open(lfile[0]) as fileobject:
            for line in fileobject:
                if line[0]=='#' or line[0]=='!': continue
                row_len = len(line.split())
                break

        print 'line length ',row_len

        vals = np.zeros(row_len)

#  loop through data and accumulate mean variance weighted by snr**2 to deduce trend relation
        yval = np.zeros(numsnrbin)
        ny = np.zeros(numsnrbin)
        k=0
        with open(lfile[0]) as fileobject:
            for line in fileobject:
                if line[0]=='#' or line[0]=='!': continue
                vals = [float(lineobject) for lineobject in line.split()]
                snr = vals[10]
                var = vals[20]
                weight = vals[4]
                fitclass = vals[5]
                if var>0 and snr>0 and (fitclass==0 or fitclass==-9) and (weight>0 or snr<20):
                    lsnr = np.log(snr)
                    snrbin = int((numsnrbin-1)*(lsnr-minsnr)/(maxsnr-minsnr))
                    if snrbin<0: snrbin=0
                    if snrbin>=numsnrbin: snrbin=numsnrbin-1
                    yval[snrbin] += np.log(var)
                    ny[snrbin] += 1
#                    print snr,var,var*snr*snr
#  object count
                k+=1


#  smooth in 1D snr
#	print 'smoothing 1d'
        syval = np.zeros(numsnrbin)
        sny = np.zeros(numsnrbin)
        for snrbin in range(numsnrbin):
            for dm in range(yvalsmoothing):
                mbin = snrbin + dm - yvalsmoothing/2
                if mbin>=0 and mbin<numsnrbin:
                    syval[snrbin] += yval[mbin]
                    sny[snrbin] += ny[mbin]
            if sny[snrbin]>0:
                syval[snrbin] /= sny[snrbin]
#            print snrbin, np.exp(syval[snrbin]), syval[snrbin], yval[snrbin], sny[snrbin], ny[snrbin]

#  force the trend to be non-zero and not too small
        logthresh = np.log(0.01)
        for snrbin in range(numsnrbin):
            if syval[snrbin] < logthresh: syval[snrbin] = logthresh

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
        maxemod = np.zeros([numsnrbin,dtrange])

        print '! correcting ',lfile[0],' -> ',ofile
        output = open(ofile,'w')

#  loop through data again and accumulate mean weights etc
        tol = 0.01
        k=0
        with open(lfile[0]) as fileobject:
            for line in fileobject:
#  if reading header lines, write out a new header
                if line[0]=='#' or line[0]=='!': 
                    if k==5:
                        output.write( '#  5  corrected weight\n')
                    else:
                        output.write( '%s' % line ) 
                    k += 1
                    continue
#  otherwise read the catalogue values
                vals = [float(lineobject) for lineobject in line.split()]
#  bin in log(snr)
                snr = vals[10]
                e1 = vals[2]
                e2 = vals[3]
                emod = np.sqrt(e1*e1+e2*e2)
                fitclass = vals[5]
                weight = vals[4]
                var = vals[20]
#  only select galaxies with good measurements
                if (fitclass==0 or fitclass==-9) and weight>0 and var>0 and snr>0 and abs(emod-efunc_emax)>tol:
#  bin in log(SNR)
                    lsnr = np.log(snr)
                    snrbin = int((numsnrbin-1)*(lsnr-minsnr)/(maxsnr-minsnr))
                    if snrbin<0: snrbin=0
                    if snrbin>=numsnrbin: snrbin=numsnrbin-1
#  average the detrended measurement variance
                    tw = syval[snrbin]
                    var = np.log(var) - tw
#  use uncorrected ellipticity for variance correction binning
                    e1 -= vals[22]
                    e2 -= vals[23]
                    emod = np.sqrt(e1*e1+e2*e2)
                    angle = np.arctan2(e2,e1)
                    e1bin = int(emod/de)
                    e2bin = int((np.pi+angle)/(2.*np.pi*dt))
                    fdata[snrbin,e1bin,e2bin] += var
                    fw[snrbin,e1bin,e2bin] += 1.
                    if emod < emodmax:
                        totw[snrbin] += var
                        swm[snrbin] += 1
                     
# pad in emod out to the maximum value
        for snrbin in range(numsnrbin):
            for e1bin in range(drange):
                for e2bin in range(dtrange):
                    if fw[snrbin,e1bin,e2bin]>0:
                        maxemod[snrbin,e2bin] = e1bin
            for e1bin in range(drange):
                for e2bin in range(dtrange):
                    ibin = maxemod[snrbin,e2bin]
                    if e1bin>ibin:
                        fw[snrbin,e1bin,e2bin] = fw[snrbin,ibin,e2bin]
                        fdata[snrbin,e1bin,e2bin] = fdata[snrbin,ibin,e2bin]


# reflect in emod axis
        for snrbin in range(numsnrbin):
            for e1bin in range(drange):
                for e2bin in range(dtrange):
                    fdata[snrbin,2*drange-e1bin-1,e2bin] = fdata[snrbin,e1bin,e2bin]
                    fw[snrbin,2*drange-e1bin-1,e2bin] = fw[snrbin,e1bin,e2bin]

#        save_fits(fdata,fdpath)

# smooth using scipy in ellipticity
        for snrbin in range(numsnrbin):
            sp.filters.gaussian_filter(fdata[snrbin],smoothing,
                                      output=fsdata[snrbin],mode='wrap')
            sp.filters.gaussian_filter(fw[snrbin],smoothing,
				      output=fsw[snrbin],mode='wrap')
#            sp.filters.uniform_filter(fdata[snrbin],size=smoothing,
#                                      output=fsdata[snrbin],mode='constant',cval=0.)
#            sp.filters.uniform_filter(fw[snrbin],size=smoothing,
#				      output=fsw[snrbin],mode='constant',cval=0.)

#  smooth in 1D snr
#	print 'smoothing 1d'
        for snrbin in range(numsnrbin):
            for dm in range(snrsmoothing):
                mbin = snrbin + dm - snrsmoothing/2
                if mbin>=0 and mbin<numsnrbin:
                    fssdata[snrbin] += fsdata[mbin]
                    fssw[snrbin] += fsw[mbin]
                    stotw[snrbin] += totw[mbin]
                    sswm[snrbin] += swm[mbin]

# normalise
        for snrbin in range(numsnrbin):
            if sswm[snrbin]>0:
                stotw[snrbin] /= sswm[snrbin]

        thresh = 0.5/(2*np.pi*smoothing*smoothing)
        for snrbin in range(numsnrbin):
            for n in range(2*drange):
                for j in range(dtrange):
                    if fw[snrbin,n,j]>0:
                        fdata[snrbin,n,j] /= fw[snrbin,n,j]
                    if fssw[snrbin,n,j] > thresh:
                        fssdata[snrbin,n,j] /= fssw[snrbin,n,j]
                    else:
                        fssdata[snrbin,n,j] = 0.

#  output one example 
        if i==0:
            save_fits(fssdata,fpath)
#            fssw *= smoothing*smoothing
#            rpath = 'fnum.fits'
#            save_fits(fssw,rpath)

#  write new header line for extra column for original weight
        ostring = "# %2d  original weight\n" % (k)
        output.write('%s' % ostring)

#  work through data and apply correction

#  write out data
        with open(lfile[0]) as fileobject:
            for line in fileobject:
                if line[0]=='#' or line[0]=='!': 
#                    output.write( '%s' % line ) 
                    continue
                vals = [float(lineobject) for lineobject in line.split()]
#  bin in log(snr)
                snr = vals[10]
                lsnr = minsnr
                if snr > 0.: 
                    lsnr = np.log(snr)
                snrbin = int((numsnrbin-1)*(lsnr-minsnr)/(maxsnr-minsnr))
                if snrbin<0: snrbin=0
                weight = oldweight = vals[4]
                factor = 1
                if snrbin < numsnrbin:
                    e1 = vals[2]
                    e2 = vals[3]
                    emod = np.sqrt(e1*e1+e2*e2)
#  remove the autocalibration correction except when this is invalid at |e|=efunc_emax
                    if abs(emod-efunc_emax) > tol:
                        e1 -= vals[22]
                        e2 -= vals[23]
                    emod = np.sqrt(e1*e1+e2*e2)
                    angle = np.arctan2(e2,e1)
                    e1bin = int(emod/de)
                    e2bin = int((np.pi+angle)/(2.*np.pi*dt))
                    fitclass = vals[5]
                    var = vals[20]
                    if weight>0:
                        factor = np.exp( stotw[snrbin]-fssdata[snrbin,e1bin,e2bin] )
                        if emod>=efunc_emax and factor<1: factor=1
                        var *= factor
                        if var<tintervalsq: var = tintervalsq
                        weight = 2*(maxmoment-var)/(var*maxmoment + priormoment*(maxmoment-var))
                        if weight>maxweight: weight=maxweight
                        if weight<maxweight/1000: weight=0.
                        vals[4] = weight
                # fill up output string with column 5 overwritten with new weight
                output_string = "%10.6lf %10.6lf %7.4f %7.4f %7.4f %2d %7.4f %6.4f %8.1f %8.1f %8.6f %6.2f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %5.2f %5.2f %5.2f %3d %8d" % (vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7], vals[8], vals[9], vals[10], vals[11], vals[12], vals[13], vals[14], vals[15], vals[16], vals[17], vals[18], vals[19], vals[20], vals[21], vals[22], vals[23], vals[24], vals[25], vals[26], vals[27], vals[28])
                # add on the extra columns of PSF ellipticity
                if row_len>29:
                    for j in range(row_len-29):
                        output_string += ' %7.4f' % vals[j+29]
                # write the old weight as an extra column on the end
                output_string += ' %7.4f\n' % oldweight
                # write out the string
                output.write( '%s' % output_string ) 
        output.close()

