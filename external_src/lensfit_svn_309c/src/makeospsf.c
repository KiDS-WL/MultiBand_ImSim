/*-----------------------------------------------------------------

makeospsf 

Code to create a surface fit to the PSF, reusing code from previous
stacking algorithms, makepsf and makeglobalpsf.  Here's the preamble 
and history for makepsf and makeglobalpsf:

-------------------------------------------------------------------

[Code to create PSFs suitable for lensfit, by stacking stars in an
image.  Requires a list of input star positions plus the matching
image.  Stars with max brightness more than 50 percent of the CCD
saturation level are not used.  Magnitude ranges and peak S/N may also be
specified.  A list of galaxy positions is also used to eliminate
neighbouring galaxies.  Stars are sub-pixel registered onto the
stacked PSF in a two-step interative process.  Individual stars are
cross-correlated with the stacked PSF and any that are a poor match
are thrown out.  The PSF is assumed to be position dependent, and is
defined in a series of boxes dividing up the input image.  Stars are
accumulated into the PSF for each box if they fall within some
distance of that box.  The sampling box size and smoothing are both
specified on the command line.  Output is a three-dimensional FITS
image, the first two dimensions are x and y position on the image,
the next dimension is used for each 2D PSF stacked image.
There is also a table of the x and y positions of the stars that have
been used for each PSF.]

The last part is replaced by a polynomial fit, which in this code is
a fit across all chips in an exposure, with optional chip-dependent
variation

History:

original makepsfsurface first created Lance Miller. 20/01/2007
modified to crosscorrelate initially with a delta function to optimise centering
LM 27/02/2007
modified to write out single FITS file containing all PSFs, plus smoothing
algorithm and additional options LM, 17/1/2008
command-line input modified to force S/N ratio to be compulsory, also
change to output information and time estimate included LM, 30/1/2008
addition of background value into extractdata requires change to makepsf
also, needed to get correct saturation value for bright stars LM, 30/1/2008
modification to getimagesize to return CCD PA (not used in makepsf
but needed for function compatibility) LM, 4/3/2008
Change to fit polynomial surface to pixel values LM, 7/5/2008
Also read bad pixel mask LM 12/5/2008
Improve cross-correlation method  LM 14/7/2008
Made compatible with new extractdata routine LM 23/7/08
Moved to extractpostagestamp routine LM 23/9/08
Added in option of WCS coords LM 26/9/08
Allowed option of selecting additional stars to be used LM sep 08
Modified memory allocation for varylpthresholdf arrays LM Dec 08
Modified to allow swarp correction of astrometric distortion, for
consistency with lensfit v5, LM 16/01/2009
Now writes info on numbers of stars used in each psf box into the 
coefficients file, so that"bad" regions may be flagged out in use.
LM 21/1/2009
This version make global polynomial fit to the ensemble of chips
in an exposure, with optional chip-dependence of low-order 
coefficients
LM Jan 2009
Corrected inconsistencies between swarp and non-swarp versions
LM 4/3/09
Various minor changes for consistency with lensfit modifications, LM  Apr-Oct 09
Corrected bug in output RA,dec values (stars were being muddled-up) LM 4 Nov 2009
ported to SVN system LM 5th August 2010
aligned with May 7 verion, LM 7th August 2010
corrected the delimiter bug and aligned with v7, LM 17 Nov 2010

----------------------------------------------------------------------------

makeospsf

New version fits oversampled PSF model.  Model fitting proceeds by shifting model
onto stars by (linear) interpolation, rather than the other way round.  Polynomial
pixel model coefficients are jointly fit.  Image gain is used to get better estimate
of pixel noise and thus enable a chi-squared fit criterion to be employed.

Revision $Rev: 253 $ 
last changed by $Author: miller $
on $LastChangedDate: 2014-10-28 13:42:26 +0000 (Tue, 28 Oct 2014) $

---------------------------------------------------------------------------*/

#define DATE_STRING "$LastChangedDate: 2014-10-28 13:42:26 +0000 (Tue, 28 Oct 2014) $"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

/* swarp wcs include files (but note custom version of fitscat.h to avoid conflict with cfitsio */
#include "define.h"
#include "types.h"
#include "globals.h"
#include "./fitscat.h"
#include "fitswcs.h"
#include "data.h"
#include "field.h"
#include "header.h"
#include "interpolate.h"
#include "prefs.h"
#include "projapprox.h"
#include "resample.h"
#include "weight.h"
#include "wcs/wcs.h"
#include "prefs.h"

#include "lensutil.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// declare GSL random-number generator
const gsl_rng_type *T;
gsl_rng *r;

// machine memory limit (MB) (not tested at present)
// #define MEMORY_LIMIT 14000
// limit above which individual memory allocations will be reported (MB)
#define MEMORY_REPORT_LIMIT 1000

int extractpostagestamp(float *, float *, int *, double *, double *,
                        float *, float *, float, float, float,
                        float, float, int, int, float, int *, float);
void svdfit2dsqc(double *, double *, double *, double *, int, int, int,
                 double *, double **, double **, double *);
int tbswarpextract(wcsstruct *, float *, float *, int *, double *, double *,
                   int *, float, float, float, int, int, ikernelstruct *,
                   float, int *, float, double *, double *, double *,
                   double *, double *, double **, double *, double *,
                   double *, double *, double **, double **);
int tbswarpextract_nodistortion(wcsstruct *, float *, float *, float *, int *,
                                double *, double *, float, float,
                                float, int, int, float,
                                int *, double *, float, double *, double *,
                                double *, double *, int);
int writeglobalcoeffs(char *, int, double **, int, int, int, int, double *, int,
                      int, int, int, int *, int, int, int, double, double, float);
int writeoscoeffs(char *, int, int*, double*, int, int, int,
		  int, int, double*, int, int, int, int, int,
		  int*, int, int, int, double, double, float);
int writepsf_fits(char *, double **, int, int, int, int, int *, float **,
                  float **, int, int);
void psfmoments(double *, int, int, double *, double *, double *);
int readpsfcatsize(char *);     // read size of catalogue
/* choose appropriate routine for catalogue format */
//int readpsfcat (char*, float*, float*, float*);
//int readpsfcat2_wcs (char*, double*, double*, float*);
int readpsfcat_simple_wcs(char *, double *, double *, float *);
int getimagesize(char *, char *, int *, float *, float *, float *, float *,
                 int *, int *, float);
void getdata_badpix(char *, char *, int *, float *, float *, float *);
int weightfilter(float *, float *, float *, float *, int *, int);
void dilatemask(float *, int *);
void globalsvdfit(double *, double *, double *, double *, int *, int, int,
                  int, int, int, int, int, double *, double **, double **,
                  double *);
void globalreconstruct(double, double, int, int, int, int, int, int, double *,
                       double *);
int updateglobalcoeffs(char *, int, int, int, int, double *, double *);
int updateoscoeffs(char *, int, int, int, int, double, double *, double *);
int varylpthresholdf2(float *, int *, float, int, int, int, int, int, float,
                      float, int);
void mediansub(float *, float *, int *, float, float *);
int wcs_to_raw(wcsstruct *, double *, double *);
void load_dfield(fieldstruct *, catstruct *, char *, int, int);
ikernelstruct *init_ikernel(interpenum *, int);
char *create_psf_name(const char *, char *, const char *);
char *create_ospsf_name(const char *, char *, const char *);
void get_poly_order(int, int, int, int, int *, int *);
void crosscorrelate(int, double*, double*, int, double*, double*, double*);
double chisquared(double*, int, double*, double*, double*);
void downsample();
int initialisePSF();
int reinitialisePSF();
void enablechipvariation(int);
void randomisesoln();
int smatrix(double*);
int writepsfstamp(int, int, double*, char*);
void reconstructosPSF(double, double, int);
double solvePSF();
void my_fdf (const gsl_vector*, void*, double*, gsl_vector*);
double my_f (const gsl_vector*, void*);
void my_df (const gsl_vector*, void*, gsl_vector*);
int open3Dfits(char*, fitsfile**, int);
int write3Dfits(fitsfile*, int, double*);
double nonlinearity();
double storepixelvalues();
double reset_weights();
void downsamplePSFmodel(char*, int, int, int, 
			int, int, int*, double*,
			double, double);

static int compare(const void *ii, const void *jj)
{
    float *i, *j;
    i = (float *) ii;
    j = (float *) jj;
    if (*j > *i)
        return -1;
    if (*j < *i)
        return 1;
    return 0;
}

/* verbose flag */
#define VERBOSE 0

/* set whether or not to subtract a constant median background 
   (0=no 1=yes) */
#define SUBTRACT_MEDIAN 0

/* set whether or not to apply the median erosion filter to the weight image
   (same routine as weightimagefilter, used on CFHTLS weight images) */
#define FILTER_WEIGHTS 1

/* specify that weight (bad pixel) images are required */
#define WEIGHTS_REQUIRED 1

/* define whether input galaxy catalogue positions are in WCS (WCS=1)
   or pixel (WCS=0) coordinates */
#define WCS 1

/* define whether or not to use swarp astrometric corrections
   requires WCS 1 */
#define USE_SWARP 1

/* define whether to fit systematic astrometric error (1=normal use) or set these to zero 
   (0=assume that astrometric errors are zero) */
#define FIT_ASTROMETRIC_ERRORS 1

/* define that the code is being used on simulations such as GREAT3 in which a number of operations
   needed for real data should be switched off (1=assume GREAT3, otherwise assume normal data) */
#define GREAT3 0

/* set whether or not to apply the internal swarp distortion correction (1=yes) 
   requires WCS=1 and USE_SWARP=1 also */
#define CORRECT_DISTORTION 0

/* define whether input files are gzipped fits or plain fits
   0=not zipped, 1=weight files are zipped, 2=both data and weight files zipped */
#define GZIPPED 0

/* define minimum accepted exposure time (must be the same as in lensfit) */
#define minexposuretime 0.

/* Set the instrument here by defining the flag that is checked later
   in the code KIDS, CFHT, SUPRIME or SINGLE */
#define KIDS (1)
/*---------------------------------------------------------------------------*/

/* global variables */
float poserror;
int correct_distortion;
int pwidth, pheight, fwidth, fheight; // sizes of native and oversampled arrays
int oversampling=3; // final oversampling factor (must be odd integer, not too large)
int ofactor; // temporary oversampling factor
int current_star=0; // identifies which star is being sent to function chisquared
double rmax=6.; // maximum radius of region to be oversampled
double sigma=1.5; // sigma of initial trial gaussian PSF (used to get initial star centroids)
double wsigma=2.0; // sigma of weight function gaussian for measuring astrometric shifts
// wsigma should be fairly broad because we should match it to what was used 
// to measure star centroids in the input Sextractor catalogue
int dsize;
int nfit, nfitmax, nchip, *rchip, *nstar, order, ncoeffs, *pindex, *oldpindex, *npix, **spix, *zpix; 
int chipvariation, chiporder;
double **xval, **yval, **sval, **zval, *psf, *dpsf, *respsf, *fracrespsf, *refpsf, ***dpix, ***ppix;
double *xchipval, *ychipval;
double ***shift, ***bestshift, ***initshift;
double **chipxval, **chipyval;
double *soln, *oldsoln, *basesoln;
float *gain, *noise;
int ncall_f, ncall_df, ncall_fdf;
double **u, **v, *w, *avals;
double *xfit, *yfit, *zfit, *wfit, *zsort;
double bval;
double **starweight, **chisqval;
double nullchisq, bestchisq;
double deltachisqlimit;
double systematic=0.0; // NO fractional error to be added as systematic noise term
int fsize;  // size of oversampled PSF image
int writeout=1;
int printout;
int os_iter, chipvariation_loop;

void outputdstamp(char *fname)
{
  int i, m, x, y, nset;
  double pos[2], *testpsf;

  testpsf = (double*)calloc(pwidth*pheight, sizeof(double));

  reconstructosPSF(0.,0.,nchip/2); // Outputs reconstructed (oversampled) PSF array psf[] at centre of field (given global focal-plane position and chip number)
  pos[0]=pos[1]=0.;
  nset = smatrix(pos); // Specify that the PSF should be centred by giving arguments of 0 (representing the amount the PSF should be shifted by on the postage stamp)
  downsample(); // Take oversampled PSF psf[] and create downsampled PSF dpsf[]
    for (i=0; i<pwidth*pheight; i++)
      {
	m = i;
	y = m/pwidth;
	x = m - y*pwidth;
	x += pwidth/2;
	y += pwidth/2;
	if (x>=pwidth) x -= pwidth;
	if (y>=pheight) y -= pheight;
	m = y*pwidth + x;
	testpsf[m] = dpsf[i];
      }
    remove(fname);
    writepsfstamp(pwidth,pheight,testpsf,fname);
    free (testpsf);
}

void outputfstamp(char *fname)
{
  int i, m, x, y;
  double *testpsf;

  testpsf = (double*)calloc(fwidth*fheight, sizeof(double));

  reconstructosPSF(0.,0.,nchip/2);
    for (i=0; i<fwidth*fheight; i++)
      {
	m = i;
	y = m/fwidth;
	x = m - y*fwidth;
	x += fwidth/2;
	y += fwidth/2;
	if (x>=fwidth) x -= fwidth;
	if (y>=fheight) y -= fheight;
	m = y*fwidth + x;
	if (pindex[i]>=0)
	  testpsf[m] = psf[pindex[i]];
	else
	  testpsf[m] = 0.;
      }
    remove(fname);
    writepsfstamp(fwidth,fheight,testpsf,fname);
    free (testpsf);
}


void findmodelcentroid(double *centroid, double msigma)
{
  // function to recontruct a full-resolution psf model and iteratively
  // measure the weighted centroid (isotropic gaussian weight function the same as the 
  // initial psf model)

  int i, x, y, iter, niter;
  double xx, yy, xcen, ycen, sumw, pval, xxcen, yycen;
  double rsqos, sigmasq, w, ofactorsq, dist, distsq;

  // squared oversampling factor
  ofactorsq = ofactor*ofactor;

  // sigma squared of Gaussian weight function
  sigmasq = msigma*msigma;

  // distance limit for moments summations
  dist = (double)fwidth/2;
  distsq = dist*dist;

  // (iteratively) find the weighted model centroid
  niter = 1;
  centroid[0] = centroid[1] = 0.;
  for (iter=0; iter<niter; iter++)
    {
      xcen = ycen = sumw = 0.; // initialise summations for this iteration
      for (y=0; y<fheight; y++)
	{
	  yy = y>fheight/2 ? y-fheight : y; // allow for swapped quadrants
	  yycen = yy - centroid[1]; // centre the weight function on the current best guess centroid
	  for (x=0; x<fwidth; x++)
	    {
	      i = y*fwidth + x;
	      if (pindex[i]>=0)
		{
		  xx = x>fwidth/2 ? x-fwidth : x; // allow for swapped quadrants
		  dist = xx*xx + yy*yy;
		  if (dist < distsq) // force outer boundary to be circular centred on zero
		    {
		      xxcen = xx - centroid[0]; // centre the weight function on the current best guess centroid
		      rsqos = xxcen*xxcen + yycen*yycen; // radius squared defined in oversampled units relative to best-guess centroid
		      w = exp(-rsqos/2./sigmasq/ofactorsq); // weight function
		      pval = psf[pindex[i]];
		      xcen += xx*pval*w;
		      ycen += yy*pval*w;
		      sumw += pval*w;
		    }
		}
	    }
	}
      if (sumw>0.)
	{
	  centroid[0] = xcen/sumw;
	  centroid[1] = ycen/sumw;
	}
      else
	{
	  fprintf(stderr," failed to find centroid, all values zero\n");
	  exit(EXIT_FAILURE);
	}
    }

}


unsigned long int random_seed()
{

 unsigned int seed;
 struct timeval tv;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
   //printf("Got seed %u from gettimeofday()\n",seed);
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
   //printf("Got seed %u from /dev/random\n",seed);
   fclose(devrandom);
 }

 return(seed);

}


static int
dcompare(const void *ii, const void *jj)
{
    double *i, *j;
    i = (double *) ii;
    j = (double *) jj;
    if (*j > *i)
        return -1;
    if (*j < *i)
        return 1;
    return 0;
}



/* main program */
int main(int argc, char *argv[])
{
    int halfpwidth, halfpheight, dim[2], badccd, imageid;
    int padwidth, padheight, halfpadwidth, cchwidth;
    int chipnumber;
    int flagbad, nset, range;
    int i, j, ii; //, p, q;
    int x, y;
    int xmin, xmax, ymin, ymax;
    int imagesize; // size of input CCD
    int len, nobj, nobjt, nobj2, nrefpsf, ncount;
    int psfsize;
    int *numstar;
    int tot_nstar;
    int ix, iy, pixel, xx, yy;
    int num_off, num_wrong_mag, num_saturated, num_close, num_snratio;
    int num_badpix;
    int close;
    int d;
    int ngoodpix, imax_intensity;
    int nn;
    int *region, *objpix;
    int nobjects;
    int lower_area_limit, lower_merging_limit, arealimit;
    int *offscale;
    int naxis = 2;
    int narg = 0;
    int fmheight = 20;  // half-height of median filter in weightfilter
    int nchanged;
    int k, kimage, refimage, refchip;
    int image, xchipsampling, ychipsampling, big_gap, nchip_expected,
        nxchip, nychip;
    int xchip, ychip, *chips;
    int l;
    int worstpix;
    int os_niter, oversampling_loop, bval_loop;
    int newchiporder, targetchiporder;
    int copied, randomise_loop;
    int sorder, schipvariation, schiporder, scrossterm, sncoeffs;

    unsigned long int seedval;

    time_t t1, t2, t3, t4, t5;
    float memory, memloop;

    double norm_star, norm_model;
    double *xfits, *yfits, *wfits;
    double *xshiftfits, *yshiftfits, *bxshiftfits, *byshiftfits, *ow, *rx, *ry;
    double xsize, ysize, hxsize, hysize;
    double *radegs, *decdegs;   // ra,dec in input catalogue
    double *rdegs, *ddegs;    // ra,dec arrays for selected stars on each chip
    double **chiprdegs, **chipddegs;    // ra,dec arrays for selected stars on each chip
    double scalefactor[1];
    double cosmicraylimit;
    double ressum, worstpixval, worstpixratio, nval;
    double cent[2], tempshift[2];
    double offset, chisq_expected;
    double xrms, yrms;
    double pos[2], psfe[2], moments[3], xs, ys, rdx, rdy;
    double sum, sumx, sumy, sumxsq, sumysq;
    double *xcoeffs, *ycoeffs;
    double *currentsoln;

    float maxlevel;
    double noisevalue;
    float *opix, *changed, *weightfilterarray;
    float *apix, *temp, *badpix, *badtemp;
    float *sortarray, **sn;
    // float snref=50.;  // reference S/N ratio for star weighting function
    float fintensity_limit, fmax_intensity, scale;
    float *objx, *objy, *mag, *rmag;
    float brightmag, faintmag, snratio;
    float satlev, satlimit, arcperpix, angle;
    float ximageoffset, yimageoffset;
    double *dobjx, *dobjy, *refpsfx, *refpsfy, *refpsfx_glob, *refpsfy_glob;
    int crossterm, maxncoeffs, refpsfnset;
    int ref_px, ref_x, ref_y, ref_i;

    double **acoeffs;
    double ***dbadpix;
    double maxvalue, distantmaxvalue;
    double psfnorm;
    double lchisq, lchisq0, amp, modelsum, cross_sum;
    double rawpos[2], wcspos[2], wcsneg[2], wcscentre[2], *kern[2], refpsfpos[2];

    char readstring[2000];
    char readdelims[] = " \t";
    char *evariable, *headevariable, *badpixdir, *psfdir, *readitem = NULL;
    char *prefsname;
    char *catname;
    char **image_file;
    char *imagename, *headerfile;
    char *ospsfname, *coeffname, *fittedpsfname, *imageweightname;
    char *firstname, *badname, **shiftsname, **starcubefile, **psfcubefile, **rescubefile, **fracrescubefile, **refpsfcubefile, **ellname, *lockfilename;
    char **argkey, **argval;
    char *pstr, delimiter[5], headdelimiter[5], *satname, *weight_suffix;

    double *chisq; // model chisq array
    double bestamp, chisqmin, pchisq, deltachisq, previousdelta;

    // Stuff for PSF ellipticities
    double modelpsfe[2], modelcentroid[2], starpsfe[2], starcentroid[2], astrometrycentroid[2];
    double starmeanmoments[3],modelmeanmoments[3];

    int iter;
    int niter=200;
    int final_nitermin;

    fitsfile *starcubepointer, *psfcubepointer, *rescubepointer, *refpsfcubepointer, *fracrescubepointer; 
    FILE *shiftsfile, *filep, *ellfile, *refpsfdata, *lockfile;

#ifdef KIDS
    char delims[] = "_";
    char delims2[] = "O";
#endif
#ifdef CFHT
    char delims[] = "_C";
    char delims2[] = "_C";
#endif
#ifdef SUPRIME
    char delims[] = "_O";
    char delims2[] = "_O";
#endif
#ifdef SINGLE
    char delims[] = "_C";
    char delims2[] = "_C";
#endif
    char *item = NULL;

    // GSL random number generator 
    T = gsl_rng_ranlux;
    r = gsl_rng_alloc(T);
    seedval = random_seed();
    gsl_rng_set(r,seedval);

    /*
     * variables for output fits file of stars used 
     */
    int status = 0;
    /*
     * table to identify each object 
     */
    char dotdelimiter[] = { "." };

    /*
     * swarp wcs variables 
     */
    catstruct *rawcat;
    rawcat = (catstruct *) NULL;
    fieldstruct *rawfield;
    rawfield = (fieldstruct *) NULL;
    wcsstruct *wcs_raw;
    wcs_raw = (wcsstruct *) NULL;
    ikernelstruct *ikernel;

    t1 = time(NULL);

    printf("\n %s development version ",argv[0]);
    char datedelims[]="()";
    char *dstring = NULL;
    char date_string[200];
    bzero(date_string,200);
    strcpy(date_string, DATE_STRING);
    dstring  = strtok(date_string, datedelims);
    if (dstring != NULL) dstring = strtok(NULL, datedelims);
    if (dstring != NULL) printf(" %s ",dstring);
    printf("\n");

    if (argc != 5 && argc != 7) {
        printf
            (" %s <imagelist> <global order> <chip order> <snratio> [<bright mag> <faint mag>] \n",
             argv[0]);
        exit(EXIT_FAILURE);
    }

    if (USE_SWARP == 1) {
        printf(" swarp astrometric corrections will be applied \n");
        if (WCS == 0) {
            printf
                (" WARNING:  WCS flag not set, using pixel object coordinates \n");
            printf
                (" even though WCS information is needed for the swarp correction \n");
        }

    if (CORRECT_DISTORTION == 1) {
      printf(" distortion corrections will be applied to the PSF\n");
      correct_distortion = 1;
        if (USE_SWARP == 0 || WCS == 0) {
            fprintf(stderr,
                    " CORRECT_DISTORTION is set but USE_SWARP or WCS are not\n");
            exit(EXIT_FAILURE);
        }
    } else {
      printf(" distortion corrections will not be applied to the PSF\n");
	correct_distortion = 0;
    }
        /*
         * read the swarp preferences file 
         */
        narg = 0;
        naxis = 2;
        QMALLOC(argkey, char *, 1);
        QMALLOC(argval, char *, 1);
        prefsname = getenv("SWARP_CONFIG");
        if (prefsname != NULL) {
            strcpy(prefs.prefs_name, prefsname);
            readprefs(prefs.prefs_name, argkey, argval, narg);
            prefs.verbose_type = QUIET; // turn off verbose reporting
            useprefs();
        } else {
            fprintf(stderr,
                    "configuration file environment variable SWARP_CONFIG not set \n");
            exit(EXIT_FAILURE);
        }
        free(argkey);
        free(argval);
        /*
         * initialise the convolution kernel 
         */
        ikernel = init_ikernel(prefs.resamp_type, naxis);
        for(i = 0; i < 2; i++)
            kern[i] = calloc(INTERP_MAXKERNELWIDTH, sizeof(double));
    } else {
        printf(" no swarp astrometric corrections will be applied \n");
    }
    /*
     * initialise the scalefactor 
     */
    *scalefactor = 0.;

    /*
     * test whether magnitude limits for PSF stars have been specified and read values
     * if so, or else set default values.  NB PSF star magnitude values need to be
     * given in the input PSF catalogue in the first case! 
     */

    if (argc == 8) {
        brightmag = atof(argv[6]);
        faintmag = atof(argv[7]);
        printf(" selecting PSF stars with %8.2f < m < %8.2f \n", brightmag,
               faintmag);
    } else {
        brightmag = -1.e10;
        faintmag = 1.e10;
    }

    order = atoi(argv[2]);

    if (strcmp(argv[3], "none") == 0) {
        chipvariation = 0;
        chiporder = -1;
        printf(" no chip variation allowed \n");
    } else {
        chipvariation = 1;
        chiporder = atoi(argv[3]);
        printf(" chip variation order %d \n", chiporder);
    }

    snratio = atof(argv[4]);
    if (snratio <= 0.)
        snratio = 0.;
    printf(" selecting stars with peak S/N ratio > %8.1f \n", snratio);
    if (snratio < 20.)
        printf(" WARNING: low S/N ratio not advised, choose S/N > 20 \n");

    bval = 0.;
    // printf(" input non-linearity parameter bval = %f \n", bval);

    /*
    if (WCS == 1) {
        printf(" reading WCS coordinates from input catalogue\n");
    } else {
        printf(" reading xy coordinates from input catalogue\n");
    }
    */

    /*
     * set some dimensions 
     */

    /*size of each individual postage stamp (data and psf) */
    if ((int)GREAT3 < 1)
      pwidth = 32; // make usual choice 32 (this can/will be expanded later in lensfit by padding)
                   // and a smaller size here speeds things up
    else
      pwidth = 48; // set postage stamp size to 48 for GREAT3

    pheight = pwidth;           /* postage stamps must be square */
    halfpwidth = pwidth / 2 + 1;
    halfpheight = pheight / 2;
    psfsize = pheight * halfpwidth;
    padheight = 50 * pheight;
    padwidth = 50 * pwidth;
    halfpadwidth = 1 + padwidth / 2;
    cchwidth = pwidth / 4;

    if (padwidth <= pwidth || padheight <= pheight) {
        printf
            (" padded cross-correlation array size must be bigger than input size \n");
        exit(EXIT_FAILURE);
    }

    /*
     * initialise memory counters 
     */
    memory = 0.;
    memloop = 0;

    xfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "xfit");
    yfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "yfit");
    zfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "zfit");
    wfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "wfit");
    zsort = ml_calloc(pwidth * pheight, sizeof(double), &memory, "wfit");
    objpix = ml_calloc(pwidth * pheight, sizeof(double), &memory, "objpix");

    /*
     * set tolerance - 
     * (i) if any other object lies within this radius
     * the object will be rejected inside extractpostagestamp/swarpextract.
     * (ii) if any bad pixels lie within this radius the star will be rejected
     * (iii) star shifts larger than this amount will cause stars to be
     * rejected
     */

    // tolerance radius in pixels
    poserror = 6;

    evariable = ml_calloc(PATH_MAX, sizeof(char), &memory, "evariable");
    headevariable = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "headevariable");
    badpixdir = ml_calloc(PATH_MAX, sizeof(char), &memory, "badpixdir");
    psfdir = ml_calloc(PATH_MAX, sizeof(char), &memory, "psfdir");
    prefsname = ml_calloc(PATH_MAX, sizeof(char), &memory, "prefsname");
    catname = ml_calloc(PATH_MAX, sizeof(char), &memory, "catname");
    imagename = ml_calloc(PATH_MAX, sizeof(char), &memory, "imagename");
    headerfile = ml_calloc(PATH_MAX, sizeof(char), &memory, "headerfile");
    fittedpsfname = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "fittedpsfname");
    firstname = ml_calloc(PATH_MAX, sizeof(char), &memory, "firstname");
    lockfilename = ml_calloc(PATH_MAX, sizeof(char), &memory, "lockfilename");
    imageweightname = ml_calloc(PATH_MAX, sizeof(char), &memory,
                                "imageweightname");
    badname = ml_calloc(PATH_MAX, sizeof(char), &memory, "badname");
    weight_suffix = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "weigh_suffix");
    satname = ml_calloc(PATH_MAX, sizeof(char), &memory, "satname");

    /*
     * create output PSF filenames
     */

    psfdir = getenv("PSF_DIR");
    if ((int)GREAT3 < 1)
      strcpy(firstname, argv[1]);
    else
      strcpy(firstname, "starfield_image-000-0");
    // PSF coefficients file at native sampling, old lensfit format
    coeffname = create_psf_name(psfdir, firstname, dotdelimiter);
    // PSF coefficients file with oversampling, new format
    if ((int)GREAT3 < 1)
      strcpy(firstname, argv[1]);
    else
      strcpy(firstname, "starfield_image-000-0");
    ospsfname = create_ospsf_name(psfdir, firstname, dotdelimiter);

    // write a lock file to the specified location
    lockfilename = getenv("LOCKFILENAME");
    if (lockfilename != NULL) 
      {
	if ((lockfile = fopen(lockfilename, "w")) == NULL) 
	  {
	    fprintf(stderr, " failed to make lockfile %s \n", lockfilename);
	    exit(EXIT_FAILURE);
	  }
	fprintf(lockfile," %s \n",firstname);
	fclose(lockfile);
	printf(" created lock file %s \n",lockfilename);
      }

        /*
         * define the string delimiter used to allow the code to
         * create the filename for the .head file and weight file, 
	 * in the case where
         * USE_SWARP is defined to be 1.  The code will take the last
         * instance of this delimiter in the input filename string,
         * strip off the delimiter and following characters, and add
         * on the head file suffix defined in the swarp configuration
         * file (normally ".head").  e.g.  715231p_18C.sub.fits will
         * have head filename 715231p_18.head if delimiter is set to
         * "C".  If delimiter is set to "." the expected head file
         * would be 715231p_18C.sub.head The code will exit if this
         * head file is not found in the same path as the data file.
         */

    if (USE_SWARP == 1)
      {
#ifdef KIDS
        strcpy(weight_suffix, "S.weight.fits");
	strcpy(delimiter,"S");
	strcpy(headdelimiter,"O");
#endif
#ifdef CFHT
        strcpy(weight_suffix, "C.weight.fits");
	strcpy(delimiter,"C");
	strcpy(headdelimiter,"C");
#endif
#ifdef SINGLE
        strcpy(weight_suffix, "C.weight.fits");
	strcpy(delimiter,"C");
	strcpy(headdelimiter,"C");
#endif
#ifdef SUPRIME
	strcpy(weight_suffix, "OFCS.weight.fits");
	strcpy(delimiter,"O");
	strcpy(headdelimiter,"O");
#endif
      }
    else
      {
	strcpy(weight_suffix, ".weight.fits");
	strcpy(delimiter,".");
      }

    /*
     * set up the chip mosaic geometry
     */
#ifdef KIDS
    xchipsampling = 2040 + 70;
    ychipsampling = 4090 + 70;
    big_gap = 0;
    nxchip = 8;
    nychip = 4;
#endif
#ifdef CFHT
    xchipsampling = 2048 + 70;
    ychipsampling = 4612 + 70;
    big_gap = 425 - 70;
    nxchip = 9;
    nychip = 4;
#endif
#ifdef SUPRIME
    xchipsampling = 2013 + 70;  // total x size covered by chip plus gap
    ychipsampling = 4085 + 70;  // ditto for y
    big_gap = 0;                // allowance for extra gaps between some rows
    nxchip = 5;                 // number of chips on x axis
    nychip = 2;                 // number of chips on y axis
#endif
#ifdef SINGLE
  /* single chip geometry */
  xchipsampling = 4096;  // set this to the x size of the single image
  ychipsampling = 4096;  // same for y size
  big_gap = 0; // no gap between chips (because there's only one)
  nxchip = 1;  // one chip on x-axis
  nychip = 1;  // one chip on y-axis
#endif

    nchip_expected = nxchip * nychip;   // total number of chips
    hxsize = (double) (nxchip * xchipsampling) / 2.;    // half the total pixels on x axis
    hysize = (double) (nychip * ychipsampling + 2 * big_gap) / 2.;      // half the total pixels on y axis

    image_file = calloc(nchip_expected, sizeof(char *));

    /*
     * read the input file list 
     */
    filep = fopen(argv[1], "r");
    if (filep == NULL) {
        fprintf(stderr, " Error opening file %s \n", argv[1]);
        exit(EXIT_FAILURE);
    }

    nchip = 0;
    while (!feof(filep) && !ferror(filep)) {
        if (fscanf(filep, "%s", firstname) != 0) {
            if (!feof(filep) && !ferror(filep)) {
                if (nchip >= nchip_expected) {
                    fprintf(stderr, " too many images in file \n");
                    exit(EXIT_FAILURE);
                }
                image_file[nchip] = calloc(PATH_MAX, sizeof(char));
                strcpy(image_file[nchip], firstname);
                nchip++;
            }
        }
    }
    fclose(filep);
    status = 0;
    bzero(firstname, PATH_MAX);
    printf(" filenames read for %d input images \n", nchip);

    if (nchip<=0)
      {
	fflush(stdout);
	fprintf(stderr," no input images read from input list \n");
	exit(2);
      }

    if (nchip<nchip_expected)
      {
	printf(" WARNING: %d chips read from list but %d were expected \n",nchip,nchip_expected);
	fflush(stdout);
      }

    /*
     * get the order of the polynomial fit.  negative order is taken to mean
     * that cross-terms should be included 
     */
    get_poly_order(order, chiporder, chipvariation, nchip, &ncoeffs,
                   &crossterm);


    /*
     * define order of fit functions for x, y shifts 
     */

    sorder = 1;
    schipvariation = 1;
    schiporder = 1;

    if (sorder >= 0) {
        /*
         * number of coeffs for global fit 
         */
        sncoeffs = (1 + sorder) * (2 + sorder) / 2;
        /*
         * number of extra coeffs for local fits of 3 coeffs each 
         */
        if (schipvariation == 1) {
            if (schiporder > 0)
                sncoeffs +=
                    (1 + schiporder) * (2 + schiporder) * (nchip - 1) / 2;
            else
                sncoeffs += (1 - schiporder) * (1 - schiporder) * (nchip - 1);
        }
        /*
         * no global crossterms 
         */
        scrossterm = 0;
        printf
            (" fitting polynomial surface to shifts of order %d and %d coefficients \n",
             sorder, sncoeffs);
        if (schipvariation == 1)
            printf(" including chip-dependent terms of order %d \n",
                   schiporder);
    } else {
        sorder = fabs(sorder);
        sncoeffs = (1 + sorder) * (1 + sorder);
        /*
         * number of extra coeffs for local fits of 3 coeffs each 
         */
        if (schipvariation == 1) {
            if (schiporder > 0)
                sncoeffs +=
                    (1 + schiporder) * (2 + schiporder) * (nchip - 1) / 2;
            else
                sncoeffs += (1 - schiporder) * (1 - schiporder) * (nchip - 1);
        }
        /*
         * allow global crossterms 
         */
        scrossterm = 1;
        printf
            (" fitting polynomial surface to shifts of order %d with global cross-terms and %d coefficients \n",
             sorder, sncoeffs);
        if (schipvariation == 1)
            printf(" including chip-dependent terms of order %d \n",
                   schiporder);
    }

    if (sncoeffs >= 500) {
        fprintf(stderr, " total number of astrometry coefficients is too high \n");
        exit(2);
    }

    maxncoeffs = ncoeffs;
    if (sncoeffs > maxncoeffs) maxncoeffs = sncoeffs;
    if (maxncoeffs < 6)
        maxncoeffs = 6;

    u = ml_calloc((1 + maxncoeffs), sizeof(double *), &memory, "u");
    for(i = 0; i <= maxncoeffs; i++)
        u[i] = ml_calloc((1 + maxncoeffs), sizeof(double), &memory, "u[i]");

    v = (double **) calloc((1 + maxncoeffs), sizeof(double *));
    memory += (float) (1 + maxncoeffs) * sizeof(double);
    for(i = 0; i <= maxncoeffs; i++) {
        memory += (float) (1 + maxncoeffs) * sizeof(double);
        v[i] = (double *) calloc((1 + maxncoeffs), sizeof(double));
    }

    xcoeffs = (double *) calloc(maxncoeffs, sizeof(double));
    ycoeffs = (double *) calloc(maxncoeffs, sizeof(double));
    memory += 2. *maxncoeffs * sizeof(double);

    avals = (double *) calloc(1 + maxncoeffs, sizeof(double));
    w = (double *) calloc(1 + maxncoeffs, sizeof(double));
    memory += 2. * (1 + maxncoeffs) * sizeof(double);

    acoeffs = (double **) calloc((pwidth * pheight), sizeof(double *));
    memory += (float) (pwidth * pheight) * sizeof(double *);
    if (acoeffs == NULL) {
        printf("Memory allocation error for sub-image pointers\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < (pwidth * pheight); i++) {
        acoeffs[i] = (double *) calloc((ncoeffs), sizeof(double));
        memory += (float) (ncoeffs) * sizeof(double);
        if (acoeffs[i] == NULL) {
            printf("Memory allocation error for sub-image \n");
            exit(EXIT_FAILURE);
        }
    }

    /*
     * temp array is used as spare storage inside extractdata 
     * for quadrant swapping when the postage stamps are extracted 
     */
    temp = ml_calloc(pwidth * pheight, sizeof(float), &memory, "temp");
    badtemp = ml_calloc(pwidth * pheight, sizeof(float), &memory, "badtemp");

    /*
     * get environment variables holding directory names 
     */

    evariable = getenv("DATA_DIR");
    headevariable = getenv("HEAD_DIR");
    badpixdir = getenv("BADPIX_DIR");

    /*
     * create the psfstars catalogue filename 
     */
    catname = getenv("CATALOGUE_STARS");
    if (catname == NULL) {
        fprintf(stderr,
                " catalogue name not set, use environment variable CATALOGUE_STARS \n");
        exit(EXIT_FAILURE);
    }

    if (access(catname, F_OK) != 0) {
        printf(" Can't read catalogue %s \n", catname);
        exit(EXIT_FAILURE);
    }

    /*
     * Read in reference positions (x,y on a chip) at which PSF model should be output
     */
    ncount = nrefpsf = 0;
    if ( access("psfrefpos.asc", F_OK) != 0 ) {
        printf(" Can't read catalogue of PSF reference positions: psfrefpos.asc \n");
    } else {
      nrefpsf = readpsfcatsize("psfrefpos.asc"); // I guess this will work for the reference PSF catalogue format of x, y as well
      refpsfx = ml_calloc(nrefpsf, sizeof(double), &memory, "refpsfx");
      refpsfy = ml_calloc(nrefpsf, sizeof(double), &memory, "refpsfy");
      refpsfdata = fopen("psfrefpos.asc", "r");
      if (refpsfdata == NULL) {
        printf(" failed to read file psfrefpos.asc \n");
        exit(2);
      }
      while (fgets(readstring, 2000, refpsfdata) != NULL) {
        if (strncmp(readstring, "#", 1) != 0) {
	  readitem = strtok(readstring, readdelims);
	  if ( readitem == NULL ) {
	    fprintf(stderr," error reading 1st column (x) from psfrespos.asc \n");
	    exit(2);
	  }
	  refpsfx[ncount] = atof(readitem);
	  readitem = strtok(NULL, readdelims);
	  if ( readitem == NULL ) {
	    fprintf(stderr," error reading 2nd column (y) from psfrespos.asc \n");
	    exit(2);
	  }
	  refpsfy[ncount] = atof(readitem);
	  ncount++;
        }
      }
      fclose(refpsfdata);
      if (ncount != nrefpsf) {
        printf(" error reading catalogue psfrespos.asc \n");
        printf(" nrefpsf = %d, ncount = %d \n", nrefpsf, ncount);
        exit(EXIT_FAILURE);
      }
    }

    /*
     * read catalogue of objects 
     */

    nobjt = readpsfcatsize(catname);

    objx = ml_calloc(nobjt, sizeof(float), &memory, "objx");
    objy = ml_calloc(nobjt, sizeof(float), &memory, "objy");
    mag = ml_calloc(nobjt, sizeof(float), &memory, "mag");

    if (WCS == 0) {
        fprintf(stderr,
                " x,y star coordinate input not supported in this version \n");
        exit(EXIT_FAILURE);
    } else {
        /*
         * case where world coordinates have been supplied 
         */
        printf(" reading world coordinates from input catalogue\n");

        dobjx = ml_calloc(nobjt, sizeof(double), &memory, "dobjx");
        dobjy = ml_calloc(nobjt, sizeof(double), &memory, "dobjy");
        radegs = ml_calloc(nobjt, sizeof(double), &memory, "radegs");
        decdegs = ml_calloc(nobjt, sizeof(double), &memory, "decdegs");
        rdegs = ml_calloc(nobjt, sizeof(double), &memory, "radegs");
        ddegs = ml_calloc(nobjt, sizeof(double), &memory, "decdegs");
        offscale = ml_calloc(nobjt, sizeof(int), &memory, "offscale");
        rmag = ml_calloc(nobjt, sizeof(float), &memory, "rmag");

        /*
         * read in coordinates assuming in decimal WCS 
         */
        nobj2 = readpsfcat_simple_wcs(catname, radegs, decdegs, rmag);
        if (nobj2 != nobjt) {
            printf(" error reading catalogue \n");
            printf(" nobjt = %d, nobj2 = %d \n", nobjt, nobj2);
            exit(EXIT_FAILURE);
        }

        if (nobjt > 0) {
            /*
             * check coords are in sensible range 
             */
            for(i = 0; i < nobjt; i++) {
                if (radegs[i] < 0. || radegs[i] > 360. || decdegs[i] < -90.
                    || decdegs[i] > 90.) {
                    fprintf(stderr,
                            " these don't look like world coordinates \n");
                    fprintf(stderr, " %lf %lf \n", radegs[i], decdegs[i]);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    if (nobjt > 0) {
        printf(" %d objects read from PSF catalogue %s \n", nobjt, catname);
    } else {
      fflush(stdout);
      fprintf(stderr," no objects read from PSF catalogue %s \n", catname);
      exit(EXIT_FAILURE);
    }


    /*
     * memory for fitting 
     */

    chips = (int *) calloc(nobjt, sizeof(int));
    wfits = (double *) calloc(nobjt, sizeof(double));
    ow = (double *) calloc(nobjt, sizeof(double));
    xfits = (double *) calloc(nobjt, sizeof(double));
    yfits = (double *) calloc(nobjt, sizeof(double));
    xshiftfits = (double *) calloc(nobjt, sizeof(double));
    yshiftfits = (double *) calloc(nobjt, sizeof(double));
    bxshiftfits = (double *) calloc(nobjt, sizeof(double));
    byshiftfits = (double *) calloc(nobjt, sizeof(double));
    rx = (double *) calloc(nobjt, sizeof(double));
    ry = (double *) calloc(nobjt, sizeof(double));
    if (chips == NULL || wfits == NULL || xfits == NULL || yfits == NULL
        || xshiftfits == NULL || yshiftfits == NULL || rx == NULL || ry == NULL) {
        fflush(stdout);
        fprintf(stderr, " memory allocation error, fits arrays \n");
        fflush(stderr);
        exit(2);
    }

    /*
     * create the filename for the psf catalogue file 
     */

    strcpy(firstname, "starfield_image-000-0");
    if (!(pstr = strrchr(firstname, *dotdelimiter)))
        pstr = firstname + strlen(firstname);
    sprintf(pstr, "%s", ".surfacepsf.fits");

    if (psfdir != NULL) {
        strcpy(fittedpsfname, psfdir);
        len = strlen(fittedpsfname);
        if (strncmp(&fittedpsfname[len - 1], "/", 1) != 0) {
            strncat(fittedpsfname, "/", 1);
        }
        if (access(fittedpsfname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr, " Can't write psf files to %s \n", fittedpsfname);
	  exit(EXIT_FAILURE);
        }
        len = strlen(firstname);
        strncat(fittedpsfname, firstname, len);
    } else {
        strcpy(fittedpsfname, firstname);
    }

    //      printf(" creating psf name %s \n",fittedpsfname);

    if (access(fittedpsfname, F_OK) == 0) {
        printf(" *** PSF file exists: removing \n");
        remove(fittedpsfname);
    }


    numstar = calloc(nchip, sizeof(int));
    if (numstar == NULL) {
      fflush(stdout);
      fprintf(stderr, " error allocating memory for numstar \n");
      exit(EXIT_FAILURE);
    }

    /*
     * got through each input image and accumulate the star postage stamps etc 
     */

    tot_nstar = 0;

    /*
     * define a reference image for the scale factor from the middle
     * of the field. The scalefactor for all the images in this
     * exposure will be set from this reference.  Note that other
     * exposures being processed elsewhere will have a different
     * scalefactor
     */
    refchip = nxchip / 2 + (nychip / 2) * nxchip;
    // set a default halfway through the list (shouldn't be needed as
    // long as the ref chip is in the list)
    refimage = nchip / 2;
    // if more than one chip, then try to find the reference chip 
    if (nchip>1)
      {
	for(image = 0; image < nchip; image++) {
	  /*
	   * strip out the chip number and look for the refimage chip number 
	   */
	  bzero(firstname, PATH_MAX);
	  strcpy(firstname, image_file[image]);
	  item = strtok(firstname, delims);
	  item = strtok(NULL, delims2);
	  i = atoi(item) - 1;
	  if (i == refchip) {
            refimage = image;
            break;
	  }
	}
      }

    nstar = ml_calloc(nchip, sizeof(int), &memloop, "nstar array");
    dpix = ml_calloc(nchip, sizeof(double **), &memloop,
			"Memory allocation error for sub-image pointers");
    ppix = ml_calloc(nchip, sizeof(double **), &memloop,
			"Memory allocation error for sub-image pointers");
    dbadpix = ml_calloc(nchip, sizeof(double **), &memloop,
			"Memory allocation error for sub-image pointers");    
    sn = ml_calloc(nchip, sizeof(float*), &memory, "sn pointers");
    starweight = ml_calloc(nchip, sizeof(double*), &memory, "starweight pointers");
    chisqval = ml_calloc(nchip, sizeof(double*), &memory, "chisqval pointers");
    // pointers to arrays of global x,y coords of objects on each chip
    xval = ml_calloc(nchip, sizeof(double*), &memory, "xval pointers");
    yval = ml_calloc(nchip, sizeof(double*), &memory, "yval pointers");
    // pointers to arrays of local x,y coords of objects on each chip
    chipxval = ml_calloc(nchip, sizeof(double*), &memory, "xval pointers");
    chipyval = ml_calloc(nchip, sizeof(double*), &memory, "yval pointers");
    // array of x,y midpoints of each chip
    xchipval = ml_calloc(nchip, sizeof(double), &memory, "xchipval");
    ychipval = ml_calloc(nchip, sizeof(double), &memory, "ychipval");
    // array of RA, dec values for each chip
    chiprdegs = ml_calloc(nchip, sizeof(double), &memory, "chiprdegs");
    chipddegs = ml_calloc(nchip, sizeof(double), &memory, "chipddegs");

    starcubefile = ml_calloc(nchip, sizeof(char*), &memory, "starcubefile");
    psfcubefile = ml_calloc(nchip, sizeof(char*), &memory, "psfcubefile");
    rescubefile = ml_calloc(nchip, sizeof(char*), &memory, "rescubefile");
    fracrescubefile = ml_calloc(nchip, sizeof(char*), &memory, "fracrescubefile");
    refpsfcubefile = ml_calloc(nchip, sizeof(char*), &memory, "refpsfcubefile");
    shiftsname = ml_calloc(nchip, sizeof(char*), &memory, "shiftsname");
    ellname = ml_calloc(nchip, sizeof(char*), &memory, "ellname");

    noise = ml_calloc(nchip, sizeof(float), &memory, "noise[]");
    gain = ml_calloc(nchip, sizeof(float), &memory, "gain[]");

    for (kimage = 0; kimage < nchip; kimage++) {
        // start with the reference image (only relevant if distortion correction enabled)
        image = kimage + refimage;
        if (image >= nchip) image -= nchip;

        /*
         * strip out the chip number and remember it 
         */
        bzero(firstname, PATH_MAX);
        strcpy(firstname, image_file[image]);

	// work out the chip ID 
	if (nchip>1)
	  {
	    // if multiple chips, attempt to read the chip number from the filename
	    item = strtok(firstname, delims);
	    item = strtok(NULL, delims2);
	    i = atoi(item) - 1;
	    chipnumber = i;
	  }
	else
	  {
	    // only one chip expected, so set chip ID to 0
	    chipnumber = i = 0;
	  }
        ychip = i / nxchip;
        xchip = i - ychip * nxchip;
        printf("\n image %s chip %d %d \n", image_file[image], xchip + 1,
               ychip + 1);
        if (i < 0 || i >= nchip_expected) {
	  fflush(stdout);
	  fprintf(stderr," error reading chip number for image %s \n",image_file[image]);
	  fprintf(stderr," chip number %d out of %d x %d \n",i+1,nxchip,nychip);
	  exit(EXIT_FAILURE);
        }
	nstar[chipnumber] = 0;

        /*
         * strip off .fits extension on first argument if supplied 
         */
        len = strlen(image_file[image]);
        if (len <= 0) {
            fprintf(stderr, " error reading arguments from command line \n");
            exit(EXIT_FAILURE);
        }
        if (strncmp(image_file[image] + len - 5, ".fits", 5) == 0
            || strncmp(image_file[image] + len - 5, ".FITS", 5) == 0) {
            len = len - 5;
        }

        bzero(firstname, PATH_MAX);
        strncpy(firstname, image_file[image], len);

	shiftsname[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "shiftsname");
        bzero(shiftsname[chipnumber], PATH_MAX);
        strcpy(shiftsname[chipnumber], firstname);
        strcat(shiftsname[chipnumber], "_shifts.log");

	ellname[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "ellname");
        bzero(ellname[chipnumber], PATH_MAX);
        strcpy(ellname[chipnumber], firstname);
        strcat(ellname[chipnumber], "_ellipticities.log");

        starcubefile[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "starcubefile[chipnumber]");
        bzero(starcubefile[chipnumber], PATH_MAX);
        strcpy(starcubefile[chipnumber], firstname);
        strcat(starcubefile[chipnumber], "_stars.fits");

	psfcubefile[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "psfcubefile[chipnumber]");
        bzero(psfcubefile[chipnumber], PATH_MAX);
        strcpy(psfcubefile[chipnumber], firstname);
        strcat(psfcubefile[chipnumber], "_psf.fits");

        rescubefile[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "rescubefile[chipnumber]");
        bzero(rescubefile[chipnumber], PATH_MAX);
        strcpy(rescubefile[chipnumber], firstname);
        strcat(rescubefile[chipnumber], "_residuals.modelamp.fits");

        fracrescubefile[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "fracrescubefile[chipnumber]");
        bzero(fracrescubefile[chipnumber], PATH_MAX);
        strcpy(fracrescubefile[chipnumber], firstname);
        strcat(fracrescubefile[chipnumber], "_fracresiduals.fits");
        
        refpsfcubefile[chipnumber] = ml_calloc(PATH_MAX, sizeof(char), &memory, "refpsfcubefile[chipnumber]");
        bzero(refpsfcubefile[chipnumber], PATH_MAX);
        strcpy(refpsfcubefile[chipnumber], firstname);
        strcat(refpsfcubefile[chipnumber], "_refpsf.fits");

        bzero(imageweightname, PATH_MAX);

        /*
         * get the image weight filename 
         */
        if (badpixdir != NULL) {
            // printf(" Path to good/bad pixel images %s \n",badpixdir);
            strcpy(imageweightname, badpixdir);
            len = strlen(imageweightname);
            if (strncmp(&imageweightname[len - 1], "/", 1) != 0) {
                strncat(imageweightname, "/", 1);
                len++;
            }
            len = strlen(image_file[image]);
            strncat(imageweightname, image_file[image], len);
        } else {
            strcpy(imageweightname, image_file[image]);
        }

        /*
         * strip off end of input image name so that weight suffix can be added 
         */
        if (!(pstr = strrchr(imageweightname, *delimiter)))
            pstr = imageweightname + strlen(imageweightname);
        sprintf(pstr, "%s", weight_suffix);

        if (GZIPPED > 0) {
            strcat(imageweightname, ".gz");
        }

        printf(" weight file %s \n", imageweightname);

        if (evariable != NULL) {
            strcpy(imagename, evariable);
            len = strlen(imagename);
            if (strncmp(&imagename[len - 1], "/", 1) != 0) {
                strncat(imagename, "/", 1);
            }
            len = strlen(firstname);
            strncat(imagename, firstname, len);
        } else {
            strcpy(imagename, firstname);
        }
        strncat(imagename, ".fits", 5);
	if (VERBOSE == 1)
	  printf(" reading image data from %s \n", imagename);

        /*
         * by default, look for the head files in the same place as
         * the data files but if the head environment variable is set,
         * use that location instead
         */

        bzero(headerfile, PATH_MAX);

        if (USE_SWARP == 1) {
            if (headevariable == NULL && evariable != NULL) {
                headevariable = evariable;
            }

            if (headevariable != NULL) {
                strcpy(headerfile, headevariable);
                len = strlen(headerfile);
                if (strncmp(&headerfile[len - 1], "/", 1) != 0) {
                    strncat(headerfile, "/", 1);
                }
                strcat(headerfile, image_file[image]);
            } else {
                strcpy(headerfile, image_file[image]);
            }
            strcpy(satname, prefs.sat_keyword);
        } else {
            strcpy(satname, "SATLEVEL");
        }

        /*
         * use fitsio routines to get the image info and check the
         * "bad ccd" keyword (we could do this with the swarp fits
         * routines instead but leave this as it is just now).
	 Note gain & satlev are assumed to be the same for all input images.
         */
        if (getimagesize
            (imagename, satname, dim, &gain[chipnumber], &satlev, &arcperpix, &angle,
             &badccd, &imageid, (float) minexposuretime)) {
            fflush(stdout);
            fprintf(stderr, " error checking input image %s \n", imagename);
            exit(EXIT_FAILURE);
        }

	if (imageid <= 0) imageid = 1+xchip+ychip*nxchip;

        if (imageid != (1 + xchip + ychip * nxchip)) {
            fflush(stdout);
            fprintf(stderr, "error in identifying chip number, %s. \n",
                    imagename);
            fflush(stderr);
            exit(EXIT_FAILURE);
        }

        if (badccd == 1) {
            printf(" %s flagged as bad, skipping this image \n", imagename);
            continue;
        }

	// define centres of chips in global coordinates (used in randomisesoln)
	xchipval[chipnumber] = (dim[0]/2 + xchip * xchipsampling) / hxsize - 1.;
	ychipval[chipnumber] = dim[1]/2 + ychip * ychipsampling;
	if (ychip > 0)
	  ychipval[chipnumber] += big_gap;
	if (ychip > 2)
	  ychipval[chipnumber] += big_gap;
	ychipval[chipnumber] = ychipval[chipnumber]/ hysize - 1.;
		      	

        if (USE_SWARP == 1) {
            /*
             * create headerfile name needed for the swarp code 
             */
            if (!(pstr = strrchr(headerfile, *headdelimiter)))
                pstr = headerfile + strlen(headerfile);
            sprintf(pstr, "%s", prefs.head_suffix);
            /*
             * test it can be opened 
             */
            filep = fopen(headerfile, "r");
            if (filep == NULL) {
                fprintf(stderr, "error opening scamp file %s \n", headerfile);
                printf(" error opening scamp file %s \n", headerfile);
                printf
                    (" without these files no further astrometric correction can be applied \n");
                printf
                    (" and yet USE_SWARP correction has been specified in the code \n");
                printf
                    (" either define USE_SWARP 0 or supply a set of scamp head files, one per image \n");
                fflush(stdout);
                exit(EXIT_FAILURE);
            }
            fclose(filep);
        }

        imagesize = dim[0] * dim[1];

	if (VERBOSE == 1)
	  {
	    printf(" image dimensions = %d x %d \n", dim[0], dim[1]);
	    printf (" CCD gain = %6.1f e/ADU \n",gain[chipnumber]);
	    printf(" CCD saturation level = %8.0f ADU \n", satlev);
	  }

        /*
         * set maximum pixel value as 50 percent of CCD saturation level
         * in case of CCD non-linearity 
         */
        satlimit = 0.5 * satlev;

	/* set the maximum fraction of the total flux which is allowed to be in the 
	   maximum pixel - designed to cut out any cosmic rays that have made it through the
	   various filtering processes */
	cosmicraylimit = 0.25;  // 0.5 is essentially no cut: for CFHTLenS the value was 0.15

        /*
         * set maximum intensity used for object detection 
         */
        fmax_intensity = satlimit;

	if (VERBOSE == 1)
	  printf(" trying to allocate %d MB for images and arrays\n",
		 (int) (15. * imagesize * 4. / 1024. / 1024.));

        region = ml_calloc(imagesize, sizeof(int), &memloop,
                           "Error allocating memory for region mask \n");

        /*
         * allocate memory for input image 
         */
        apix = ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for image \n");
        badpix = ml_calloc(imagesize, sizeof(float), &memloop,
                           " error allocating memory for image \n");

	if (SUBTRACT_MEDIAN == 1)
	  sortarray = ml_calloc(imagesize, sizeof(float), &memloop,
				" error allocating memory for sortarray\n");

	if (FILTER_WEIGHTS == 1)
	  {
	    opix = (float*)ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for image \n");
	    changed = (float*)ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for changed array \n");
	    weightfilterarray = (float*)ml_calloc(1+2*fmheight, sizeof(float), &memloop,
                         "Error allocating memory for weightfilterarray \n");
	  }

        /*
         * use swarp code to get the full wcs info 
         */
        rawcat = read_cat(imagename);
        QCALLOC(rawfield, fieldstruct, 1);
        load_dfield(rawfield, rawcat, headerfile, 0, 0);        // modified version of load_field
        wcs_raw = copy_wcs(rawfield->wcs);      // get the wcs info

        /*
         * read image data and measure noise level 
         */

         //in GREAT3 : re-calibrated the noise function returned by the following
        getdata_badpix(imagename, imageweightname, dim, apix, badpix, &noise[chipnumber]);

	/* deal with possible problems in the bad pixel mask by stripping out individual bad pixels
	   to generate a new bad pixel mask.  This will be used in conjunction with the original
	   mask when data are extracted around stars, in order to make sure that the peaks of
	   stars have not been flagged as bad */
	if ((int)FILTER_WEIGHTS == 1  && (int)GREAT3 < 1)
	  {
	    // create new bad pixel mask in which individual bad pixels are stripped out (but not
	    // if they are part of a bad column)
	    nchanged = weightfilter(badpix, opix, changed, weightfilterarray, dim, fmheight);
	    if (VERBOSE==1)
	      printf(" %d bad pixels changed to good \n",nchanged);
	    // dilate the new bad pixel mask by one pixel to return it to the original state minus
	    // the individual bad pixels
	    dilatemask(opix,dim);
	  }
	else
	  {
	    opix = badpix;
	    nchanged = 0;
	  }

	// dilate the original bad pixel mask by one pixel in all directions to take care of 
	// left-over unmasked parts (do this regardless of whether we have created the new
	// mask above, or not)
	if ((int)GREAT3 < 1)
	  dilatemask(badpix,dim);

	if (VERBOSE == 1)
	  {
	    printf(" estimated background noise in image = %f ADU \n", noise[chipnumber]);
	  }

        if (noise[chipnumber] <= 0.) {
	  fflush(stdout);
	  fprintf(stderr, " Erroneous noise value \n");
	  exit(EXIT_FAILURE);
        }

        /*
         * subtract assumed constant median background if specified
         */
        maxlevel = satlev / 2.;

        if (SUBTRACT_MEDIAN == 1) {
            mediansub(apix, badpix, dim, maxlevel, sortarray);
        }

        /*
         * identify objects within this image 
         */

        /*
         * set area limits 
         */
        lower_area_limit = 1;   // try to eliminate single pixels eg cosmic rays
        lower_merging_limit = 5;        // anything bigger than this will be kept unmerged
        arealimit = 100000;
        /*
         * set lower intensity threshold 
         */
        fintensity_limit = 3. * noise[chipnumber];
        /*
         * maximum object detection pixel value has already been set to half of
         * the sat level.  Now allocate arrays large enough to hold the 
         * digitised pixel values.  This must match with assumed array sizes
         * inside varylpthresholdf 
         */
        if (noise[chipnumber] > 0.) {
            scale = 2. / noise[chipnumber];
            imax_intensity = (int) (fmax_intensity * scale);
        } else {
            printf(" noise value < 0 \n");
            exit(EXIT_FAILURE);
        }

        if ((int)GREAT3 < 1)
	  {
	    nobjects = varylpthresholdf2(apix, region, noise[chipnumber], dim[0], dim[1],
                                     lower_area_limit, lower_merging_limit,
                                     arealimit, fintensity_limit,
                                     fmax_intensity, imax_intensity);
	    if (VERBOSE == 1)
	      printf(" %d objects identified within image \n", nobjects);
	  }

        xsize = (double) dim[0];
        ysize = (double) dim[1];
        /*
         * convert WCS to pixel positions using swarp WCS library 
         */
        for(ii = 0; ii < nobjt; ii++) {
            dobjx[ii] = -999.;
            dobjy[ii] = -999.;
            offscale[ii] = 1;
            wcspos[0] = radegs[ii];
            wcspos[1] = decdegs[ii];

	    // check if star exists on this chip (remove for GREAT3 as star positions are arbitrary)
            if (wcs_to_raw(wcs_raw, wcspos, rawpos) == RETURN_OK) {
	      if ((int)GREAT3 == 1) 
		{
		  // GREAT3 accept all stars
		  dobjx[ii] = rawpos[0];
		  dobjy[ii] = rawpos[1];
		  offscale[ii] = 0;
		} else {
		// otherwise check if star is on-chip
                if (rawpos[0] > 0. && rawpos[0] < (double) rawfield->width) {
                    if (rawpos[1] > 0. && rawpos[1] < (double) rawfield->height) {
                        dobjx[ii] = rawpos[0];
                        dobjy[ii] = rawpos[1];
                        offscale[ii] = 0;
		    }
		}
		//  printf(" pos check %lf %lf %lf %lf \n",radegs[ii],decdegs[ii],dobjx[ii],dobjy[ii]);
	      }
	    }
	}

        /*
         * eliminate objects off the image 
         */
        /*
         * nobjt is the total number of stars in the catalogue
         * nobj is the number of stars that lie on this particular chip 
         */
        nobj = 0;
        for(i = 0; i < nobjt; i++) {
            if (offscale[i] == 0) {
                rdegs[nobj] = radegs[i];
                ddegs[nobj] = decdegs[i];
                objx[nobj] = (float) dobjx[i];
                objy[nobj] = (float) dobjy[i];
                mag[nobj] = rmag[i];
                // printf(" %lf %lf %f %f \n",radegs[i],decdegs[i],objx[nobj],objy[nobj]);
                nobj++;
            } else {
                //printf(" %lf %lf off-image \n",radegs[i],decdegs[i]);
            }
        }

	// allocate memory for postage stamp images 
        dpix[chipnumber] = ml_calloc(nobj, sizeof(double *), &memloop,
                         "Memory allocation error for sub-image pointers\n");
        ppix[chipnumber] = ml_calloc(nobj, sizeof(double *), &memloop,
                         "Memory allocation error for sub-image pointers\n");
        for(i = 0; i < nobj; i++)
	  {
            dpix[chipnumber][i] = ml_calloc(pwidth * pheight, sizeof(double), &memloop,
                                "Memory allocation error for sub-image \n");
            ppix[chipnumber][i] = ml_calloc(pwidth * pheight, sizeof(double), &memloop,
                                "Memory allocation error for sub-image \n");
	  }

        dbadpix[chipnumber] = ml_calloc(nobj, sizeof(double *), &memloop,
                            "Memory allocation error for sub-image pointers\n");
        for(i = 0; i < nobj; i++)
            dbadpix[chipnumber][i] = ml_calloc(pwidth * pheight, sizeof(double), &memloop,
                                   "Memory allocation error for sub-image \n");

	sn[chipnumber] = ml_calloc(nobj, sizeof(float), &memory, "sn");
	starweight[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "starweight");
	chisqval[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "chisqval");
	xval[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "xval");
	yval[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "yval");
	chipxval[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "xval");
	chipyval[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "yval");
        chiprdegs[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "ra");
        chipddegs[chipnumber] = ml_calloc(nobj, sizeof(double), &memory, "dec");

        /*
         * set limits on portion of image that a PSF must occupy (excludes
         * a narrow band around the edges) 
         */
	if ((int)GREAT3 == 1)
	  {
	    // apply the GREAT3 postage stamp shifts
	    xmin = (pwidth / 2) - 2;
	    ymin = (pheight / 2) -2;
	    xmax = (dim[0] - pwidth / 2) + 5;
	    ymax = (dim[1] - pheight / 2) + 5;
	  } else {
	  // assume usual convention for postage stamp centering
	  xmin = pwidth / 2;
	  ymin = pheight / 2;
	  xmax = dim[0] - pwidth / 2;
	  ymax = dim[1] - pheight / 2;
	}

	if (VERBOSE == 1)
	  printf(" total memory allocated %d MB \n",
		 (int) ((memory + memloop) / 1024. / 1024.));

        /*
         * loop through each star from the catalogue, find which psf
         * it belongs to: then extract the data
         */

        t2 = time(NULL);
        //  printf (" Initial setup time was %g secs \n", difftime(t2, t1) );

	num_snratio = 0;
	num_saturated = 0;
	num_off = 0;
	num_wrong_mag = 0;
	num_close = 0;
	num_badpix = 0;

	for(i = 0; i < nobj; i++) {

	  close = 0;
	  if ((int)GREAT3 < 1)
	    {
	      // test if close to another star in PSF catalogue
	      for(d = 0; d < nobj; d++) {
		if (d != i) {
		  if (fabs(objx[i] - objx[d]) < pwidth / 2
		      && fabs(objy[i] - objy[d]) < pheight / 2) {
		    close = 1;  // has a close neighbour
		  }
		}
	      }
	      if (close == 1)
		num_close++;
	    }

	  // if not near another star in the input catalogue, extract its postage stamp
	  if (close == 0) {
	    
	    if ((int)GREAT3 == 1)
	      {
		// fudge the star positions for GREAT3
		x = (48.0*i + 24.0) + 1;
		y = (24.0) + 1;
	      } else {
	      // use the actual star positions (normal use)
	      x = (int) (objx[i] + 0.5) - 1;
	      y = (int) (objy[i] + 0.5) - 1;
	    }

	    if (x >= xmin && x < xmax && y >= ymin && y < ymax ) {
	      if (mag[i] >= brightmag && mag[i] <= faintmag) {
		ximageoffset = 0.;
		yimageoffset = 0.;
		maxvalue = 0.;
		maxlevel = satlev / 2.;
		if (USE_SWARP == 0)
		  ngoodpix = extractpostagestamp
		    (apix, badpix, dim, ppix[chipnumber][nstar[chipnumber]], dbadpix[chipnumber][nstar[chipnumber]],
		     temp, badtemp, objx[i], objy[i],
		     poserror, ximageoffset, yimageoffset,
		     pwidth, pheight, noise[chipnumber], region,
		     maxlevel);
		else {
		  /*
		   * if this is the first time
		   * swarpextract has been called, the
		   * scalefactor will be zero, so report
		   * the current image.  Thereafter the
		   * scalefactor will be set to the
		   * central value in this image for all
		   * subsequent calls to swarpextract
		   */
		  if (scalefactor[0] <= 0.)
		    printf (" using image %s for the scalefactor reference value\n",image_file[image]);
		      // use the swarp astrometry, but decide whether to also correct the shear distortion or not
		  if (CORRECT_DISTORTION == 1) {
		    ngoodpix = tbswarpextract
		      (wcs_raw, apix, badpix, dim, ppix[chipnumber][nstar[chipnumber]],
		       dbadpix[chipnumber][nstar[chipnumber]], objpix, objx[i], objy[i],
		       poserror, pwidth, pheight, ikernel,
		       noise[chipnumber], region, maxlevel, scalefactor,
		       rawpos, wcspos, wcsneg, wcscentre,
		       kern, xfit, yfit, zfit, wfit, u, v);
		    // printf(" thread %d image %d object %d ngoodpix %d \n",nt,i,ii,ngoodpix); fflush(stdout);
		  } else {
		    if ((int)GREAT3 < 1)
		      flagbad = 1; // specify that background objects are to be flagged (for real data)
		    else
		      flagbad = 0; // specify that background objects are not flagged (for GREAT3)
		    ngoodpix = tbswarpextract_nodistortion
		      (wcs_raw, apix, badpix, opix, dim, ppix[chipnumber][nstar[chipnumber]],
		       dbadpix[chipnumber][nstar[chipnumber]], objx[i], objy[i],
		       poserror, pwidth, pheight,
		       noise[chipnumber], region, scalefactor, maxlevel,
		       rawpos, wcspos, wcsneg, wcscentre, flagbad);
		    // printf(" thread %d image %d object %d ngoodpix %d \n",nt,i,ii,ngoodpix); fflush(stdout);
		  }
		}

		// printf(" %f %f %lf %d \n",objx[i],objy[i],dpix[chipnumber][nstar[chipnumber]][0],ngoodpix);
		if (ngoodpix > 0) {
		  psfnorm = noisevalue = maxvalue = distantmaxvalue = 0.;
		  for(iy = 0; iy < pheight; iy++) {
		    yy = iy > pheight / 2 ? iy - pheight : iy;
		    for(ix = 0; ix < pwidth; ix++) {
		      xx = ix >
			pwidth / 2 ? ix - pwidth : ix;
		      if (xx * xx + yy * yy <= poserror*poserror) 
			{
			  pixel = ix + iy * pwidth;
			  psfnorm += ppix[chipnumber][nstar[chipnumber]][pixel];
			  noisevalue += noise[chipnumber] * noise[chipnumber];
			  if (ppix[chipnumber][nstar[chipnumber]][pixel] > maxvalue)
			    maxvalue = ppix[chipnumber][nstar[chipnumber]][pixel];
			}
		      else
			{
			  // check maximum pixel outside star area
			  pixel = ix + iy * pwidth;
			  if (ppix[chipnumber][nstar[chipnumber]][pixel] > distantmaxvalue)
			    distantmaxvalue = ppix[chipnumber][nstar[chipnumber]][pixel];
			}
		    }
		  }
		  if (noisevalue <= 0.) {
		    fflush(stdout);
		    fprintf(stderr,
			    " negative noise value \n");
		    exit(EXIT_FAILURE);
		  }
		  noisevalue = sqrt(noisevalue);
		  /*
		   * printf(" maxvalues %lf %lf %f %f %d %f %f \n",
		   * rdegs[i],ddegs[i],objx[i],objy[i],chipnumber+1,maxvalue/psfnorm,psfnorm/noisevalue);
		   */
		  if (psfnorm > (double) (snratio * noisevalue)) {
		    // check star isn't saturated nor a likely cosmic ray hit
		    //  and also that its maximum is significantly larger than any value outside poserror
		    if (maxvalue <= satlimit
			&& maxvalue < cosmicraylimit * psfnorm
			&& maxvalue > 2.*distantmaxvalue) {

		      // calculate star S/N ratio and weight for each star (downweight high S/N in chi-squared
		      // calculation so that the brightest stars above S/N=snref are all given equal weight) 
		      sn[chipnumber][nstar[chipnumber]] = psfnorm / noisevalue;
		      starweight[chipnumber][nstar[chipnumber]] = 1.;

		      // define global x,y position in range -1 < val < 1 
		      xval[chipnumber][nstar[chipnumber]] = (objx[i] + xchip * xchipsampling) / hxsize - 1.;
		      yval[chipnumber][nstar[chipnumber]] = objy[i] + ychip * ychipsampling;
		      chipxval[chipnumber][nstar[chipnumber]] = objx[i];
		      chipyval[chipnumber][nstar[chipnumber]] = objy[i];
                      chiprdegs[chipnumber][nstar[chipnumber]] = rdegs[i];
                      chipddegs[chipnumber][nstar[chipnumber]] = ddegs[i];
		      if (ychip > 0)
                        yval[chipnumber][nstar[chipnumber]] += big_gap;
		      if (ychip > 2)
                        yval[chipnumber][nstar[chipnumber]] += big_gap;
		      yval[chipnumber][nstar[chipnumber]] = yval[chipnumber][nstar[chipnumber]]/ hysize - 1.;
		      
		      nstar[chipnumber]++;
		    } else {
		      num_saturated++;
		    }
		  } else {
		    num_snratio++;
		  }
		} else {
		  num_badpix++;
		}
	      } else {
		num_wrong_mag++;
	      }
	    } else {
	      num_off++;
	    }
	  }
	}

	printf(" excluded stars: \n");
	printf(" %d stars outside image area \n", num_off);
	printf(" %d stars outside mag range \n", num_wrong_mag);
	printf(" %d stars with no good data \n", num_badpix);
	if (snratio > 0.)
	  printf(" %d stars with peak S/N ratio < %8.1f \n",
		 num_snratio, snratio);
	printf(" %d stars with peak > %8.0f ADU \n", num_saturated,
	       satlimit);
	printf(" %d stars in a close pair\n\n", num_close);
	printf(" %d stars matching selection criteria \n", nstar[chipnumber]);

        t3 = time(NULL);
  
	tot_nstar += nstar[chipnumber];

        free(region);
        free(apix);
        free(badpix);
	if (FILTER_WEIGHTS == 1)
	  {
	    free(opix);
	    free(changed);
	    free(weightfilterarray);
	  }
	if (SUBTRACT_MEDIAN == 1) free(sortarray);

	// free unused postage stamp memory
	if (nstar[chipnumber] < nobj) 
	  {
	    for(i = nstar[chipnumber]; i < nobj; i++)
	      {
		//printf(" free %d \n",i); fflush(stdout);
		free(ppix[chipnumber][i]);
		free(dbadpix[chipnumber][i]);
	      }
	  }

	// free swarp field structure
        free_cat(&rawcat, 1);
        free(rawfield);

        //printf(" memory freed ready for next image\n");

        memloop = 0.;

	// next image 
    }
    
    if (tot_nstar <= 0) {
      fprintf(stderr," no useful stars found in mosaic, no PSF created \n");
      exit(EXIT_FAILURE);
    }

    // remove the lock file
    if (lockfilename != NULL) 
      {
	remove(lockfilename);
      }

    /*
    for (image=0; image<nchip; image++)
      {
	for (i=0; i<nstar[image]; i++)
	  printf ("%d %d %g %g \n",image+1,i+1,ppix[image][i][0],dbadpix[image][i][0]);
      }
    */

    // set size of oversampled model PSF array
    ofactor = oversampling;
    printf(" model PSF oversampling factor %d \n",ofactor);
    if (ofactor<1 || ofactor%2 == 0)
      {
	fprintf(stderr,"invalid oversampling factor, must be odd integer\n");
	exit(EXIT_FAILURE);
      }
    fwidth = pwidth*ofactor;
    fheight = pheight*ofactor;
    fsize = fwidth*fheight;
    // allocate local PSF model memory
    psf = ml_calloc(fsize, sizeof(double), &memloop, "local model PSF");
    dpsf = ml_calloc(pwidth*pheight, sizeof(double), &memloop, "downsampled model PSF");
    respsf = ml_calloc(pwidth*pheight, sizeof(double), &memloop, "PSF residuals");
    fracrespsf = ml_calloc(pwidth*pheight, sizeof(double), &memloop, "Fractional PSF residuals");
    refpsf = ml_calloc(pwidth*pheight, sizeof(double), &memloop, "PSF at reference position");
    chisq = ml_calloc(fsize, sizeof(double), &memloop, "chisq array");
    // allocate memory for lookup table of selected pixels to be fitted
    pindex = ml_calloc(fsize, sizeof(int), &memloop, "forward pixel lookup table");
    oldpindex = ml_calloc(fsize, sizeof(int), &memloop, "forward pixel lookup table");

    // initialise model PSF and set number of pixels to be fitted, nfit
    // where the selected pixels are defined by the pindex array
    soln = ml_calloc(fsize*ncoeffs, sizeof(double), &memloop, "solution vector");    
    currentsoln = ml_calloc(fsize*ncoeffs, sizeof(double), &memloop, "solution vector");    
    oldsoln = ml_calloc(fsize*ncoeffs, sizeof(double), &memloop, "solution vector");    
    basesoln = ml_calloc(fsize*ncoeffs, sizeof(double), &memloop, "solution vector");    

    // call initialisePSF with nominal oversampling to fix array sizes correctly
    ofactor = oversampling;
    nfitmax = nfit = initialisePSF();
    printf(" %d oversampled model PSF pixels to be fitted \n",nfit);
    if (nfit>fsize)
      {
	fprintf(stderr," error in number of fittable pixels returned by initialisePSF\n");
	fflush(stderr);
	exit(EXIT_FAILURE);
      }

    // allocate mixing matrix and matrix arrays
    npix = ml_calloc((pwidth*pheight), sizeof(int), &memloop, "mixing matrix number of elements");    
    zpix = ml_calloc((pwidth*pheight), sizeof(int), &memloop, "zero pixel number of elements");    
    spix = ml_calloc((pwidth*pheight), sizeof(int*), &memloop, "mixing matrix image element pointer");    
    sval = ml_calloc((pwidth*pheight), sizeof(double*), &memloop, "mixing matrix image element value pointer");    
    zval = ml_calloc((pwidth*pheight), sizeof(double*), &memloop, "zero pixel image element value pointer");    
    for (i=0; i<pwidth*pheight; i++)
      {
	spix[i] = ml_calloc(nfit, sizeof(int), &memloop, "mixing matrix nfit elements");    
	sval[i] = ml_calloc(nfit, sizeof(double), &memloop, "mixing matrix image element values");    
	zval[i] = ml_calloc(nfit, sizeof(double), &memloop, "zero pixel image element values");    
      }
    dsize = nfit*ncoeffs;

    rchip = ml_calloc(nchip, sizeof(int), &memloop, "random chip numbers");

    // array for star shift values
    shift = ml_calloc(nchip, sizeof(double**), &memloop, "shift pointers");
    bestshift = ml_calloc(nchip, sizeof(double**), &memloop, "shift pointers");
    initshift = ml_calloc(nchip, sizeof(double**), &memloop, "shift pointers");
    for (k=0; k<nchip; k++)
      {
	if (nstar[k]>0)
	  {
	    shift[k] = ml_calloc(nstar[k], sizeof(double*), &memloop, "shift[k] pointers");
	    bestshift[k] = ml_calloc(nstar[k], sizeof(double*), &memloop, "shift[k] pointers");
	    initshift[k] = ml_calloc(nstar[k], sizeof(double*), &memloop, "shift[k] pointers");
	    for (i=0; i<nstar[k]; i++)
	      {
		shift[k][i] = ml_calloc(2, sizeof(double), &memloop, "shift values");
		bestshift[k][i] = ml_calloc(2, sizeof(double), &memloop, "shift values");
		initshift[k][i] = ml_calloc(2, sizeof(double), &memloop, "shift values");
	      }
	  }
      }

    printf(" total stars %d \n",tot_nstar);
    printf(" approx total fitted pixels %d \n",tot_nstar*nfit/ofactor/ofactor);
    printf(" number of free parameters %d \n",nfit*ncoeffs);

    // start the iterations with no oversampling to prevent solution becoming
    // stuck in a false minimum.  Call initialisePSF with no oversampling and
    // no chip variation allowed, to 
    // set the initial solution vector
    chipvariation = 0;
    ncoeffs = (1+order)*(2+order)/2; 
    for (i=0; i<fsize; i++) pindex[i] = 0;
    for (i=0; i<fsize*ncoeffs; i++) soln[i] = basesoln[i] = 0.;
    ofactor = 1;
    fwidth = pwidth*ofactor;
    fheight = pheight*ofactor;
    fsize = fwidth*fheight;
    nfit = initialisePSF();
    if (nfit>fsize)
      {
	fprintf(stderr," error in number of fittable pixels returned by initialisePSF %d\n",nfit);
	fflush(stderr);
	exit(EXIT_FAILURE);
      }
    dsize = nfit*ncoeffs;
    //outputfstamp("initialpsf.fits");

    // copy data pixel values into array for safe keeping
    bestchisq = nullchisq = chisqmin = storepixelvalues();
    printf(" null model chi-squared = %g \n", nullchisq);

    // measure the initial position shifts with this starting model
    // cross-correlate current model with stars, correct for non-zero model centroid 
    // and fill array of shift values to be used when reconstructing the PSF for each star
    range = 1;
    lchisq = 0.;
    tempshift[0]=tempshift[1]=cent[0]=cent[1]=0.;
    printf("chip star  x  y nset used shifts  cross-corr shifts  centroids  amp chisq\n");
    for (k=0; k<nchip; k++)
      {
	if (nstar[k]>0)
	  {
	    for (i=0; i<nstar[k]; i++)
	      {
		printout=0;
		// print out some detail of the fitting for the shift corrections
		// for one star
		// if (k==14 && i==40) printout=1;
		// reconstruct oversampled, centred PSF from coefficients
		reconstructosPSF(xval[k][i], yval[k][i], k); // printf("recon\n"); fflush(stdout);
		// cross-correlate
		crosscorrelate(range, chisq, dpix[k][i], k, initshift[k][i], &bestamp, &lchisq0); // printf("cross\n"); fflush(stdout);
		// generate mixing matrix for this position shift
		nset = smatrix(initshift[k][i]); // printf("smatrix\n"); fflush(stdout);
		downsample(); // printf("ds\n"); fflush(stdout);
		// get amplitude, chi-squared for this star (actually here this should yield same values as returns from crosscorrelate)
		lchisq0 = chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum); // printf("chisq\n"); fflush(stdout);
		lchisq += starweight[k][i]*lchisq0;
		if (VERBOSE == 1)
		  printf("initial shifts %d %d %f %f %d %f %f %f %f %f %f %f %f \n",
			 k+1,i+1,xval[k][i],yval[k][i],nset,initshift[k][i][0],initshift[k][i][1],
			 tempshift[0],tempshift[1],cent[0],cent[1],amp,lchisq0); fflush(stdout);
		shift[k][i][0] = initshift[k][i][0];
		shift[k][i][1] = initshift[k][i][1];
		offset = sqrt(pow(shift[k][i][0],2)+pow(shift[k][i][1],2));
		if (offset > poserror)
		  {
		    shift[k][i][0]=shift[k][i][1]=initshift[k][i][0]=initshift[k][i][1]=0.;
		    starweight[k][i]=0.;
		  }
	      }
	  }
      }
    printf(" initial chi-squared %f\n",lchisq);

    // set which variables will be included in the iteration loops
    oversampling_loop = 1;  // allow oversampling 
    chipvariation_loop = 2; // allow chip variation 
    randomise_loop = 3; // allow randomisation to improve global minimisation
    bval_loop = 18;  // allow bval testing after loop 18
    os_niter = 4;  // total number of passes
    final_nitermin = 1; // minimum number of passes in the final loop (to try to guarantee convergence)
    targetchiporder = chiporder; // final chiporder required
    newchiporder = chiporder = -1; // initial chiporder to be used
    copied = 0;

    printf(" %d passes will be made \n",os_niter);

    t4 = time(NULL);

    // make multiple passes, the first with unit oversampling and no chip variation, the second with
    // final oversampling and then chip variation as specified
    for (os_iter=0; os_iter<os_niter; os_iter++)
      {

	// enable oversampling if specified 
	if (os_iter==oversampling_loop && oversampling>1)
	  {
	    ofactor = oversampling;
	    fwidth = pwidth*ofactor;
	    fheight = pheight*ofactor;
	    fsize = fwidth*fheight;
	    nfit = reinitialisePSF();
	    if (nfit>fsize)
	      {
		fprintf(stderr," error in number of fittable pixels returned by reinitialisePSF %d\n",nfit);
		fflush(stderr);
		exit(EXIT_FAILURE);
	      }
	    dsize = nfit*ncoeffs;
	    if (oversampling_loop > chipvariation_loop || newchiporder >= targetchiporder)
	      {
		for (i=0; i<dsize; i++)
		  {
		    basesoln[i] = soln[i];
		    soln[i] = currentsoln[i] = 0.;
		  }
	      }
	    printf(" oversampling %d enabled\n",oversampling); fflush(stdout);
	    // test the new value of chi-squared (should be the same as at the end of the previous iteration
	  }

	// enable chip variation 
	if (os_iter == chipvariation_loop  && newchiporder < targetchiporder)
	  {
	    newchiporder = targetchiporder; // jump to desired chip order in one leap
	    enablechipvariation(newchiporder);
	    printf(" discontinuous variations between chips enabled: chip order = %d\n",chiporder);
	    if (chipvariation_loop >= oversampling_loop || oversampling <= 1)
	      {
		for (i=0; i<dsize; i++)
		  {
		    basesoln[i] = soln[i];
		    soln[i] = currentsoln[i] = 0.;
		  }
	      }
	  }

	if (os_iter >= randomise_loop && os_iter <= bval_loop) 
	  {
	    // replace the soln vector to be randomised 
	    if (copied>=1)
	      {
		for (i=0; i<dsize; i++)
		  {
		    soln[i] = currentsoln[i];
		  }
		for (k=0; k<nchip; k++)
		  {
		    if (nstar[k]>0)
		      {
			for (i=0; i<nstar[k]; i++)
			  {
			    shift[k][i][0] = bestshift[k][i][0];
			    shift[k][i][1] = bestshift[k][i][1];
			  }
		      }
		  }
	      }
	    // loop to randomse the chip variations between chips in an attempt to step out of
	    // local minima and obtain a better global minimum
	  }
	if (os_iter >= randomise_loop && os_iter < bval_loop) 
	  {
	    randomisesoln();
	    printf(" solution vector randomised\n");
	  }

	// if allowing non-linearity to be fitted
	if (os_iter >= bval_loop)
	  {
	    bval += 0.2e-6;
	    // calculate the non-linearity coefficient and correct the data
	    chisqmin = storepixelvalues();
	    // chisqmin = nonlinearity();
	    deltachisq = pchisq-chisqmin;
	    printf(" bval = %g initial chisq, deltachisq = %f %f \n", bval, chisqmin, deltachisq);
	  }

	// calculate chi-squared immediately after initialising the solution vector
	lchisq = 0.;
	for (k=0; k<nchip; k++)
	  {
	    for (i=0; i<nstar[k]; i++)
	      {
		// create shifted, downsampled model
		reconstructosPSF(xval[k][i], yval[k][i], k);
		nset = smatrix(shift[k][i]);
		downsample();
		// get amplitude, chi-squared for this star
		lchisq += starweight[k][i]*chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum);
	      }
	  }
	printf(" chi-squared on entry = %f\n",lchisq); fflush(stdout);

	// set limit for halting iterations
	if (os_iter==0)
	  deltachisqlimit = (nfit*ncoeffs)/100.;
	else
	  deltachisqlimit = (nfit*ncoeffs)/1000.;

	printf(" outer loop iteration %d delta-chi-squared limit %f \n",os_iter+1, deltachisqlimit);

	deltachisq = deltachisqlimit + 10.;
	previousdelta = deltachisq;
	iter = 0;
	// keep iterating until either two successive delta-chisq values are below threshold or number limit reached
	while (iter<niter && 
	       ( deltachisq > deltachisqlimit 
		 || previousdelta > deltachisqlimit 
		 || (os_iter==os_niter-1 && iter<final_nitermin)) 
	       )
	  {
	    previousdelta = deltachisq;
	    pchisq = chisqmin;
	    printf(" position iteration %d \n",iter+1); fflush(stdout);

	    // gradient solve
	    chisqmin = solvePSF();
	    deltachisq = pchisq-chisqmin;
	    printf (" solvePSF chi-squared = %f delta chi squared = %f \n",chisqmin, deltachisq);
	    fflush(stdout);

	    // work through each star and cross-correlate the model with the star
	    // to determine its best fit shift - only after a first minimisation with
	    // the current set of positions has been obtained
	    printf(" measuring/correcting positions\n"); fflush(stdout);

	    // cross-correlate current model with stars, correct for non-zero model centroid 
	    // and fill array of shift values to be used when reconstructing the PSF for each star
	    lchisq = 0.;
	    if (VERBOSE == 1)
	      printf("os_iter iter chip star   x   y  nset  used   shifts  cross-corr shifts  centroids  amp chisq\n");
	    for (k=0; k<nchip; k++)
	      {
		if (nstar[k]>0)
		  {
		    for (i=0; i<nstar[k]; i++)
		      {
			current_star = i+1;
			printout=0;
			// print out some detail of the fitting for the shift corrections
			// for one star
			// if (k==14 && i==40) printout=1;
			// reconstruct oversampled, centred PSF from coefficients
			reconstructosPSF(xval[k][i], yval[k][i], k); // printf("recon\n"); fflush(stdout);
			// find model centroid
			// here, we want to keep the PSF peak centred, so the weight function should be narrow
			findmodelcentroid(cent, sigma); // printf("centroid\n"); fflush(stdout);
			// correct current star offset by model centroid offset in attempt to keep models centred.
			tempshift[0] = shift[k][i][0] + cent[0];
			tempshift[1] = shift[k][i][1] + cent[1];
			// cross-correlate using the current shift as an input starting point
			crosscorrelate(range, chisq, dpix[k][i], k, tempshift, &bestamp, &lchisq0); 
			// check that measured shift has not gone out of range: if it has, reset it to the initial value
			if (tempshift[0]<=-pwidth/2 || tempshift[0]>=pwidth/2 || 
			    tempshift[1]<=-pheight/2 || tempshift[1]>=pheight/2)
			  {
			    shift[k][i][0] = initshift[k][i][0];
			    shift[k][i][1] = initshift[k][i][1];
			    printf(" WARNING: shift out of range for star %d %d, %f %f\n",k+1,i+1,tempshift[0],tempshift[1]);
			  }
			else
			  {
			    // correct star offset by model centroid offset in attempt to keep models centred.
			    shift[k][i][0] = tempshift[0] - cent[0];
			    shift[k][i][1] = tempshift[1] - cent[1];
			  }
			// generate mixing matrix for this position shift
			nset = smatrix(shift[k][i]); // printf("smatrix\n"); fflush(stdout);
			downsample(); // printf("ds\n"); fflush(stdout);
			// get amplitude, chi-squared for this star
			lchisq0 = chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum); // printf("chisq\n"); fflush(stdout);
			chisqval[k][i] = lchisq0;
			lchisq += starweight[k][i]*lchisq0;
			if (VERBOSE == 1)
			  printf("shifts %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f \n",os_iter,iter,
				 k+1,i+1,xval[k][i],yval[k][i],shift[k][i][0],shift[k][i][1],
				 tempshift[0],tempshift[1],cent[0],cent[1],starweight[k][i],amp,lchisq0); fflush(stdout);
			current_star = 0;
		      }
		  }
	      }
	    printf(" chi-squared after star shifts %f delta chi-sq %f \n",lchisq,pchisq-lchisq); fflush(stdout);

	    // update minimum chi-squared valaue and delta-chi-squared after position shifts applied
	    chisqmin = lchisq;
	    deltachisq = pchisq-chisqmin;

	    // keep best-fit solution and shifts
	    if (os_iter==chipvariation_loop || os_iter>=randomise_loop)
	      {
		if (chisqmin < bestchisq)
		  {
		    printf(" reset best chi-squared from %f to %f at outer loop %d\n",bestchisq,chisqmin,os_iter+1); fflush(stdout);
		    copied = 1;
		    bestchisq = chisqmin;
		    for (i=0; i<dsize; i++) 
		      {
			currentsoln[i] = soln[i];
		      }
		    for (k=0; k<nchip; k++)
		      {
			if (nstar[k]>0)
			  {
			    for (i=0; i<nstar[k]; i++)
			      {
				bestshift[k][i][0] = shift[k][i][0];
				bestshift[k][i][1] = shift[k][i][1];
			      }
			  }
		      }
		  }
	      }
		       
	    // next iteration
	    iter++;
	  }

	// if we're into the final stages, test the chi-squared values and set weight to zero
	// for any star with a big chi-squared value or with a big shift value
	if (os_iter >= chipvariation_loop)
	  {
	    chisq_expected = pwidth*pheight;  // approx expected value
	    for (k=0; k<nchip; k++)
	      {
		if (nstar[k]>0)
		  {
		    for (i=0; i<nstar[k]; i++)
		      {
			// find position offset 
			offset = sqrt(pow(shift[k][i][0],2)+pow(shift[k][i][1],2));
			// eliminate star if chisq is too large or position offset too large
			if (chisqval[k][i] > 4.*chisq_expected || offset > poserror)
			  {
			    if (starweight[k][i] > 0.) 
			      {
				printf(" downweighted star %d %d chisq %f\n",k+1, i+1, chisqval[k][i]);
				if (chisqmin <= bestchisq)
				  {
				    bestchisq -= chisqval[k][i];
				  }
				chisqmin -= chisqval[k][i];
			      }
			    starweight[k][i] = 0.;
			  }
			else if (starweight[k][i] < 1. && chisqval[k][i] < 3.*chisq_expected && offset < poserror)
			  {
			    // include this star again if chisq is now OK
			    starweight[k][i] = 1.;
			    printf(" upweighted star %d %d chisq %f\n",k+1, i+1, chisqval[k][i]);
			    if (chisqmin <= bestchisq)
			      {
				bestchisq += chisqval[k][i];
			      }
			    chisqmin += chisqval[k][i];
			  }
		      }
		  }
	      }
	    if (chisqmin < bestchisq) bestchisq = chisqmin;
	  }

	if (os_iter >= bval_loop)
	  {
	    printf(" bval = %g final chisq, deltachisq = %f %f \n", bval, chisqmin, deltachisq);
	  }

	// next outer loop
      }

    t5 = time(NULL);
    printf (" Iterations time was %g secs \n", difftime(t5, t4) );

    // iterations have finished
    // update basesoln and zero the other vectors
    for (i=0; i<dsize; i++) 
      {
	basesoln[i] += currentsoln[i];
	soln[i] = currentsoln[i] = 0.;
      }


    // count numbers of used stars
    for (k=0; k<nchip; k++)
      {
	numstar[k] = 0;
	if (nstar[k]>0)
	  {
	    for (i=0; i<nstar[k]; i++)
	      {
		if (starweight[k][i] > 0.) 
		  {
		    numstar[k]++;
		  }
	      }
	  }
      }

    // create and write out downsampled PSF model suitable for lensfit
    downsamplePSFmodel(coeffname, xchipsampling, ychipsampling, big_gap, 
		       nxchip, nychip, numstar, scalefactor,
		       hxsize, hysize);

    printf(" downsampled PSF model written to %s \n", coeffname); fflush(stdout);

    // write oversampled solution vector
    writeoscoeffs(ospsfname, correct_distortion, pindex, basesoln, nfit, ncoeffs, 
		  order, chipvariation, chiporder, scalefactor,
		  fwidth, fheight, oversampling, nxchip, nychip, numstar, xchipsampling, ychipsampling,
		  big_gap, hxsize, hysize, poserror);

    printf(" oversampled PSF model written to %s \n", ospsfname); fflush(stdout);



    /*************************************************************
     output diagnostic images, residuals etc
    *************************************************************/

    // create output FITS cubes of downsampled psf models and stars and write shifts to log file

    lchisq = 0.;
    for (k=0; k<nchip; k++)
      {

       // Write downsampled, centred PSF at reference positions
       remove(refpsfcubefile[k]);
       if (nrefpsf>0) 
	 {
	   open3Dfits(refpsfcubefile[k], &refpsfcubepointer, nrefpsf);
	   refpsfx_glob = ml_calloc(nrefpsf, sizeof(double), &memory, "refpsfx_glob");
	   refpsfy_glob = ml_calloc(nrefpsf, sizeof(double), &memory, "refpsfy_glob");
	 }
       ychip = k / nxchip;         // e.g. for chip k = 11 (which is really chip 12...), ychip = 2
       xchip = k - ychip * nxchip; // and xchip = 3
       for ( i=0; i<nrefpsf; i++ ) {
         // Convert the chip coordinates in refpsfx and refpsfy into global coordinates in range -1 < x,y < 1
         refpsfx_glob[i] = ( refpsfx[i] + xchip * xchipsampling ) / hxsize - 1.;
         refpsfy_glob[i] = refpsfy[i] + ychip * ychipsampling;
         if ( ychip > 0 || ychip > 2 ) 
             refpsfy_glob[i] += big_gap;
         refpsfy_glob[i] = refpsfy_glob[i] / hysize - 1.;
         // Reconstruct the PSF at each of those points for given chip number (this is stored in psf[])
         reconstructosPSF(refpsfx_glob[i], refpsfy_glob[i], k);
         // Make sure PSF is centred
         refpsfpos[0] = refpsfpos[1] = 0.;
         refpsfnset = smatrix(refpsfpos);
         // Take oversampled psf[] and downsample it (storing it in dpsf[])
         downsample();
         // Make sure order of array is correct for output
         for ( ref_px=0; ref_px<pwidth*pheight; ref_px++ ) {
           ref_y = ref_px / pwidth;
           ref_x = ref_px - ref_y*pwidth;
           ref_x += pwidth/2;
           ref_y += pwidth/2;
           if ( ref_x >= pwidth  ) ref_x -= pwidth;
           if ( ref_y >= pheight ) ref_y -= pheight;
           ref_i = ref_y*pwidth + ref_x;
           refpsf[ref_px] = dpsf[ref_px];
         }
	 if (VERBOSE == 1)
	   printf(" Chip %d (xchip=%d and ychip=%d): writing PSF at reference point %f %f (global coordinates %f %f) \n", k, xchip, ychip, refpsfx[i], refpsfy[i], refpsfx_glob[i], refpsfy_glob[i]);
         write3Dfits(refpsfcubepointer, i, dpsf);
       }
       if (nrefpsf>0) fits_close_file(refpsfcubepointer, &status);

	if (nstar[k]>0)
	  {
	    remove(starcubefile[k]);
	    open3Dfits(starcubefile[k], &starcubepointer, nstar[k]);
	    remove(psfcubefile[k]);
	    open3Dfits(psfcubefile[k], &psfcubepointer, nstar[k]);
	    remove(rescubefile[k]);
	    open3Dfits(rescubefile[k], &rescubepointer, nstar[k]);
            remove(fracrescubefile[k]);
            open3Dfits(fracrescubefile[k], &fracrescubepointer, nstar[k]);
	    if ((shiftsfile = fopen(shiftsname[k], "w")) == NULL) {
	      fprintf(stderr, " failed to make file %s \n", shiftsname[k]);
	    }
            if ((ellfile = fopen(ellname[k], "w")) == NULL) {
              fprintf(stderr, "failed to make file %s \n", ellname[k]);
            }
	    for (i=0; i<nstar[k]; i++)
	      {
		current_star = i+1;
		// create shifted, downsampled model
		reconstructosPSF(xval[k][i], yval[k][i], k);
		// measure shifts without centroid correction
		crosscorrelate(range, chisq, dpix[k][i], k, shift[k][i], &bestamp, &lchisq0); 
		// create shifted downsampled PSF
		nset = smatrix(shift[k][i]);
		downsample();
		// get amplitude, chi-squared for this star
		lchisq = chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum);
		// print out to log file
		//            fprintf(shiftsfile,"%d \n", iter);
		fprintf(shiftsfile,"%d %d %f %f %f %f %f %f %f %f %f \n",
			k+1, i+1, chiprdegs[k][i], chipddegs[k][i], starweight[k][i],
			chipxval[k][i],chipyval[k][i],shift[k][i][0],shift[k][i][1],amp,lchisq);
		// Calculate ellipticity of PSF-model at this point
                psfmoments(dpsf, pheight, pwidth, modelpsfe, modelcentroid, modelmeanmoments);

                // Calculate ellipticity of star
                psfmoments(dpix[k][i], pheight, pwidth, starpsfe, starcentroid, starmeanmoments);

                // Write ellipticities to file
                fprintf(ellfile,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", i+1, chiprdegs[k][i], chipddegs[k][i], 
			modelcentroid[0], modelcentroid[1], modelpsfe[0], modelpsfe[1], modelmeanmoments[0],modelmeanmoments[1],modelmeanmoments[2], 
			starcentroid[0], starcentroid[1], starpsfe[0], starpsfe[1], starmeanmoments[0],starmeanmoments[1],starmeanmoments[2], starweight[k][i]);

                // at this point the star image is in array dpix[k][i][] and the shifted, downsampled PSF is in dpsf[]
		// (i.e. dpsf is  1D array that gets overwritten each time downsample is called)
		// now write the stars and PSFs out to an output FITS cube
		write3Dfits(psfcubepointer, i, dpsf);
		// write out the original star array, without the nonlinearity correction
		write3Dfits(starcubepointer, i, ppix[k][i]);
		// Also output the PSF residuals (dpsf is normalised so multiply by amplitude)
		// But first find normalisation factors (normalise so sum of all pixels is 1).
		norm_model = 0.;
		norm_star  = 0.;
		for (l=0; l<pwidth*pheight; l++) {
		  norm_star  = norm_star  + dpix[k][i][l];
		  norm_model = norm_model + dpsf[l];
		}
		worstpixval = 0.;
		ressum      = 0.;
		for (l=0; l<pwidth*pheight; l++) {
		  //                respsf[l] = dpix[k][i][l]/norm_star - dpsf[l]/norm_model;
		  // normalise residuals by the pixel noise
		  nval = pow(noise[k],2);
		  if (dpix[k][i][l] > 0.) nval += dpix[k][i][l]/gain[k]; // shot noise squared in ADU 
		  nval = sqrt(nval);
		  respsf[l]     = (dpix[k][i][l] - amp*dpsf[l]) / nval;
		  if (fabs(amp*dpsf[l]) > fabs(dpix[k][i][l] - amp*dpsf[l]))
		    fracrespsf[l] = (dpix[k][i][l] - amp*dpsf[l]) / (amp*dpsf[l]);
		  else
		    fracrespsf[l] = 0.;
		  // record the sum of residuals, and the worst pixel
		  ressum = ressum + respsf[l];
		  if(pow(respsf[l],2)>pow(worstpixval,2)) {
		    worstpix      = l;
		    worstpixval   = respsf[l];
		    worstpixratio = nval*respsf[l]/dpix[k][i][l];
		  }
		  /*
		    if(l<10) {
                    fprintf(shiftsfile,"%d %f %f %f %f %f %f \n", l, dpix[k][i][l], dpsf[l], amp, respsf[l], norm_star, norm_model);
		    }
		  */
		}
                // Write residuals (star - psfmodel) to FITS cube
		write3Dfits(rescubepointer, i, respsf);
                // Write fractional residuals ((star - psfmodel) / psfmodel) to FITS cube
                write3Dfits(fracrescubepointer, i, fracrespsf);
		fprintf(shiftsfile, " Residuals: %d %f %f %d %f %f %lf %lf \n", i+1, starweight[k][i], ressum, worstpix, worstpixval, worstpixratio, chiprdegs[k][i], chipddegs[k][i]);
	      }
	    fprintf(shiftsfile, " Iterations: %d \n", iter);
	    fits_close_file(starcubepointer, &status);
	    fits_close_file(psfcubepointer, &status);
	    fits_close_file(rescubepointer, &status);
            fits_close_file(fracrescubepointer, &status);
	    fclose(shiftsfile);
	  }
      }




    /*****************************************************************************************
      create table of fitted shift values (replacement for globalshifts program, the following
      code is copied from that)
    ******************************************************************************************/

    // recalculate star shifts 
    if (VERBOSE == 1)
      printf("os_iter iter chip star   x   y     shifts   \n");
    nn = 0;
    for (k=0; k<nchip; k++)
      {
	if (nstar[k]>0)
	  {
	    for (i=0; i<nstar[k]; i++)
	      {
		printout=0;
		// reconstruct oversampled, centred PSF from coefficients
		reconstructosPSF(xval[k][i], yval[k][i], k); 
		// find oversampled model centroid with same weight function as wehn creating models
		findmodelcentroid(modelcentroid, sigma); 
		// find oversampled model centroid with broader weight function for astrometry
		// which should roughly match the Sextractor weight function used for the input catalogue
		findmodelcentroid(astrometrycentroid, wsigma); 
		pos[0] = pos[1] = 0.;
		// generate mixing matrix for zero position shift
		nset = smatrix(pos); 
		downsample(); // downsampled PSF model is now in array dpsf
		// find the centroid position from the downsampled model using the same routine as assumed in lensfit
		psfmoments(dpsf, pheight, pwidth, psfe, cent, moments);
		// evaluate the extraction pixelisation shift
		xs = chipxval[k][i] - (int)(chipxval[k][i] + 0.5);
		ys = chipyval[k][i] - (int)(chipyval[k][i] + 0.5);
		// remove original model centroid correction
		bestshift[k][i][0] += modelcentroid[0];
		bestshift[k][i][1] += modelcentroid[1];
		// generate shifts with low-resolution, broadly weighted centroid correction applied
		// flip the sign to match the lensfit convention!
		shift[k][i][0] = -bestshift[k][i][0] + cent[0] - xs;
		shift[k][i][1] = -bestshift[k][i][1] + cent[1] - ys;
		// generate shifts with the full-resolution model astrometry centroid correction applied
		// flip the sign to match the lensfit convention!
		bestshift[k][i][0] = -bestshift[k][i][0] + astrometrycentroid[0] - xs;
		bestshift[k][i][1] = -bestshift[k][i][1] + astrometrycentroid[1] - ys;
		if (VERBOSE == 1)
		  {
		    printf(" %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f \n",
			   k,i,chiprdegs[k][i],chipddegs[k][i],starweight[k][i],shift[k][i][0],shift[k][i][1],
			   bestshift[k][i][0],bestshift[k][i][1],
			   cent[0],cent[1],modelcentroid[0],modelcentroid[1],xs,ys);
		  }
		// accumulate values for astrometric residuals fitting
		if (starweight[k][i] > 0.) 
		  {
		    chips[nn] = k;
		    // modified S/N weighting not to give too much weight to individual bright stars
		    wfits[nn] = starweight[k][i]/(1. + 50.*50./(sn[k][i]*sn[k][i]));
		    xfits[nn] = xval[k][i];
		    yfits[nn] = yval[k][i];
		    xshiftfits[nn] = shift[k][i][0];
		    yshiftfits[nn] = shift[k][i][1];
		    bxshiftfits[nn] = bestshift[k][i][0];
		    byshiftfits[nn] = bestshift[k][i][1];
		    rx[nn] = chiprdegs[k][i];
		    ry[nn] = chipddegs[k][i];
		    //printf(" %d %d %d %f %d %lf %lf \n",nn,k,i,wfits[nn],chips[nn],chiprdegs[k][i],chipddegs[k][i]);
		    // count total number of useful values 
		    nn++;
		  }
	      }
	  }
      }

    if (VERBOSE == 1)
      printf(" shift arrays filled with %d values \n", nn); fflush(stdout);


    if (FIT_ASTROMETRIC_ERRORS == 1)
      {

    // store original weights 
    for(i = 0; i < nn; i++)
        ow[i] = wfits[i];

    xrms = yrms = 1.;

    for(iter = 0; iter < 3; iter++) {
        printf(" star shift fitting, iteration %d \n", iter + 1);
        /*
         * fit polynomials to x & y shifts across mosaic (use same routine as for PSF pixel values 
         */
        globalsvdfit(xfits, yfits, xshiftfits, wfits, chips, nn,
                     nchip, sorder, scrossterm, schipvariation, schiporder,
                     sncoeffs, avals, u, v, w);
        // put fit coefficients into global pixel array 
        for(j = 0; j < sncoeffs; j++) {
            xcoeffs[j] = avals[j + 1];  /* shift by 1 from NR routine */
        }
        globalsvdfit(xfits, yfits, yshiftfits, wfits, chips, nn,
                     nchip, sorder, scrossterm, schipvariation, schiporder,
                     sncoeffs, avals, u, v, w);
        // put fit coefficients into global pixel array 
        for(j = 0; j < sncoeffs; j++) {
            ycoeffs[j] = avals[j + 1];  /* shift by 1 from NR routine */
        }

        sum = sumx = sumy = sumxsq = sumysq = 0.;

        /*
         * test if any stars look way off, set their weights to zero and refit 
         */
        for(i = 0; i < nn; i++) {
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              xcoeffs, &xs);
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              ycoeffs, &ys);
            rdx = xs - xshiftfits[i];
            rdy = ys - yshiftfits[i];
            /*
             * reject any star that's out by more than 0.3 pixels or 3 sigma
             */
            if (fabs(rdx) > 0.3 || fabs(rdx) > 3. * xrms || fabs(rdy) > 0.3
                || fabs(rdy) > 3. * yrms) {
                if (wfits[i] > 0.) {
                    printf(" rejected star %lf %lf %d %f %f %f %f %f \n",
                           rx[i], ry[i], chips[i], wfits[i], xshiftfits[i], xs,
                           yshiftfits[i], ys);
                }
                wfits[i] = 0.;
            } else {
                /*
                 * restore weights back to original values 
                 */
                wfits[i] = ow[i];
                sum++;
                sumx += rdx;
                sumy += rdy;
                sumxsq += rdx * rdx;
                sumysq += rdy * rdy;
                // if (iter==2) printf (" good star %f %f %f %f \n",xshiftfits[i],xs,yshiftfits[i],ys);
            }
        }
        if (sum > 3) {
            xrms = sqrt((sumxsq - sumx * sumx / sum) / (sum - 1.));
            yrms = sqrt((sumysq - sumy * sumy / sum) / (sum - 1.));
        } else {
            xrms = 1.;
            yrms = 1.;
        }
    }

    printf
        ("\n good stars ra, dec, chip, weight, measured xshift, fitted xshift, measured yshift, fitted yshift\n");
    fflush(stdout);

    for(i = 0; i < nn; i++) {
        if (wfits[i] > 0.) {
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              xcoeffs, &xs);
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              ycoeffs, &ys);
            printf(" %lf %lf %d %f %f %f %f \n", rx[i], ry[i], chips[i],
                   xshiftfits[i], xs, yshiftfits[i], ys);
        }
    }

      }
    else
      {
	printf(" not fitting astrometric errors \n");
      }

    updateglobalcoeffs(coeffname, sncoeffs, sorder, schipvariation,
                       schiporder, xcoeffs, ycoeffs);

    printf(" coefficients FITS file %s updated \n", coeffname);
    fflush(stdout);

    /*****************************************************************************************
     repeat for the full-resolution shifts
    ******************************************************************************************/

    if (FIT_ASTROMETRIC_ERRORS == 1)
      {

    // restore original weights 
    for(i = 0; i < nn; i++)
        wfits[i] = ow[i];

    xrms = yrms = 1.;

    for(iter = 0; iter < 3; iter++) {
        printf(" full-resolution star shift fitting, iteration %d \n", iter + 1);
	fflush(stdout);
        /*
         * fit polynomials to x & y shifts across mosaic (use same routine as for PSF pixel values 
         */
        globalsvdfit(xfits, yfits, bxshiftfits, wfits, chips, nn,
                     nchip, sorder, scrossterm, schipvariation, schiporder,
                     sncoeffs, avals, u, v, w);
        // put fit coefficients into global pixel array 
        for(j = 0; j < sncoeffs; j++) {
            xcoeffs[j] = avals[j + 1];  /* shift by 1 from NR routine */
        }
        globalsvdfit(xfits, yfits, byshiftfits, wfits, chips, nn,
                     nchip, sorder, scrossterm, schipvariation, schiporder,
                     sncoeffs, avals, u, v, w);
        // put fit coefficients into global pixel array 
        for(j = 0; j < sncoeffs; j++) {
            ycoeffs[j] = avals[j + 1];  /* shift by 1 from NR routine */
        }

        sum = sumx = sumy = sumxsq = sumysq = 0.;

        /*
         * test if any stars look way off, set their weights to zero and refit 
         */
        for(i = 0; i < nn; i++) {
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              xcoeffs, &xs);
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              ycoeffs, &ys);
            rdx = xs - bxshiftfits[i];
            rdy = ys - byshiftfits[i];
            /*
             * reject any star that's out by more than 0.3 pixels or 3 sigma
             */
            if (fabs(rdx) > 0.3 || fabs(rdx) > 3. * xrms || fabs(rdy) > 0.3
                || fabs(rdy) > 3. * yrms) {
                if (wfits[i] > 0.) {
                    printf(" rejected star %lf %lf %d %f %f %f %f %f \n",
                           rx[i], ry[i], chips[i], wfits[i], bxshiftfits[i], xs,
                           byshiftfits[i], ys);
                }
                wfits[i] = 0.;
            } else {
                /*
                 * restore weights back to original values 
                 */
                wfits[i] = ow[i];
                sum++;
                sumx += rdx;
                sumy += rdy;
                sumxsq += rdx * rdx;
                sumysq += rdy * rdy;
                // if (iter==2) printf (" good star %f %f %f %f \n",bxshiftfits[i],xs,byshiftfits[i],ys);
            }
        }
        if (sum > 3) {
            xrms = sqrt((sumxsq - sumx * sumx / sum) / (sum - 1.));
            yrms = sqrt((sumysq - sumy * sumy / sum) / (sum - 1.));
        } else {
            xrms = 1.;
            yrms = 1.;
        }
    }

    printf
        ("\n good stars ra, dec, chip, weight, measured xshift, fitted xshift, measured yshift, fitted yshift\n");
    fflush(stdout);
    for(i = 0; i < nn; i++) {
        if (wfits[i] > 0.) {
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              xcoeffs, &xs);
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chips[i], nchip,
                              ycoeffs, &ys);
	    printf(" %lf %lf %d %f %f %f %f \n", rx[i], ry[i], chips[i],
		   bxshiftfits[i], xs, byshiftfits[i], ys);
        }
    }

      }

    updateoscoeffs(ospsfname, sncoeffs, sorder, schipvariation,
		   schiporder, wsigma, xcoeffs, ycoeffs);

    printf(" o/s coefficients FITS file %s updated \n", ospsfname);
    fflush(stdout);


    printf(" Memory in use at end of program = %f MB\n", memory/pow(1024,2));
    fflush(stdout);

    return 0;
}  



double solvePSF()
{
  // GSL gradient method

       size_t iter = 0;
       int i, status, iterationlimit;
       double lchisq, pchisq, deltachisq;
       double stepsize, tol;
       char *fname;
     
       const gsl_multimin_fdfminimizer_type *T;
       gsl_multimin_fdfminimizer *s;
     
       /* set some parameters (dummy values not used) */
       double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };
     
       gsl_vector *x;
       gsl_multimin_function_fdf my_func;

       my_func.n = (int)dsize;
       my_func.f = my_f;
       my_func.df = my_df;
       my_func.fdf = my_fdf;
       my_func.params = par;

       fname = (char*)calloc(100,sizeof(char));
     
       // max no. iterations per cycle before positions are recalculated
       iterationlimit = 500;

       /* Starting point */
       x = gsl_vector_alloc (dsize);
       for (i=0; i<dsize; i++)
	 gsl_vector_set (x, i, soln[i]);

       //outputfstamp("start.fits");

       // use the Fletcher-Reeves conjugate gradient method
       // (this seems to work better than Polak-Ribiere provided frequent resets are made)
       // T = gsl_multimin_fdfminimizer_conjugate_fr;

       // use the BFGS method
       T = gsl_multimin_fdfminimizer_vector_bfgs2;
       s = gsl_multimin_fdfminimizer_alloc (T, dsize);
     
       stepsize = 0.1;  // 0.01;
       tol = 0.01;  // 1.e-5;
       gsl_multimin_fdfminimizer_set (s, &my_func, x, stepsize, tol);

       status = GSL_CONTINUE;

       ncall_f = ncall_df = ncall_fdf = 0;
       lchisq = nullchisq;
       do {
           iter++;
	   pchisq = lchisq;
           status = gsl_multimin_fdfminimizer_iterate (s);
	   if (status) break;
	   status = gsl_multimin_test_gradient (s->gradient, 1e-5);
	   lchisq = s->f;
	   deltachisq = pchisq-lchisq;
	   // if (VERBOSE == 1)
	     // printf(" minimisation iteration %d %g %f \n",(int)iter, s->f, deltachisq); fflush(stdout);
	   if (status == GSL_SUCCESS) printf ("Minimum found\n");
     	   for (i=0; i<dsize; i++)
	     {
	       soln[i] = gsl_vector_get(s->x,i);
	     }
	   //sprintf(fname,"iter.%d.fits",(int)iter);
	   //outputfstamp(fname);
       } while (status == GSL_CONTINUE && deltachisq > deltachisqlimit && (int)iter < iterationlimit);

       printf(" number of function calls to f,df,fdf: %d %d %d in %d iterations\n",
	      ncall_f,ncall_df,ncall_fdf,(int)iter); fflush(stdout);
     
       gsl_multimin_fdfminimizer_free (s);
       gsl_vector_free (x);
       free(fname);

       return lchisq;

}


void reconstructosPSF(double xval, double yval, int chip)
{
  // reconstruct the valid PSF region at this location on the mosaic
  // from the polynomial coefficients.

  // the psf[] array returned here has elements which are matched to
  // the solution vector - i.e. elements are not yet ordered by location
  // on the image.  The placement of elements onto a downsampled image
  // is done by functions smatrix and downsample
  
  int i, n, nq, qx, qy, ichip;

  for (i=0; i<nfit; i++)
    {
      nq=0;
      psf[i] = 0.;
      for (qx=0; qx<=order; qx++)
	{
	  for (qy=0; qy<=(order-qx); qy++)
	    {
	      // decide whether this order is chip-dependent or global
	      if (chipvariation==1 && qx+qy <= chiporder)
		{ // chip dependent
		  for (ichip=0; ichip<nchip; ichip++)
		    {
		      n = i*ncoeffs + nq;
		      if (ichip==chip)
			{ // select the appropriate coefficient for this chip
			  psf[i] += (basesoln[n]+soln[n])*pow(xval,qx)*pow(yval,qy);
			}
		      nq++;
		    }
		}
	      else
		{ // not chip dependent
		  n = i*ncoeffs + nq;
		  psf[i] += (basesoln[n]+soln[n])*pow(xval,qx)*pow(yval,qy);
		  nq++;
		}
	    }
	}
    }

}


void downsamplePSFmodel(char *psfname, int xchipsampling, int ychipsampling, int big_gap, 
			int nxchip, int nychip, int *numstar, double *scalefactor,
			double hxsize, double hysize)
{
  // reconstruct the valid PSF region at this location on the mosaic
  // from the polynomial coefficients.

  // the psf[] array returned here has elements which are matched to
  // the solution vector - i.e. elements are not yet ordered by location
  // on the image.  The placement of elements onto a downsampled image
  // is done by functions smatrix and downsample
  
  int i, n, nq, qx, qy, ichip, nset, nqmax, x, y ,p;
  double pos[2];
  double **psfcoeffs;

  // allocate memory for psf coeffs array
  psfcoeffs = (double**)calloc(fsize, sizeof(double*));
  for (i=0; i<fsize; i++)
    {
      psfcoeffs[i] = (double*)calloc(ncoeffs, sizeof(double));
    }

  // create mixing matrix for a centred PSF model
  pos[0]=pos[1]=0.;
  nset = smatrix(pos);

  // split out the solution coefficients (based on reconstructosPSF)
  for (i=0; i<nfit; i++)
    {
      nq=0;
      for (qx=0; qx<=order; qx++)
	{
	  for (qy=0; qy<=(order-qx); qy++)
	    {
	      // decide whether this order is chip-dependent or global
	      if (chipvariation==1 && qx+qy <= chiporder)
		{ // chip dependent
		  for (ichip=0; ichip<nchip; ichip++)
		    {
		      n = i*ncoeffs + nq;
		      psfcoeffs[i][nq] = (basesoln[n]+soln[n]);
		      nq++;
		    }
		}
	      else
		{ // not chip dependent
		  n = i*ncoeffs + nq;
		  psfcoeffs[i][nq] = (basesoln[n]+soln[n]);
		  nq++;
		}
	      if (nq > ncoeffs)
		{
		  fprintf(stderr," error in downsamplePSF, nq = %d ncoeffs = %d \n",nq,ncoeffs);
		  exit(EXIT_FAILURE);
		}
	    }
	}
    }

  // for each coefficient in turn create a pseudo-PSF array
  nqmax = nq;
  if (nqmax != ncoeffs)
    {
      fprintf(stderr, " mismatch in coeffcients count, %d and %d \n",nqmax,ncoeffs);
      exit(EXIT_FAILURE);
    }
  for (nq=0; nq<nqmax; nq++)
    {
      for (i=0; i<nfit; i++)
	{      
	  psf[i] = psfcoeffs[i][nq];
	}
      // apply the downsample algorithm to generate a coefficients image across the postage stamp
      downsample();
      // refill this plane of the psfcoeffs array with quadrant swapping to match lensfit
      for (y=0; y<pheight; y++)
	{
	  //	  yy = y+pheight/2;
	  //	  if (yy>=pheight) yy -= pheight;
	  for (x=0; x<pwidth; x++)
	    {
	      //	      xx = x+pwidth/2;
	      //	      if (xx>=pwidth) xx -= pwidth;
	      p = y*pwidth+x;
	      //	      pp = yy*pwidth+xx;
	      psfcoeffs[p][nq] = dpsf[p];
	    }
	}
    }

  // write out PSF coefficients to FITS cube
  printf(" writing PSF file to %s \n",psfname); fflush(stdout);
 
  correct_distortion = (int)CORRECT_DISTORTION;

  writeglobalcoeffs(psfname, correct_distortion, psfcoeffs, nqmax, order, chipvariation, chiporder, scalefactor,
		    pwidth, pheight, nxchip, nychip, numstar, xchipsampling, ychipsampling,
		    big_gap, hxsize, hysize, poserror);

  // free psfcoeffs memory
  for (i=0; i<fsize; i++)
    {
      free(psfcoeffs[i]);
    }
  free(psfcoeffs);

}


double my_f (const gsl_vector *x, void *params)
{
  //double *fp = (double *)params;
  int k, i, nset;
  double lchisq, amp, modelsum, cross_sum;
      
  ncall_f++;

  for (i=0; i<dsize; i++)
    soln[i] = gsl_vector_get(x, i);

  lchisq = 0.;
  for (k=0; k<nchip; k++)
    {
      for (i=0; i<nstar[k]; i++)
	{
	  // create shifted, downsampled model
	  reconstructosPSF(xval[k][i], yval[k][i], k);
	  nset = smatrix(shift[k][i]);
	  downsample();
	  // get amplitude, chi-squared for this star
	  lchisq += starweight[k][i]*chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum);
	}
     }

  return lchisq;
}
     
/* The gradient of f, df = (df/dx, df/dy). */
void my_df (const gsl_vector *x, void *params, gsl_vector *df)
{
  //double *fp = (double *)params;
  int k, i, j, m, nq, qx, qy, n, p, nset, psize, ichip;
  double noisesq, lchisq, *grad, amp, modelsum, cross_sum;

  ncall_df++;
  //printf (" my_df call %d\n",ncall_df);fflush(stdout);

  psize = pwidth*pheight;
  grad = (double*)calloc(dsize, sizeof(double));

  for (i=0; i<dsize; i++)
    soln[i] = gsl_vector_get(x, i);

  lchisq = 0.;
  for (k=0; k<nchip; k++)
    {
      for (i=0; i<nstar[k]; i++)
	{
	  // create shifted, downsampled model
	  reconstructosPSF(xval[k][i], yval[k][i], k);
	  nset = smatrix(shift[k][i]);
	  downsample();
	  // get amplitude, chi-squared for this star
	  lchisq += starweight[k][i]*chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum);
	  //if (ncall_df==1) printf("chip %d star %d %f %f \n",k,i,modelsum,cross_sum);
	  // use the downsampled model to calculate the gradient
	  for (p=0; p<psize; p++)
	    {
	      noisesq = noise[k]*noise[k] + dpix[k][i][p]/gain[k];
	      noisesq += pow(dpix[k][i][p]*systematic, 2); // add a small fractional error in quadrature
	      // set gradients
	      for (j=0; j<npix[p]; j++)
		{
		  m = spix[p][j];
		  nq=0;
		  for (qx=0; qx<=order; qx++)
		    {
		      for (qy=0; qy<=(order-qx); qy++)
			{
			  /* accumulate gradient */
			  // decide whether this order is chip-dependent or global
			  if (chipvariation==1 && qx+qy <= chiporder)
			    { // chip dependent
			      for (ichip=0; ichip<nchip; ichip++)
				{
				  n = m*ncoeffs + nq;
				  if (ichip==k)
				    { // select the appropriate coefficient for this chip
				      grad[n] += 
					-2.*starweight[k][i]*amp*sval[p][j]*pow(xval[k][i],qx)*pow(yval[k][i],qy)*
					( (dpix[k][i][p]-amp*dpsf[p]/modelsum)*(2.-modelsum)/noisesq +
					  cross_sum*(1.-modelsum)/modelsum );
				    }
				  nq++;
				}
			    }
			  else
			    { // not chip dependent
			      n = m*ncoeffs + nq;
			      grad[n] += 
				-2.*starweight[k][i]*amp*sval[p][j]*pow(xval[k][i],qx)*pow(yval[k][i],qy)*
				( (dpix[k][i][p]-amp*dpsf[p]/modelsum)*(2.-modelsum)/noisesq +
				  cross_sum*(1.-modelsum)/modelsum );
			      nq++;
			    }
			}
		    }
		}
	    }
	}
    }

  for (i=0; i<dsize; i++)
    gsl_vector_set(df, i, grad[i]);

  free(grad);

}
     
/* Compute both f and df together (required by GSL routine but not used!) */
void my_fdf (const gsl_vector *x, void *params, 
             double *lchisq, gsl_vector *df) 
{
  ncall_fdf++;

  *lchisq = my_f(x, params);
  my_df(x, params, df);

}



int smatrix(double *pos)
{
  // create mixing matrices for a given star offset pos[0],pos[1] (native coordinates)
  // allowing linear subpixel interpolation

  int x,xx,xxx,y,yy,yyy,jj,tfsize;
  int xmin,xmax,ymin,ymax,p,fp,nset,setpindex;
  double rxmin,rxmax,rymin,rymax,fx1,fx2,fy1,fy2;

  tfsize = fwidth*fheight;

  nset = 0;
  for (y=0; y<pheight; y++)
    {
      rymin = ofactor*(y+pos[1]-0.5);
      rymax = ofactor*(y+pos[1]+0.5);
      ymin = (int)(rymin+fheight+0.5)-fheight;
      ymax = (int)(rymax+fheight+0.5)-fheight;
      for (x=0; x<pwidth; x++)
	{
	  rxmin = ofactor*(x+pos[0]-0.5);
	  rxmax = ofactor*(x+pos[0]+0.5);
	  xmin = (int)(rxmin+fwidth+0.5)-fwidth;
	  xmax = (int)(rxmax+fwidth+0.5)-fwidth;
	  p = y*pwidth + x;
	  npix[p] = 0;
	  for (yy=ymin; yy<=ymax; yy++)
	    {
	      fy1 = yy+0.5-rymin;
	      if (fy1>1.) fy1=1.;
	      if (fy1<0.) fy1=0.;
	      fy2 = rymax-yy+0.5;
	      if (fy2>1.) fy2=1.;
	      if (fy2<0.) fy2=0.;
	      yyy = yy;
	      if (yyy<0) yyy += fheight;
	      if (yyy>=fheight) yyy -= fheight;
	      for (xx=xmin; xx<=xmax; xx++)
		{
		  fx1 = xx+0.5-rxmin;
		  if (fx1>1.) fx1=1.;
		  if (fx1<0.) fx1=0.;
		  fx2 = rxmax-xx+0.5;
		  if (fx2>1.) fx2=1.;
		  if (fx2<0.) fx2=0.;
		  xxx = xx;
		  if (xxx<0) xxx += fwidth;
		  if (xxx>=fwidth) xxx -= fwidth;		  
		  if (xxx>=0 && xxx<fwidth && yyy>=0 && yyy<fheight)
		    {
		      fp = yyy*fwidth + xxx;
		      if (pindex[fp] >= 0)
			{
			  if (pindex[fp]>=tfsize)
			    {
			      fprintf(stderr," smatrix pindex error at %d, pindex = %d\n",fp,pindex[fp]);
			      fflush(stderr);
			      exit(EXIT_FAILURE);
			    }
			  // determine whether this high-res pixel has been accessed already 
			  // for this low-res pixel, and increment the weight if so
			  setpindex = 0;
			  for (jj=0; jj<npix[p]; jj++)
			    {
			      if (spix[p][jj]==pindex[fp])
				{
				  sval[p][jj] += fx1*fx2*fy1*fy2;
				  setpindex = 1;
				  break;
				}
			    }
			  // if this high-res pixel hasn't been used then add it to the list
			  if (setpindex==0)
			    {
			      spix[p][npix[p]] = pindex[fp];
			      sval[p][npix[p]] = fx1*fx2*fy1*fy2;
			      npix[p]++;
			      nset++;
			    }
			}
		    }
		}
	    }
	}
    }

  return (nset);

}


void crosscorrelate(int range, double *chisq, double *starimage, int chip,
		    double *pos, double *bestamp, double *chisqmin)
{
  int x, y, xmax, ymax, bestp, previousp, dx, dy, xt[2], p, nset;
  int fncoeffs, forder, fcrossterm, i, j, nt, fitrange, bestdx, bestdy;
  double amp, rpos[2], denom, modelsum, cross_sum, lchisq, refchisq;

  // on input pos is the current best guess of the star offset
  // on exit pos is the updated estimate
  // floating point values measured in native pixel coordinates

  // test if the input position is already at an extreme value - if it is, return with no changes
  if (pos[0] <= -pwidth/2 || pos[0] >= pwidth/2 || pos[1] <= -pheight/2 || pos[1] >= pheight/2)
    {
      xmax=ymax=0;
      rpos[0] = pos[0];
      rpos[1] = pos[1];
      nset = smatrix(rpos);
      downsample();
      lchisq = chisqmin[0] = chisquared(starimage, chip, bestamp, &modelsum, &cross_sum);
      //if (printout==1)
      printf(" *** pos out of range, star %d returned values %d %d %d %f %f %f %f\n",
	     current_star,ofactor,xmax,ymax,chisqmin[0],lchisq,pos[0],pos[1]);fflush(stdout);
    } else { // input position is within sensible range

  // set nominal position at pos[0],pos[1] in centre of array
  bestp = fheight*fwidth/2 + fwidth/2;
  previousp = -1; // force iterations to start
  chisqmin[0] = nullchisq;
  for (y=0; y<fheight; y++)
    {
      for (x=0; x<fwidth; x++)
	{
	  p = y*fwidth+x;
	  chisq[p] = 2*nullchisq;
	}
    }
  // iterate until minimum has been found
  while (bestp != previousp)
    {
      previousp = bestp;
      ymax = bestp/fwidth;
      xmax = bestp - ymax*fwidth;
      ymax -= fheight/2;
      xmax -= fwidth/2;
      // printf(" xmax, ymax %d %d \n",xmax,ymax); fflush(stdout);
      // shift to grid of points around nominal position and find lowest chi-squared
      for (dy=-range; dy<=range; dy++)
	{
	  xt[1] = ymax+dy;
	  for (dx=-range; dx<=range; dx++)
	    {
	      xt[0] = xmax+dx;
	      //printf(" xt %d yt %d \n",xt,yt);
	      if (xt[0]>=-fwidth/2 && xt[0]<fwidth/2 && xt[1]>=-fheight/2 && xt[1]<fheight/2)
		{
		  p = (xt[1]+fheight/2)*fwidth + (xt[0]+fwidth/2);
		  //if (printout==1)
		  //printf(" dx, xt, p %d %d %d %d %d \n",dx,dy,xt[0],xt[1],p); fflush(stdout);
		  if (chisq[p] >= nullchisq)
		    {
		      rpos[0] = (double)xt[0]/ofactor + pos[0];
		      rpos[1] = (double)xt[1]/ofactor + pos[1];
		      nset = smatrix(rpos);
		      //if (printout==1) printf(" smatrix\n"); fflush(stdout);
		      downsample();
		      //if (printout==1) printf(" downsample\n"); fflush(stdout);
		      if (xt[0]==0 && xt[1]==0 && writeout==0) 
			{
			  writepsfstamp(pwidth,pheight,dpsf,"teststamp.fits");
			  writeout=1;
			}
		      chisq[p] = chisquared(starimage, chip, &amp, &modelsum, &cross_sum);
		      if (amp<=0.)
			{
			  printf("chisquared chip %d star %d shift %f %f shift %f %f\n",chip+1,current_star,
				 pos[0],pos[1],rpos[0],rpos[1]);
			  fflush(stdout);
			}
		      //if (printout==1)
		      //printf(" chisq vals %d %d %d %f %f %f %f \n",p,xt[0],xt[1],rpos[0],rpos[1],chisq[p],amp);fflush(stdout);
		      if (chisq[p] < chisqmin[0])
			{
			  chisqmin[0] = chisq[p];
			  bestamp[0] = amp;
			  bestp = p;
			}
		    }
		}
	    }
	}
    }

  ymax = bestp/fwidth;
  xmax = bestp - ymax*fwidth;
  ymax -= fheight/2;
  xmax -= fwidth/2;
    
  // obtain more accurate value by fitting quadratic surface to the minimum chi-squared
  forder = 2;
  fncoeffs = 6;
  fcrossterm = 0;
  
  //if (printout==1)
  //printf(" init\n");fflush(stdout);

  for(i = 0; i <= fncoeffs; i++) {
    w[i] = 0.;
    for(j = 0; j <= fncoeffs; j++) {
      u[i][j] = 0.;
      v[i][j] = 0.;
    }
  }

  //if (printout==1)
  //printf(" fitted values\n");fflush(stdout);

  // zoom in with a higher resolution grid for the fitting
  nt=0;
  bestdx = bestdy = 0;
  fitrange = 3;
  refchisq = chisqmin[0];
  for (dy=-fitrange; dy<=fitrange; dy++)
    {
      for (dx=-fitrange; dx<=fitrange; dx++)
	{
	  rpos[0] = ((double)xmax+(double)dx/fitrange)/ofactor + pos[0];
	  rpos[1] = ((double)ymax+(double)dy/fitrange)/ofactor + pos[1];
	  nset = smatrix(rpos);
	  downsample();
	  wfit[nt] = 1.;
	  xfit[nt] = (double)dx/fitrange;
	  yfit[nt] = (double)dy/fitrange;
	  lchisq = chisquared(starimage, chip, &amp, &modelsum, &cross_sum);
	  if (amp<=0.)
	    {
	      printf("chisquared high res chip %d star %d shift %f %f shift %f %f\n",chip+1,current_star,
		     pos[0],pos[1],rpos[0],rpos[1]);
	      fflush(stdout);
	    }
	  zsort[nt] = zfit[nt] = refchisq - lchisq;
	  //if (printout==1)
	  //printf(" %f %f %f \n",rpos[0],rpos[1],zfit[nt]);fflush(stdout);
	  if (lchisq < chisqmin[0])
	    {
	      chisqmin[0] = lchisq;
	      bestdx = dx;
	      bestdy = dy;
	    }
	  nt++;
	}
    }
      
  // set default best-fit position as best-fit discrete value in the above grid search
  pos[0] += ((double)xmax+(double)bestdx/fitrange)/ofactor;
  pos[1] += ((double)ymax+(double)bestdy/fitrange)/ofactor;

  // 2D polynomial fit 
  if (nt>=6)
    {
      if (nt>9)
	{
	  // sort the values and select the nearest to minimum chi-squared
	  qsort(zsort, nt, sizeof(double), dcompare);
	  for (i=0; i<nt; i++)
	    {
	      // recenter the grid to be fitted
	      xfit[i] -= (double)bestdx/fitrange;
	      yfit[i] -= (double)bestdy/fitrange;
	      if (zfit[i] < zsort[nt-9]) 
		{
		  wfit[i]=0.;
		}
	      else
		{
		  //if (printout==1)
		  //printf(" %f %f %f \n",xfit[i],yfit[i],zfit[i]);fflush(stdout);
		}
	    }
	}
      // 2D polynomial fit
      svdfit2dsqc(xfit, yfit, zfit, wfit, nt, forder,
		  fcrossterm, avals, u, v, w);
      // test if it's a sensible fit with a maximum within the specified radius 
      denom = 4. * avals[3] * avals[6] - avals[5] * avals[5];
      if (avals[6] < 0. && avals[3] < 0. && denom > 0.) 
	{
	  // calculate position of maximum in this relative grid but using native pixel sizes
	  rpos[0] = (avals[5] * avals[2] - 2. * avals[4] * avals[3]) / denom / ofactor;
	  rpos[1] = (avals[5] * avals[4] - 2. * avals[2] * avals[6]) / denom / ofactor;
	  // warn if the calculated maximum is outside the measured range of relative position
	  if (fabs(rpos[0]) > 1. || fabs(rpos[1]) > 1.) 
	    {
	      printf(" current star %d \n",current_star);
	      printf(" crosscorrelate warning on calculated position: %f %f \n",rpos[0],rpos[1]);
	      fflush(stdout);
	    }
	  // convert relative position into absolute position in units of native pixels
	  rpos[0] += pos[0];
	  rpos[1] += pos[1];
	  // check chi-squared and update best-fit position if improved
	  nset = smatrix(rpos);
	  downsample();
	  lchisq = chisquared(starimage, chip, &amp, &modelsum, &cross_sum);
	  if (amp<=0.)
	    {
	      printf("chisquared fitted pos chip %d star %d shift %f %f shift %f %f\n",chip+1,current_star,
		     pos[0],pos[1],rpos[0],rpos[1]);
	      for (i=0; i<nt; i++)
		{
		  printf("%f %f %f %f\n",xfit[i],yfit[i],zfit[i],wfit[i]);
		}
	      printf("avals %f %f %f %f %f %f \n",avals[2],avals[3],avals[4],avals[5],avals[6],denom);
	      fflush(stdout);
  	    } 
	  if (lchisq < chisqmin[0])
	    {
	      // update the position only if the fitting has resulted in an improvement
	      // over the best grid position
	      chisqmin[0] = lchisq;
	      pos[0] = rpos[0];
	      pos[1] = rpos[1];
	    } 
	}
    }
  }

  if (printout==1)
    printf(" final values %d %d %d %f %f %f %f\n",ofactor,xmax,ymax,chisqmin[0],lchisq,pos[0],pos[1]);fflush(stdout);

}


double chisquared(double *starimage, int chip, double *amp, double *modelsum, double *cross_sum)
{
  // find chi-squared between model and data with amplitude as a free parameter
  // returns chi-squared and best-fitting amplitude value
  // NB amplitude is defined to be the amplitude in the case where the sum of
  // model values is unity.  This forces the PSF to stay normalised to unity

  int x, y ,p;
  double s1, s2, s3, sm, lchisq, nval;

  //printf(" noise, gain %g %g \n",noise,gain);

  if (printout==1) printf(" entered chisquared\n"); fflush(stdout);

  lchisq = nullchisq;
  amp[0] = -99.;
  s1 = s2 = s3 = sm = 0.;
  for (y=0; y<pheight; y++)
    {
      for (x=0; x<pwidth; x++)
	{
	  p = y*pwidth + x;
	  nval = pow(noise[chip],2) + starimage[p]/gain[chip]; // shot noise squared in ADU 
	  nval += pow(starimage[p]*systematic, 2); // add a small fractional error in quadrature
	  s1 += pow(starimage[p],2)/nval;
	  s2 += pow(dpsf[p],2)/nval;
	  s3 += dpsf[p]*starimage[p]/nval;
	  sm += dpsf[p];
	}
    }
  if (printout==1) printf(" chisquared loop\n"); fflush(stdout);
  if (s2 > 0.)
    {
      amp[0] = s3*sm/s2;
      lchisq = s1 - amp[0]*s3*(2.-sm);
    }
  else
    {
      printf(" error calculating chi-squared, s2 = %g \n",s2);
      printf(" error calculating chi-squared, s1 = %g s3 = %g sm = %g\n",s1,s3,sm);
      printf(" chip %d star %d \n",chip+1, current_star); 
      fflush(stdout);
      lchisq = s1;
      amp[0] = 0.;
    }

  *modelsum = sm;
  *cross_sum = s3;
  if (printout==1) printf(" chisquared done\n"); fflush(stdout);
  return (lchisq);

}


double storepixelvalues()
{
  // function to copy data pixel arrays
  // if bval is non-zero on input, a nonlinearity correction will be applied
  // also calculate a "null model" value of chi-squared

  int k, i, x, y, p;
  double lchisq, nval, pval;

  lchisq = 0.;
  for (k=0; k<nchip; k++)
    {
      for (i=0; i<nstar[k]; i++)
	{
	  for (y=0; y<pheight; y++)
	    {
	      for (x=0; x<pwidth; x++)
		{
		  p = y*pwidth + x;
		  pval = ppix[k][i][p];
		  dpix[k][i][p] = pval - bval*pval*pval;;
		  nval = pow(noise[k],2) + dpix[k][i][p]/gain[k]; // shot noise squared in ADU 
		  nval += pow(dpix[k][i][p]*systematic, 2); // add a small fractional error in quadrature
		  lchisq += starweight[k][i]*pow(dpix[k][i][p],2)/nval;
		}
	    }
	}
    }

  return (lchisq);

}


double nonlinearity()
{
  // function to calculate best estimate of non-linearity correction coefficient

  int i, k, nset, x, y, p;
  double y3, y4, ym, y2m, m2, nval, pval;
  double amp, modelsum, cross_sum, lchisq, linearchisq;

  // printf(" initial non-linearity coefficient %g \n", bval);

  // restore the original pixels with no nonlinearity correction
  for (k=0; k<nchip; k++)
    {
      for (i=0; i<nstar[k]; i++)
	{
	  for (y=0; y<pheight; y++)
	    {
	      for (x=0; x<pwidth; x++)
		{
		  p = y*pwidth + x;
		  dpix[k][i][p] = ppix[k][i][p];
		}
	    }
	}
    }
    
  // accumulate summations for best-fit non-linearity coefficient bval
  y3 = y4 = 0.;
  lchisq = 0.;
  for (k=0; k<nchip; k++)
    {
      for (i=0; i<nstar[k]; i++)
	{
	  // create shifted, downsampled model
	  reconstructosPSF(xval[k][i], yval[k][i], k);
	  nset = smatrix(shift[k][i]);
	  downsample();	 
	  // measure chi-squared for linear model
	  lchisq += starweight[k][i]*chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum);
	  ym = y2m = m2 = 0.;
	  for (y=0; y<pheight; y++)
	    {
	      for (x=0; x<pwidth; x++)
		{
		  p = y*pwidth + x;
		  nval = pow(noise[k],2) + dpix[k][i][p]/gain[k]; // shot noise squared in ADU 
		  nval += pow(dpix[k][i][p]*systematic, 2); // add a small fractional error in quadrature
		  y3 += pow(dpix[k][i][p],3)/nval;
		  y4 += pow(dpix[k][i][p],4)/nval;
		  m2 += pow(dpsf[p],2)/nval;
		  ym += dpsf[p]*dpix[k][i][p]/nval;
		  y2m += dpsf[p]*pow(dpix[k][i][p],2)/nval;
		}
	    }
	  if (m2 > 0.)
	    {
	      y3 -= (ym*y2m/m2);
	      y4 -= (y2m*y2m/m2);
	    }
	  else
	    {
	      fprintf(stderr, " error in nonlinearity, m2 = %g\n",m2);
	      exit(EXIT_FAILURE);
	    }
	}
    }

  if (y4 != 0.)
    {
      bval = y3/y4;
    }
  else
    {
      fprintf(stderr, " error in nonlinearity, y3 = %g y4 = %g \n",y3,y4);
      exit(EXIT_FAILURE);
    }

  // store linear-model chi-squared
  linearchisq = lchisq;

  // apply the new nonlinearity correction and calculate new chi-squared
  if (bval != 0.)
    {
      lchisq = 0.;
      for (k=0; k<nchip; k++)
	{
	  for (i=0; i<nstar[k]; i++)
	    {
	      for (y=0; y<pheight; y++)
		{
		  for (x=0; x<pwidth; x++)
		    {
		      p = y*pwidth + x;
		      pval = dpix[k][i][p];
		      dpix[k][i][p] = pval - bval*pval*pval;
		    }
		}
	      // calculate new chi-squared
	      // create shifted, downsampled model
	      reconstructosPSF(xval[k][i], yval[k][i], k);
	      nset = smatrix(shift[k][i]);
	      downsample();
	      // get amplitude, chi-squared for this star
	      lchisq += starweight[k][i]*chisquared(dpix[k][i], k, &amp, &modelsum, &cross_sum);
	    }
	}
    }
  
  printf(" non-linear correction coefficient %g, chisquared values %g %g %g \n",
	 bval, linearchisq, lchisq, linearchisq-lchisq);

  return(lchisq);

}


void downsample()
{
  int x, y, k, p;

  // sum PSF array into downsampled array using the smatrix values
  for (y=0; y<pheight; y++)
    {
      for (x=0; x<pwidth; x++)
	{
	  p = y*pwidth + x;
	  if (printout==1) printf(" ds %d %d %d \n",x,y,p); fflush(stdout);
	  dpsf[p] = 0.;
	  if (npix[p]>nfitmax)
	    {
	      fprintf(stderr," downsample error in npix for pixel %d, %d > %d\n",p,npix[p],nfitmax);
	      fflush(stderr);
	      exit(EXIT_FAILURE);
	    }
	  for (k=0; k<npix[p]; k++)
	    {
	      if (spix[p][k]<0 || spix[p][k]>=fwidth*fheight)
		{
		  fprintf(stderr," spix error at %d %d, spix = %d\n",p,k,spix[p][k]);
		  fflush(stderr);
		  exit(EXIT_FAILURE);
		}
	      dpsf[p] += sval[p][k]*psf[spix[p][k]];
	    }
	}
    }

}

int initialisePSF()
{
  // define the subregion within which the PSF is nonzero, and fill the coefficients
  // that describes that subregion with initial PSF values, assumed to be an invariant gaussian

  int px, py, dx, dy, x, y, xx, yy, n, p, ichip, osx, osy;
  double rsq, rsqos, rsqoutermax, rsqosmax, ofactorsq, sigmasq, norm, val;

  // define the region within which pixels will be fitted - defined as a radius
  // equal to the half-width of the postage stamp
  rsqoutermax = pwidth*pwidth/4;

  // define also the radius-squared of the region to be oversampled
  rsqosmax = rmax*rmax;

  // square of oversampling
  ofactorsq = ofactor*ofactor;

  // sigma squared of initial trial PSF
  sigmasq = sigma*sigma;

  // printf(" rsqmax, sigmasq %g %g %d %d \n",rsqmax,sigmasq,fheight,fwidth);

  n=0;
  norm=0.; 
  for (py=0; py<pheight; py++)
    {
      yy = py<pheight/2 ? py : py-pheight;
      for (px=0; px<pwidth; px++)
	{
	  xx = px<pwidth/2 ? px : px-pwidth;
	  rsq = xx*xx + yy*yy;
	  for (dy=0; dy<ofactor; dy++)
	    {
	      y = ofactor*py  + dy;
	      for (dx=0; dx<ofactor; dx++)
		{
		  x = ofactor*px + dx;
		  p = y*fwidth + x;
		  // printf(" %d %d %d %d %g %g \n",x,y,xx,yy,rsq,rsqmax);
		  // select pixels to be fitted inside the outer radius
		  if (rsq<=rsqoutermax)
		    {
		      // calculate radius-squared in oversampled units
		      osx = x<fwidth/2 ? x : x-fwidth;
		      osx -= ofactor/2;
		      osy = y<fheight/2 ? y : y-fheight;
		      osy -= ofactor/2;
		      rsqos = osx*osx + osy*osy;
		      // select pixels to be oversampled
		      if (rsq<=rsqosmax) // radius defined at detector sampling
			{
			  // DEFINE OVERSAMPLED REGION
			  // define lookup table of pixels
			  // to be used in fitting
			  pindex[p] = n;
			  // set the pixel values as a gaussian and accumulate the normalisation
			  val = exp(-rsqos/2./sigmasq/ofactorsq); // radius defined in oversampled units
			  norm += val;
			  // fill the constant term of the coefficients 
			  if (chipvariation==1)
			    { // allow chip variations
			      for (ichip=0; ichip<nchip; ichip++)
				{
				  soln[n*ncoeffs+ichip] = val;
				}
			    }
			  else
			    { // no chip variations
			      soln[n*ncoeffs] = val;
			    }
			  n++; // increment fittable pixel number
			}
		      else
			{
			  // DEFINE REGION SAMPLED AT DETECTOR SAMPLING
			  // define lookup table of pixels
			  // to be used in fitting
			  pindex[p] = n; // note that multiple p values will have the same n value
			  // set the pixel values as a gaussian and accumulate the normalisation
			  val = exp(-rsq/2./sigmasq); // radius defined at detector sampling
			  norm += val;
			  if (dx==ofactor-1 && dy==ofactor-1)
			    {
			      // only increment fittable pixel once (i.e. at large radii the sampling is
			      // fixed at detector sampling, at small radii it is oversampled).
			      // fill the constant term of the coefficients 
			      if (chipvariation==1)
				{ // allow chip variations
				  for (ichip=0; ichip<nchip; ichip++)
				    {
				      soln[n*ncoeffs+ichip] = val;
				    }
				}
			      else
				{ // no chip variations
				  soln[n*ncoeffs] = val;
				}
			      n++; // increment fittable pixel number
			    }
			}
		    }
		  else
		    {
		      // outside defined radius, define no pixel
		      // to be used in fitting
		      pindex[p] = -1;
		    }
		}
	    }
	}
    }

  // normalise
  for (p=0; p<n*ncoeffs; p++)
    {
      soln[p] /= norm;
    }

  // printf(" %d %g \n",n,norm);

  return (n);

}


void enablechipvariation(int newchiporder)
{
  // switch on the chip variation by propagating soln coefficients into new array

  int p, newncoeffs, m, nq, oldnq, qx, qy, n, oldn, ichip;

  if (newchiporder < 0)
    {
      printf(" no chip variation enabled\n");
      return;
    }

  for (p=0; p<dsize; p++)
    {
      oldsoln[p] = soln[p];
      soln[p] = 0.;
    }

  newncoeffs = (1+order)*(2+order)/2 + 
    (1+newchiporder)*(2+newchiporder)*(nchip-1)/2;

  for (m=0; m<nfit; m++)
    {
      nq=oldnq=0;
      for (qx=0; qx<=order; qx++)
	{
	  for (qy=0; qy<=(order-qx); qy++)
	    {
	      // decide whether this order is chip-dependent or global
	      if (qx+qy <= newchiporder)
		{ // chip dependent
		  // propagate old solution into solutions for each chip
		  // decide if this was previously chip-dependent or not
		  if (chipvariation==1 && qx+qy <= chiporder)
		    {
		      // was previously chip dependent so copy entire vector
		      for (ichip=0; ichip<nchip; ichip++)
			{
			  oldn = m*ncoeffs + oldnq;
			  oldnq++;
			  n = m*newncoeffs + nq;
			  soln[n] = oldsoln[oldn];
			  nq++;
			}
		    }
		  else
		    {
		      // was not previously chip dependent so propagate
		      oldn = m*ncoeffs + oldnq;
		      oldnq++;
		      for (ichip=0; ichip<nchip; ichip++)
			{
			  n = m*newncoeffs + nq;
			  soln[n] = oldsoln[oldn];
			  nq++;
			}
		    }
		}
	      else
		{ // not chip dependent
		  // propagate coefficients unchanged
		  oldn = m*ncoeffs + oldnq;
		  oldnq++;
		  n = m*newncoeffs + nq;
		  soln[n] = oldsoln[oldn];
		  nq++;
		}
	    }
	}
    }

  // update control variables and array sizes
  chipvariation = 1;
  chiporder = newchiporder;
  ncoeffs = newncoeffs;
  dsize = nfit*ncoeffs;

}


void randomisesoln()
{
  // randomly resample the chip variations to create a random new starting point

  int p;
  double alpha;

  /*
  if (chiporder < 0 || chipvariation==0)
    {
      printf(" no chip variation enabled\n");
      return;
    }

  // copy solution
  for (p=0; p<dsize; p++)
    {
      oldsoln[p] = soln[p];
    }
  */

  // alpha sets the magnitude of the random perturbations
  // experiment indicates alpha=1 works best 
  alpha = 1.;

  // generate a randomly perturbed solution
  for (p=0; p<dsize; p++)
    {
      soln[p] *= (1.+gsl_ran_gaussian(r, alpha));
    }


  /*
  // generate a random resampling of the chips
  for (p=0; p<nchip; p++)
    {
      rchip[p] = (int)gsl_rng_uniform_int(r, (unsigned long int)nchip);
    }


  // shuffle the chip positions in the soln vector
  for (m=0; m<nfit; m++)
    {
      nq=0;
      for (qx=0; qx<=order; qx++)
	{
	  for (qy=0; qy<=(order-qx); qy++)
	    {
	      // decide whether this order is chip-dependent or global
	      if (qx+qy <= chiporder)
		{ // shuffle chip dependent terms
		  n0 = m*ncoeffs + nq;
		  for (ichip=0; ichip<nchip; ichip++)
		    {
		      n = n0 + ichip;
		      oldn = n0 + rchip[ichip];
		      // set the delta-solution to be a fraction alpha of the previous solution
		      delta = alpha*oldsoln[oldn];
		      soln[n] = delta;
		      // if this is a linear term then correct the zero offset term
		      if (qx==1 && qy==0)
			{
			  nzero = m*ncoeffs+ichip;
			  soln[nzero] += delta*(xchipval[rchip[ichip]]-xchipval[ichip]);
			}
		      if (qy==1 && qx==0)
			{
			  nzero = m*ncoeffs+ichip;
			  soln[nzero] += delta*(ychipval[rchip[ichip]]-ychipval[ichip]);
			}
		      nq++;
		    }
		}
	      else
		{ // not chip dependent - set delta-solution to be zero
		  n = m*ncoeffs + nq;
		  soln[n] = 0.;
		  nq++;
		}
	    }
	}
    }
  */

}



int reinitialisePSF()
{
  // redefine the oversampled region and copy the previous solution into the new solution vector

  int i, px, py, dx, dy, x, y, xx, yy, n, p, oldn, oldp, psize;
  double rsq, rsqoutermax, rsqosmax, ofactorsq;

  // define the region within which pixels will be fitted - defined as a radius
  // equal to the half-width of the postage stamp
  rsqoutermax = pwidth*pwidth/4;

  // define also the radius-squared of the region to be oversampled
  rsqosmax = rmax*rmax;

  ofactorsq = ofactor*ofactor;

  // printf(" rsqmax, sigmasq %g %g %d %d \n",rsqmax,sigmasq,fheight,fwidth);

  // remember the old lookup table and solution vector
  psize = pwidth*pheight;
  for (p=0; p<psize; p++)
    {
      oldpindex[p] = pindex[p];
      pindex[p] = -1;
    }
  for (p=0; p<dsize; p++)
    {
      oldsoln[p] = soln[p];
      soln[p] = 0.;
    }

  n=0;
  for (py=0; py<pheight; py++)
    {
      yy = py<pheight/2 ? py : py-pheight;
      for (px=0; px<pwidth; px++)
	{
	  xx = px<pwidth/2 ? px : px-pwidth;
	  rsq = xx*xx + yy*yy;
	  oldp = py*pwidth + px;  // define the old array element assuming ofactor was 1
	  for (dy=0; dy<ofactor; dy++)
	    {
	      y = ofactor*py  + dy - ofactor/2;
	      if (y<0) y += pheight*ofactor;
	      for (dx=0; dx<ofactor; dx++)
		{
		  x = ofactor*px + dx - ofactor/2;
		  if (x<0) x += pwidth*ofactor;
		  p = y*fwidth + x;
		  // select pixels to be fitted inside the outer radius
		  if (rsq<=rsqoutermax)
		    {
		      // select pixels to be oversampled
		      if (rsq<=rsqosmax) // radius defined at detector sampling
			{
			  // DEFINE OVERSAMPLED REGION
			  // define lookup table of pixels
			  // to be used in fitting
			  pindex[p] = n;
			  // identify the old solution components for this low-res pixel
			  // and copy them into the high-res array, suitably rescaled
			  oldn = oldpindex[oldp];
			  for (i=0; i<ncoeffs; i++)
			    {
			      soln[n*ncoeffs+i] = oldsoln[oldn*ncoeffs+i]/ofactorsq;
			    }
			  n++; // increment fittable pixel number
			}
		      else
			{
			  // DEFINE REGION SAMPLED AT DETECTOR SAMPLING
			  // define lookup table of pixels
			  // to be used in fitting
			  pindex[p] = n; // note that multiple p values will have the same n value
			  // identify the old solution components for this low-res pixel
			  // and copy them into the high-res array, suitably rescaled
			  oldn = oldpindex[oldp];
			  if (dx==ofactor-1 && dy==ofactor-1)
			    {
			      for (i=0; i<ncoeffs; i++)
				{
				  soln[n*ncoeffs+i] = oldsoln[oldn*ncoeffs+i]/ofactorsq;
				}
			      n++; // increment fittable pixel number
			    }
			}
		    }
		  else
		    {
		      // outside defined radius, define no pixel
		      // to be used in fitting
		      pindex[p] = -1;
		    }
		}
	    }
	}

    }

  return (n);

}



char *
create_psf_name(const char *psfdir, char *firstname, const char *dotdelimiter)
{
    char *pstr;
    int len;
    char *coeffname;

    //printf(" **i : file handling : %s", psfdir);

    coeffname = (char *) calloc(PATH_MAX, sizeof(char));

    if (!(pstr = strrchr(firstname, *dotdelimiter)))
        pstr = firstname + strlen(firstname);
    sprintf(pstr, "%s", ".psfcoeffs.fits");

    if (psfdir != NULL) {
        strcpy(coeffname, psfdir);
        len = strlen(coeffname);
        if (strncmp(&coeffname[len - 1], "/", 1) != 0) {
            strncat(coeffname, "/", 1);
        }
        if (access(coeffname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr," Can't write psf files to %s \n", coeffname);
            exit(EXIT_FAILURE);
        }
        strcat(coeffname, firstname);
    } else {
        strcpy(coeffname, firstname);
    }
    printf(" creating psf name %s \n", coeffname);

    if (access(coeffname, F_OK) == 0) {
        printf(" *** PSF file exists: removing \n");
        remove(coeffname);
    }
    return coeffname;
}


char *
create_ospsf_name(const char *psfdir, char *firstname, const char *dotdelimiter)
{
    char *pstr;
    int len;
    char *oscoeffname;

    oscoeffname = (char *) calloc(PATH_MAX, sizeof(char));

    if (!(pstr = strrchr(firstname, *dotdelimiter)))
        pstr = firstname + strlen(firstname);
    sprintf(pstr, "%s", ".ospsfcoeffs.fits");

    if (psfdir != NULL) {
        strcpy(oscoeffname, psfdir);
        len = strlen(oscoeffname);
        if (strncmp(&oscoeffname[len - 1], "/", 1) != 0) {
            strncat(oscoeffname, "/", 1);
        }
        if (access(oscoeffname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr," Can't create os psf file %s \n", oscoeffname);
            exit(EXIT_FAILURE);
        }
        strcat(oscoeffname, firstname);
    } else {
        strcpy(oscoeffname, firstname);
    }
    printf(" creating oversampled psf name %s \n", oscoeffname);

    if (access(oscoeffname, F_OK) == 0) {
        printf(" *** oversampled PSF file exists: removing \n");
        remove(oscoeffname);
    }
    return oscoeffname;
}



void
getdata_badpix(char *imagename, char *badpixelname, int *dim, float *apix,
               float *badpix, float *mednoiseval)
{
    fitsfile *afptr, *bfptr;    /* FITS file pointers */
    int status = 0;             /* CFITSIO status value must be initialized to zero */
    int anaxis, bnaxis;
    long anaxes[2] = { 1, 1 }, fpixel[2] = {
    1, 1}, bnaxes[2] = {
    1, 1};
    int size, bsize;
    int i, n, ix, iy, ymin, ymax, xmin, xmax;
    int x, y, medbin;
    float sortval[10201];
    float median, quartile, noiseval[32];

/* open input image */
    fits_open_file(&afptr, imagename, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(EXIT_FAILURE);
    }

/* read dimensions */
    fits_get_img_dim(afptr, &anaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(EXIT_FAILURE);
    }

    if (anaxis != 2) {
        printf
	  ("Error: %d image axes found: images with other than 2 dimensions are not supported\n",anaxis);
        exit(EXIT_FAILURE);
    }

    if (dim[0] != anaxes[0]) {
        fprintf(stderr, " error reading image dimensions \n");
        exit(EXIT_FAILURE);
    }
    if (dim[1] != anaxes[1]) {
        fprintf(stderr, " error reading image dimensions \n");
        exit(EXIT_FAILURE);
    }

    size = dim[0] * dim[1];
    bsize = dim[0] * dim[1];

/* read input data into image array */

    if (fits_read_pix(afptr, TFLOAT, fpixel, size, NULL, apix, NULL, &status)) {
        printf(" error reading pixel data \n");
        exit(EXIT_FAILURE);
    }

/* close main image file */

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(EXIT_FAILURE);
    }


/* repeat for the bad pixel mask (which may not exist) */

/* open input images */
    fits_open_file(&bfptr, badpixelname, READONLY, &status);

    if (status == 0) {

/* read dimensions */
        fits_get_img_dim(bfptr, &bnaxis, &status);
        fits_get_img_size(bfptr, 2, bnaxes, &status);

        if (status) {
            fits_report_error(stderr, status);  /* print error message */
            exit(EXIT_FAILURE);
        }

        if (bnaxis != 2) {
	  printf("Error: %d weight image axes found: images with other than 2 dimensions are not supported\n",bnaxis);
            exit(EXIT_FAILURE);
        }

        if (dim[0] != bnaxes[0]) {
            fprintf(stderr, " error reading bad pixel image dimensions \n");
            exit(EXIT_FAILURE);
        }
        if (dim[1] != bnaxes[1]) {
            fprintf(stderr, " error reading bad pixel image dimensions \n");
            exit(EXIT_FAILURE);
        }

/* read input data into image array */

        if (fits_read_pix(bfptr, TFLOAT, fpixel, bsize, NULL, badpix,
                          NULL, &status)) {
            printf(" error reading bad pixel mask \n");
            fits_report_error(stderr, status);  /* print error message */
            exit(EXIT_FAILURE);
        }

        fits_close_file(bfptr, &status);

        if (status) {
            fits_report_error(stderr, status);  /* print error message */
            exit(EXIT_FAILURE);
        }

    } else {
        if (WEIGHTS_REQUIRED > 0) {
            fprintf(stderr, " bad pixel/weight mask not found \n");
            printf(" bad pixel/weight mask not found %s \n", badpixelname);
            exit(EXIT_FAILURE);
        } else {
            printf(" bad pixel/weight mask not in use, continuing \n");
            for(i = 0; i < bsize; i++)
                badpix[i] = 1.;
            status = 0;
        }
    }

/* sample pixels, fill up an array of values of pixels drawn from
the central region of the image, sort it
and hence find the median and 25-percentile.  If the image has been sky subtracted,
the difference divided by 0.67 will be a good estimate of the rms noise.
If the image has not been subtracted but has a uniform background, this will
still work.  If there are significant background variations, these will contribute
to the noise measurement, making it larger than just the random noise. 
repeat at a number of locations across the image and take the median of these noise
estimates.
*/

    medbin = 0;

    for(iy = 0; iy < 8; iy++) {
        ymin = dim[1] / 16 + (iy * dim[1]) / 8 - 50;
        ymax = ymin + 101;
        if (ymin < 0)
            ymin = 0;
        if (ymax >= dim[1])
            ymax = dim[1] - 1;

        for(ix = 0; ix < 4; ix++) {
            xmin = dim[0] / 8 + (ix * dim[0]) / 4 - 50;
            xmax = xmin + 101;
            if (xmin < 0)
                xmin = 0;
            if (xmax >= dim[0])
                xmax = dim[0] - 1;

            n = 0;
            for(y = ymin; y < ymax; y++) {
                for(x = xmin; x < xmax; x++) {
                    i = y * dim[0] + x;
                    if (badpix[i] > 0.) {
                        sortval[n] = apix[i];
                        n++;
                    }
                }
            }

            if (n > 0) {
                qsort(sortval, n, sizeof(float), compare);
                median = sortval[(n / 2)];
                quartile = sortval[(n / 4)];
                if (median > quartile) {
                    noiseval[medbin] = (median - quartile) / 0.67;
                    medbin++;
                }
            }
        }
    }

    // GREAT3 : get an estimate of the background pixels the rms of pixels 
    //that are on the background, i.e. not near a star.  You could do this 
    // by taking pixels from your long thin image at high y (say y>40) and l
    // ow y (say y<=8) and simply calculate the rms about zero of all those pixels


    // float pixel_sum = 0.;
    // int elem_count = 0;
    // int pheight = 48;
    // int pwidth = 48;
    // int halfwidth, halfheight = 24;
    // int xx, yy, j = 0;
    // ix, iy = 24;

    // for(y = 0; y < pheight; y++) {
    //   yy = y + iy - halfheight;
    //   for(x = 0; x < pwidth; x++) {
    //     xx = x + ix - halfwidth;

    //     j = xx + yy * dim[0];

    //     //take outer pixels
    //     //printf("%d %d %f\n", x, y, apix[j]);
    //     if (( y < 5 )  || (y > 43) && ( x < 5 )  || (x > 43)){
    //       //printf("  **i : taking point at (%d, %d) : %f\n", x, y, apix[j] );
    //       pixel_sum = pixel_sum + (apix[j] * apix[j]);
    //       elem_count++;
    //     }                 
    //   }
    // }
    // printf("\n");
    // float noise_level = 0.;

    // noise_level = sqrt((1.0/elem_count)*pixel_sum);
    // printf("  **i : noise level : %f, elem count : %d\n", noise_level, elem_count);
    // *mednoiseval = noise_level;
    // ------------------------------------------

    if (medbin > 0) {
        qsort(noiseval, medbin, sizeof(float), compare);
        *mednoiseval = noiseval[(medbin / 2)];
    } else {
        *mednoiseval = 0.;
    }

    //printf("  **i : mednoiseval : %f\n", noiseval[(medbin / 2)]);

    status = 0;
    return;

}


void dilatemask(float *badpix, int *dim)
{
  int x, y, i, j;
  
  /* dilate the bad pixel mask by one pixel in all directions - takes care
     of inaccuracies in the bad pixel mask */
  for (y=0; y<dim[1]; y++)
    {
      for (x=1; x<dim[0]; x++)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      badpix[i-1]=0.;
	    }
	}
      for (x=dim[0]-2; x>=0; x--)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      badpix[i+1]=0.;
	    }
	}
    }
  for (x=0; y<dim[0]; x++)
    {
      for (y=1; y<dim[1]; y++)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      j = x + (y-1)*dim[0];
	      badpix[j]=0.;
	    }
	}
      for (y=dim[1]-2; y>=0; y--)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      j = x + (y+1)*dim[0];
	      badpix[j]=0.;
	    }
	}
    }
}


int open3Dfits(char *psfname, fitsfile* *afptr, int num)
{
  int status = 0;
  int anaxis;
  long anaxes[3];
  int bitpix;

  // create a new FITS cube

  bitpix = DOUBLE_IMG;
  anaxis = 3;
  anaxes[0] = pwidth;
  anaxes[1] = pheight;
  anaxes[2] = num;

  printf(" creating 3D FITS imagen %s size %d %d %d \n",
	 psfname, pwidth, pheight, num);
  fflush(stdout);

  fits_create_file(afptr, psfname, &status); 
  fits_create_img(*afptr, bitpix, anaxis, anaxes, &status);

  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }

  return 0;

}


int write3Dfits(fitsfile *afptr, int num, double *f)
{
  int status = 0;
  long fpixel[3]={1,1,1};
  int psfsize;
  double *sf;
  int x, y, xx, yy, ipix, opix;

  psfsize = pwidth*pheight;

  // swap quadrants
  sf = (double*)calloc(psfsize, sizeof(double));
  for (y=0; y<pheight; y++)
    {
      yy = y + pheight/2;
      if (yy >= pheight) yy -= pheight;
      for (x=0; x<pwidth; x++)
	{
	  xx = x + pwidth/2;
	  if (xx >= pwidth) xx -= pwidth;
	  ipix = y*pwidth+x;
	  opix = yy*pwidth+xx;
	  sf[opix] = f[ipix];
	}
    }

  /* write swapped data into cube at correct location */
  fpixel[2] = num+1;

  if (fits_write_pix(afptr, TDOUBLE, fpixel, psfsize, sf, &status) )
    {
      fprintf(stderr," error writing pixel data into 3D cube \n");
    }
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    exit(2);
  }

  free(sf);
  
  return 0;

}


