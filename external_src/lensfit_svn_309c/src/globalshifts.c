/*--------------------------------------------------------------------------

globalshifts

measures offsets between stars and PSF model

based on makeglobalpsf.  LM 1 Nov 2009

bug fix and align with lensfit v7, 17 Nov 2010.

Revision $Rev: 307 $ 
last changed by $Author: miller $
on $LastChangedDate: 2017-08-11 16:28:35 +0100 (Fri, 11 Aug 2017) $

--------------------------------------------------------------------------- */

#define DATE_STRING "$LastChangedDate: 2017-08-11 16:28:35 +0100 (Fri, 11 Aug 2017) $"
#define REV_STRING "$Rev: 307 $"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

/* swarp wcs include files (but note custom version of fitscat.h to
   avoid conflict with cfitsio */
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

#include <complex.h>
#include <fftw3.h>

void psfmoments(double *, int, int, double *, double *, double *);
int updateglobalcoeffs(char *, int, int, int, int, double *, double *);
void globalsvdfit(double *, double *, double *, double *, int *, int, int,
                  int, int, int, int, int, double *, double **, double **,
                  double *);
int getcoeffs_fits_expand_v5(char *, double **, int, int, int, int, int,
                             int *);
int getglobalpsfsize(char *, int *, int *, int *, int *, int *, int *, int *, int *, int *,
                     int *, int *, int *, int *, double *, double *, float *);
int extractpostagestamp(float *, float *, int *, double *, double *, float *,
                        float *, float, float, float, float, float, int, int,
                        float, int *, float);
int swarpextract(wcsstruct *, float *, float *, int *, double *, double *,
                 float, float, float, int, int, ikernelstruct *, float, int *,
                 float, double *);
int tbswarpextract_nodistortion(wcsstruct *, float *, float *, float *, int *,
                                double *, double *, float, float,
                                float, int, int, float,
                                int *, double *, float, double *, double *,
                                double *, double *, int);
int readpsfcatsize(char *);     // read size of catalogue
/* choose appropriate routine for catalogue format */
//int readpsfcat (char*, float*, float*, float*);
//int readpsfcat2_wcs (char*, double*, double*, float*);
int readpsfcat_simple_wcs(char *, double *, double *, float *);
int getimagesize(char *, char *, int *, float *, float *, float *, float *,
                 int *, int *, float);
void getdata_badpix(char *, char *, int *, float *, float *, float *);
int weightfilter(float *, float *, float *, float *, int *, int);
void dilatemask(float *, int*);
int xcorr_measure(fftw_plan, fftw_complex *, fftw_complex *,
                  int, int, int, int,
                  fftw_complex *, double *, fftw_complex *, double *);
int shiftft(int, int, fftw_complex *, double *);
void globalreconstruct(double, double, int, int, int, int, int, int, double *,
                       double *);
int varylpthresholdf2(float *, int *, float, int, int, int, int, int, float,
                      float, int);
void mediansub(float *, float *, int *, float, float *);
int wcs_to_raw(wcsstruct *, double *, double *);
int raw_to_wcs(wcsstruct *, double *, double *);
void load_dfield(fieldstruct *, catstruct *, char *, int, int);
ikernelstruct *init_ikernel(interpenum *, int);

static int
compare(const void *ii, const void *jj)
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

#define VERBOSE 1

/* set whether or not to subtract a constant median background 
   (0=no 1=yes) */
#define SUBTRACT_MEDIAN 0

/* filter weight images (as in CFHTLS analysis) */
#define FILTER_WEIGHTS 0

/* specify that weight (bad pixel) images are required */
#define WEIGHTS_REQUIRED 1

// specify whether to filter out stars based on having a high maximum value (=1) or not (=0)
#define TEST_MAXIMUM 0

/* set minimum number of stars required for subregions to be
   considered valid (only used for warning) (regions with numbers
   equal to this are considered acceptable) 
#define PSF_STAR_LIMIT 1
*/

/* define whether input galaxy catalogue positions are in WCS (WCS=1)
   or pixel (WCS=0) coordinates */
#define WCS 1

#define USE_SWARP 1

/* define whether input files are gzipped fits or plain fits
   0=not zipped, 1=weight files are zipped, 2=both data and weight files zipped */
#define GZIPPED 0

/* define minimum accepted exposure time */
#define minexposuretime 0.

/*---------------------------------------------------------------------------*/
/* main program */

int
main(int argc, char *argv[])
{
  // output version information based on SVN info
    char version[100];
    bzero(version,100);
    strcpy(version, REV_STRING);
    int len; 
    len = strlen(version);
    strcpy(&version[len-2],"\0");
    printf("\n %s version SVN:%s",argv[0],&version[len-5]);
    char datedelims[]="()";
    char *dstring = NULL;
    char date_string[200];
    bzero(date_string,200);
    strcpy(date_string, DATE_STRING);
    dstring  = strtok(date_string, datedelims);
    if (dstring != NULL) dstring = strtok(NULL, datedelims);
    if (dstring != NULL) printf(" %s ",dstring);
    printf("\n");
  
    int flagbad;
    int nchanged;
    int halfpwidth, halfpheight, dim[2], badccd, imageid;
    int padwidth, padheight, halfpadwidth, cchwidth;
    int chipvariation, chiporder;
    int iter;
    int CORRECT_DISTORTION;
    int i, j, ii;
    int x, y;
    int xmin, xmax, ymin, ymax;
    int pheight, pwidth, imagesize;
    int nobj, nobjt, nobj2;
    int **listnum;              /* numbers of good and bad stars in each psf */
    int nstars;
    int ix, iy, pixel, xx, yy;
    int num_off, num_wrong_mag, num_saturated, num_snratio;
    int num_badpix, ncross;
    int ngoodpix, imax_intensity;
    int *region, nobjects;
    int lower_area_limit, lower_merging_limit, arealimit;
    int *offscale;
    int naxis = 2;
    int narg = 0;
    int fheight = 20;  // half-height of median filter in weightfilter
    int image, nchip, xchipsampling, ychipsampling, big_gap, nchip_expected, npsfchip;
    int *xchip, *ychip, ichip;
    int status = { 0 };
    int nn, *chip;
    int sorder, schipvariation, schiporder, scrossterm, sncoeffs;

    time_t t1;
    float memory, memloop;

    double *wfits, *xfits, *yfits, *xshiftfits, *yshiftfits, *ow;
    double *shift;
    double sum, sumx, sumy, sumxsq, sumysq, xrms, yrms, rdx, rdy;
    double psfe[2], centroid[2], moments[3];
    double psfauto, dauto, cross, norm, dnorm, pnorm;
    double xval, yval;
    double xsize, ysize, hxsize, hysize;
    double *radegs, *decdegs;   // ra,dec in input catalogue
    double *rdegs, *ddegs;      // ra,dec arrays after stars off-image have been eliminated  
    double *rx, *ry;
    double wcsneg[2], wcspos[2], rawpos[2], wcscentre[2], wtest[2], wcstest, wcstesttol;
    double scalefactor[1];
    double xs, ys;
    double **u, **v, *w, *avals;

    float poserror, maxlevel;
    float noise;
    float *opix, *changed, *weightfilterarray;
    float *apix, *temp, *badpix, *badtemp;
    float *sortarray, *sn;
    float fintensity_limit, fmax_intensity, scale;
    float *objx, *objy, *mag, *rmag;
    float brightmag, faintmag, snratio;
    float satlev, gain, satlimit, arcperpix, angle;
    float ximageoffset, yimageoffset;
    double *dobjx, *dobjy;
    int order, crossterm, ncoeffs, maxncoeffs;

    double **acoeffs, *xcoeffs, *ycoeffs;
    double *c, *shiftpix, **diff, *num_obj_used;
    double diffmin, diffmax, sx, sy;
    double *bpix, *dpix, *dbadpix;
    double maxvalue, psfnorm, noisevalue;
    fftw_complex *Bpix, *C, *Dpix, *Cpad, *shiftpixft;
    fftw_plan starplan, psfplan, iqplan, padinv;

    char *evariable, *headevariable, *badpixdir, *psfdir;
    char *prefsname;
    char *catname;
    char **image_file;
    char *imagename, *headerfile;
    char *psfname, *coeffname, *fittedpsfname, *imageweightname;
    char *firstname, *badname, *ellname, *usedstarsname, *stampname, *shiftsname;
    char **argkey, **argval;
    char *pstr, delimiter[5], headdelimiter[5], *satname, *weight_suffix;

    FILE *filep, *ellfile, *shiftsfile;

    //int NCHIPS = 1; // default number of CCDs in mosaic

    double pi = M_PI;

    // camera and waveband choices
    int kids,cfht,single,suprime,rband,iband;
    char *cameraenv, *wavebandenv;

    // get the camera and waveband environment variables to be used later for setting parameters and string values
    cameraenv = getenv("CAMERA");
    if (cameraenv != NULL)
      {
        kids = cfht = single = suprime = 0;
        if (strncmp(cameraenv,"kids",1)==0 || strncmp(cameraenv,"KIDS",1)==0) {
          kids = 1;
          printf(" setting up for KIDS (Omegacam) camera\n");
        } else {
          if (strncmp(cameraenv,"cfht",4)==0 || strncmp(cameraenv,"CFHT",4)==0 ||
              strncmp(cameraenv,"mega",4)==0 || strncmp(cameraenv,"MEGA",4)==0 ) 
            {
            cfht = 1;
            printf(" setting up for CFHT (Megacam/MegaPrime) camera\n");
          } else {
            if (strncmp(cameraenv,"single",2)==0 || strncmp(cameraenv,"SINGLE",2)==0) {
              single = 1;
              printf(" setting up for single chip camera\n");
            } else {
              if (strncmp(cameraenv,"suprime",2)==0 || strncmp(cameraenv,"SUPRIME",2)==0) {
                suprime = 1;
                printf(" setting up for Suprime camera\n");
              } else {
                fflush(stdout);
                fprintf(stderr," Camera not recognised in environment variable CAMERA: %s\n",cameraenv);
                exit(EXIT_FAILURE);
              }
            }
          }
        }
      } else {
      fflush(stdout);
      fprintf(stderr," CAMERA environment variable not set \n");
      exit(EXIT_FAILURE);
    }

    wavebandenv = getenv("WAVEBAND");
    if (wavebandenv != NULL) 
      {
	rband = iband = 0;
	if (strncmp(wavebandenv,"R",1)==0 || strncmp(wavebandenv,"r",1)==0)
	  {
	    rband = 1;
	    printf(" waveband set to R\n");
	  } else {
	  if (strncmp(wavebandenv,"I",1)==0 || strncmp(wavebandenv,"i",1)==0)
	    {
	      iband = 1;
	      printf(" waveband set to I\n");
	    } else {
	    fflush(stdout);
	    fprintf(stderr," waveband not recognised: %s \n",wavebandenv);
	    exit(EXIT_FAILURE);
	  }
	}
      } else {
      fflush(stdout);
      fprintf(stderr," WAVEBAND environment variable not set \n");
      exit(EXIT_FAILURE);
    }

    char delims[3], delims2[3];
    bzero(delims,3);
    bzero(delims2,3);
    if (kids==1)
      {
        strcpy(delims,"_");
        strcpy(delims2,"O");
      }
    if (cfht==1 || single==1)
      {
        strcpy(delims,"_C");
        strcpy(delims2,"_C");
      }
    if (suprime==1)
      {
        strcpy(delims,"_O");
        strcpy(delims2,"_O");
      }
    char *item = NULL;
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

    // specify WCS test tolerance when checking for multivalued astrometric solutions
    wcstesttol = pow((1./60.),2); // set to 1 arcmin squared
   
    t1 = time(NULL);

    if (argc != 3 && argc != 5) {
        printf(" %s <imagelist> <snratio> [<bright mag> <faint mag>] \n",
               argv[0]);
        exit(2);
    }

    if (USE_SWARP == 1) {
        printf(" swarp astrometric corrections will be applied \n");
        if (WCS == 0) {
            printf
                (" WARNING:  WCS flag not set, using pixel object coordinates \n");
            printf
                (" even though WCS information is needed for the swarp correction \n");
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
	  fflush(stdout);
            fprintf(stderr,
                    "configuration file environment variable SWARP_CONFIG not set \n");
            exit(2);
        }
        free(argkey);
        free(argval);
        /*
         * initialise the convolution kernel 
         */
        ikernel = init_ikernel(prefs.resamp_type, naxis);
    } else {
        printf(" no swarp astrometric corrections will be applied \n");
    }


    /*
     * test whether magnitude limits for PSF stars have been specified
     * and read values if so, or else set default values.  NB PSF star
     * magnitude values need to be given in the input PSF catalogue in
     * the first case!
     */

    if (argc == 5) {
        brightmag = atof(argv[3]);
        faintmag = atof(argv[4]);
        printf(" selecting PSF stars with %8.2f < m < %8.2f \n", brightmag,
               faintmag);
    } else {
        brightmag = -1.e10;
        faintmag = 1.e10;
    }

    snratio = atof(argv[2]);
    if (snratio <= 0.)
        snratio = 0.;
    printf(" selecting stars with peak S/N ratio > %8.1f \n", snratio);
    if (snratio < 20.)
        printf(" WARNING: low S/N ratio not advised, choose S/N > 20 \n");

    /*
    if (WCS == 1) {
        printf(" reading WCS coordinates from input catalogue\n");
    } else {
        printf(" reading xy coordinates from input catalogue\n");
    }
    */

    /*
     * initialise memory counters 
     */
    memory = 0.;
    memloop = 0;


    evariable = (char *) calloc(300, sizeof(char));
    headevariable = (char *) calloc(300, sizeof(char));
    badpixdir = (char *) calloc(300, sizeof(char));
    psfdir = (char *) calloc(300, sizeof(char));
    prefsname = (char *) calloc(300, sizeof(char));
    catname = (char *) calloc(500, sizeof(char));
    imagename = (char *) calloc(500, sizeof(char));
    headerfile = (char *) calloc(500, sizeof(char));
    psfname = (char *) calloc(500, sizeof(char));
    fittedpsfname = (char *) calloc(500, sizeof(char));
    usedstarsname = (char *) calloc(500, sizeof(char));
    coeffname = (char *) calloc(500, sizeof(char));
    firstname = (char *) calloc(500, sizeof(char));
    imageweightname = (char *) calloc(500, sizeof(char));
    badname = (char *) calloc(500, sizeof(char));
    stampname = (char *) calloc(500, sizeof(char));
    weight_suffix = (char *) calloc(100, sizeof(char));
    satname = (char *) calloc(30, sizeof(char));
    memory += 5300. * sizeof(char);

    if (evariable == NULL ||
        catname == NULL ||
        imagename == NULL ||
        badname == NULL || psfname == NULL || coeffname == NULL) {
	  fflush(stdout);
	  fprintf(stderr,"Memory allocation error for filename strings\n");
        exit(2);
    }

    strcpy(stampname, argv[1]);
    if (!(pstr = strrchr(stampname, *dotdelimiter)))
        pstr = stampname + strlen(stampname);
    sprintf(pstr, "%s", ".diff.fits");

    /*
     * create output PSF filename 
     */

    psfdir = getenv("PSF_DIR");

    strcpy(firstname, argv[1]);
    if (!(pstr = strrchr(firstname, *dotdelimiter)))
        pstr = firstname + strlen(firstname);
    sprintf(pstr, "%s", ".psfcoeffs.fits");

    ellname = ml_calloc(PATH_MAX, sizeof(char), &memory, "ellname");
    bzero(ellname, PATH_MAX);
    strcpy(ellname, firstname);
    strcat(ellname, "_ellipticities.log");
    if ((ellfile = fopen(ellname, "w")) == NULL) {
      fprintf(stderr, "failed to make file %s \n", ellname);
    }

    shiftsname = ml_calloc(PATH_MAX, sizeof(char), &memory, "shiftsname");
    bzero(shiftsname, PATH_MAX);
    strcpy(shiftsname, firstname);
    strcat(shiftsname, "_globalshifts.log");
    if ((shiftsfile = fopen(shiftsname, "w")) == NULL) {
      fprintf(stderr, "failed to make file %s \n", shiftsname);
    }

    
    if (psfdir != NULL) {
        strcpy(coeffname, psfdir);
        len = strlen(coeffname);
        if (strncmp(&coeffname[len - 1], "/", 1) != 0) {
            strncat(coeffname, "/", 1);
        }
        if (access(coeffname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr," Can't access %s \n", coeffname);
            exit(1);
        }
        strcat(coeffname, firstname);
    } else {
        strcpy(coeffname, firstname);
    }
    if (VERBOSE == 1)
      printf(" looking for psf name %s \n", coeffname);

    if (access(coeffname, F_OK) != 0) {
      fflush(stdout);
      fprintf(stderr," psf coefficients file %s does not exist in this location \n",coeffname);
        exit(2);
    }

    /*
     * allocate dummy pointer for getcoeffs_fits_expand 
     */
    listnum = (int **) calloc(1, sizeof(int *));

    xchip = (int*)calloc(100, sizeof(int));
    ychip = (int*)calloc(100, sizeof(int));
    npsfchip = 0;
    
    /*
     * get the information from the psf coefficients file 
     */
    getglobalpsfsize(coeffname, &CORRECT_DISTORTION, &pwidth, &pheight, &order,
                     &chipvariation, &chiporder, &ncoeffs,
                     &npsfchip, xchip, ychip, &xchipsampling, &ychipsampling,
                     &big_gap, &hxsize, &hysize, &poserror);
    if (poserror<=0.)
      {
	fflush(stdout);
	fprintf(stderr,"keyword  poserror not found or invalid in psf file \n");
	exit(EXIT_FAILURE);
      }
    /*
     * (assume survey geometry same for all chips, don't bother checking) 
     */

    if (USE_SWARP != 0 && CORRECT_DISTORTION == 1)
      {
	if (scalefactor[0] <= 0.) {
	  fflush(stdout);
	  fprintf(stderr, " error in returned scalefactor %f \n",
		  scalefactor[0]);
	  exit(2);
	}
      }

    /*
     * check number of chips agrees though 
     */
    /*
    if ((int) NCHIPS != nxchip * nychip) {
        fflush(stdout);
        fprintf(stderr,
                " mismatch in assumed survey geometry and PSF-file geometry\n");
        fprintf(stderr, " %d x %d != %d \n", nxchip, nychip, (int) NCHIPS);
        fflush(stderr);
        exit(2);
    }
    */

    /*
     * set some dimensions 
     */

    halfpwidth = pwidth / 2 + 1;
    halfpheight = pheight / 2;
    padheight = 100 * pheight;
    padwidth = 100 * pwidth;
    halfpadwidth = 1 + padwidth / 2;

    if (padwidth <= pwidth || padheight <= pheight) {
	  fflush(stdout);
        fprintf
	  (stderr," padded cross-correlation array size must be bigger than input size \n");
        exit(1);
    }

        /*
         * define the string delimiter used to allow the code to
         * create the filename for the .head file, in the case where
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
	if (kids==1) {
	  if (rband==1) {
	    strcpy(weight_suffix, "OFCSI.weight.fits");
	  } else {
	    strcpy(weight_suffix, "OFCSF.weight.fits");
	  }
	  strcpy(delimiter,"O");
	  strcpy(headdelimiter,"O");
	} else {
	  if (cfht==1) {
	    strcpy(weight_suffix, "C.weight.fits");
	    strcpy(delimiter,"C");
	    strcpy(headdelimiter,"C");
	  } else {
	    if (single==1) {
	      strcpy(weight_suffix, "C.weight.fits");
	      strcpy(delimiter,"C");
	      strcpy(headdelimiter,"C");
	    } else {
	      if (suprime==1) {
		strcpy(weight_suffix, "OFCS.weight.fits");
		strcpy(delimiter,"O");
		strcpy(headdelimiter,"O");
	      } else {
		fflush(stdout);
		fprintf(stderr," camera not set \n");
		exit(EXIT_FAILURE);
	      }
	    }
	  }
	}
      }
    else
      {
        strcpy(weight_suffix, ".weight.fits");
        strcpy(delimiter,".");
        strcpy(headdelimiter,".");
      }

    nchip_expected = npsfchip;   // total number of chips
    image_file = (char **) calloc(nchip_expected, sizeof(char *));

    /*
     * read the input file list 
     */
    filep = fopen(argv[1], "r");
    if (filep == NULL) {
	  fflush(stdout);
        fprintf(stderr, " Error opening file %s \n", argv[1]);
        exit(2);
    }

    nchip = 0;
    while (!feof(filep) && !ferror(filep)) {
        if (fscanf(filep, "%s", firstname) != 0) {
            if (!feof(filep) && !ferror(filep)) {
                if (nchip >= nchip_expected) {
		  fflush(stdout);
                    fprintf(stderr, " too many images in file \n");
                    exit(2);
                }
                image_file[nchip] = (char *) calloc(500, sizeof(char));
                strcpy(image_file[nchip], firstname);
                nchip++;
            }
        }
    }
    fclose(filep);
    status = 0;
    bzero(firstname, 500);
    if (VERBOSE == 1)
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
      }

    fflush(stdout);

    /*
     * get the order of the polynomial fit.  negative order is taken to mean
     * that cross-terms should be included 
     */
    get_poly_order(order, chiporder, chipvariation, nchip, &ncoeffs,
                   &crossterm);

    acoeffs = (double **) calloc((pwidth * pheight), sizeof(double *));
    memory += (float) (pwidth * pheight) * sizeof(double *);
    if (acoeffs == NULL) {
        printf("Memory allocation error for sub-image pointers\n");
        exit(2);
    }

    for(i = 0; i < (pwidth * pheight); i++) {
        acoeffs[i] = (double *) calloc((1 + ncoeffs), sizeof(double));
        memory += (float) (1 + ncoeffs) * sizeof(double);
        if (acoeffs[i] == NULL) {
            printf("Memory allocation error for sub-image \n");
            exit(2);
        }
    }

    /*
     * get PSF coefficients (listnum is not set so is skipped over
     * inside this routine)
     */
    getcoeffs_fits_expand_v5(coeffname, acoeffs, pwidth, pheight,
                             pwidth, pheight, ncoeffs, listnum[0]);

    /*
     * allocate memory for FFTW arrays etc 
     */

    Bpix =
        (fftw_complex *) calloc((halfpwidth * pheight), sizeof(fftw_complex));
    memory += (float) (halfpwidth * pheight) * sizeof(fftw_complex);
    Dpix =
        (fftw_complex *) calloc((halfpwidth * pheight), sizeof(fftw_complex));
    memory += (float) (halfpwidth * pheight) * sizeof(fftw_complex);
    if (Dpix == NULL) {
        printf("Memory allocation error for Dpix \n");
        exit(2);
    }
    dpix = (double *) calloc((pwidth * pheight), sizeof(double));
    memory += (float) (pwidth * pheight) * sizeof(double);
    /*
     * printf(" Allocated dpix for image %d, %d\n",i,pwidth*pheight);
     */
    if (dpix == NULL) {
        printf("Memory allocation error for sub-image \n");
        exit(2);
    }
    shiftpixft =
        (fftw_complex *) calloc((halfpwidth * pheight), sizeof(fftw_complex));
    memory += (float) (halfpwidth * pheight) * sizeof(fftw_complex);
    bpix = (double *) calloc((pwidth * pheight), sizeof(double));
    memory += (float) (pwidth * pheight) * sizeof(double);
    diff = (double **) calloc(nchip, sizeof(double *));
    memory += (float) nchip *sizeof(double *);
    num_obj_used = (double *) calloc(nchip, sizeof(double));
    memory += (float) nchip *sizeof(double);
    for(i = 0; i < nchip; i++) {
        diff[i] = (double *) calloc((pwidth * pheight), sizeof(double));
        memory += (float) (pwidth * pheight) * sizeof(double);
    }
    shiftpix = (double *) calloc((pwidth * pheight), sizeof(double));
    memory += (float) (pwidth * pheight) * sizeof(double);
    C = (fftw_complex *) calloc((halfpwidth * pheight), sizeof(fftw_complex));
    memory += (float) (halfpwidth * pheight) * sizeof(fftw_complex);
    c = (double *) calloc((padwidth * padheight), sizeof(double));
    memory += (float) (padwidth * padheight) * sizeof(double);
    Cpad =
        (fftw_complex *) calloc((halfpadwidth * padheight),
                                sizeof(fftw_complex));
    memory += (float) (halfpwidth * pheight) * sizeof(fftw_complex);
    if (Bpix == NULL || bpix == NULL || C == NULL || c == NULL
        || Cpad == NULL) {
        printf("memory allocation error for fftw arrays \n");
        exit(2);
    }

    /*
     * temp array is used as spare storage inside extractdata 
     * for quadrant swapping when the postage stamps are extracted 
     */
    temp = (float *) calloc((pwidth * pheight), sizeof(float));
    memory += (float) (pwidth * pheight) * sizeof(float);

    if (temp == NULL) {
        printf("Memory allocation error for temporary image\n");
        exit(2);
    }

    badtemp = (float *) calloc((pwidth * pheight), sizeof(float));
    memory += (float) (pwidth * pheight) * sizeof(float);
    if (badtemp == NULL) {
        printf("Memory allocation error for temporary image\n");
        exit(2);
    }

    shift = (double *) calloc(2, sizeof(double));
    memory += 2. * sizeof(double);
    if (shift == NULL) {
        printf("Memory allocation error for cross-correlation shift array\n");
        exit(2);
    }

    /*
     * create FFTW plans 
     */

    psfplan = fftw_plan_dft_r2c_2d(pwidth, pheight, bpix, Bpix, FFTW_MEASURE);
    starplan =
        fftw_plan_dft_r2c_2d(pwidth, pheight, dpix, Dpix, FFTW_MEASURE);
    iqplan =
        fftw_plan_dft_c2r_2d(pwidth, pheight, shiftpixft, shiftpix,
                             FFTW_MEASURE);
    padinv = fftw_plan_dft_c2r_2d(padwidth, padheight, Cpad, c, FFTW_MEASURE);

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
	  fflush(stdout);
        fprintf(stderr,
                " catalogue name not set, use environment variable CATALOGUE_STARS \n");
        exit(EXIT_FAILURE);
    }

    if (access(catname, F_OK) != 0) {
	  fflush(stdout);
        printf(" Can't read catalogue %s \n", catname);
        exit(EXIT_FAILURE);
    }

    /*
     * read catalogue of objects 
     */

    nobjt = readpsfcatsize(catname);

    objx = (float *) calloc(nobjt, sizeof(float));
    objy = (float *) calloc(nobjt, sizeof(float));
    mag = (float *) calloc(nobjt, sizeof(float));
    memory += (float) (3 * nobjt) * sizeof(float);

    if (WCS == 0) {
	  fflush(stdout);
        fprintf(stderr,
                " x,y star coordinate input not supported in this version \n");
        exit(2);
    } else {
        /*
         * case where world coordinates have been supplied 
         */
      if (VERBOSE == 1)
        printf(" reading world coordinates from input catalogue\n");

        dobjx = (double *) calloc(nobjt, sizeof(double));
        dobjy = (double *) calloc(nobjt, sizeof(double));
        radegs = (double *) calloc(nobjt, sizeof(double));
        decdegs = (double *) calloc(nobjt, sizeof(double));
        rdegs = (double *) calloc(nobjt, sizeof(double));
        ddegs = (double *) calloc(nobjt, sizeof(double));
        rx = (double *) calloc(nobjt, sizeof(double));
        ry = (double *) calloc(nobjt, sizeof(double));
        offscale = (int *) calloc(nobjt, sizeof(int));
        rmag = (float *) calloc(nobjt, sizeof(float));
        sn = (float *) calloc(nobjt, sizeof(float));
        memory += (float) (7 * nobjt) * sizeof(double);

        /*
         * read in coordinates assuming in decimal WCS 
         */
        nobj2 = readpsfcat_simple_wcs(catname, radegs, decdegs, rmag);
        if (nobj2 != nobjt) {
            printf(" error reading catalogue \n");
            printf(" nobjt = %d, nobj2 = %d \n", nobjt, nobj2);
            fflush(stdout);
            exit(2);
        }

        if (nobjt > 0) {
            /*
             * check coords are in sensible range 
             */
            for(i = 0; i < nobjt; i++) {
                if (radegs[i] < 0. || radegs[i] > 360. || decdegs[i] < -90.
                    || decdegs[i] > 90.) {
                  fflush(stdout);
                    fprintf(stderr,
                            " these don't look like world coordinates \n");
                    fprintf(stderr, " %lf %lf \n", radegs[i], decdegs[i]);
                    exit(2);
                }
            }
        }
    }

    if (nobjt > 0) {
        printf(" %d objects read from PSF catalogue %s \n", nobjt, catname);
    } else {
        printf(" no objects read from PSF catalogue %s \n", catname);
        exit(2);
    }

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
        sorder = abs(sorder);
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
        fprintf(stderr, " total number of coefficients is too high \n");
        exit(2);
    }

    /*
     * memory for fitting 
     */

    nn = 0;
    chip = (int *) calloc(nobjt, sizeof(int));
    wfits = (double *) calloc(nobjt, sizeof(double));
    ow = (double *) calloc(nobjt, sizeof(double));
    xfits = (double *) calloc(nobjt, sizeof(double));
    yfits = (double *) calloc(nobjt, sizeof(double));
    xshiftfits = (double *) calloc(nobjt, sizeof(double));
    yshiftfits = (double *) calloc(nobjt, sizeof(double));
    if (chip == NULL || wfits == NULL || xfits == NULL || yfits == NULL
        || xshiftfits == NULL || yshiftfits == NULL) {
        fflush(stdout);
        fprintf(stderr, " memory allocation error, fits arrays \n");
        fflush(stderr);
        exit(2);
    }

    maxncoeffs = sncoeffs;
    if (maxncoeffs < 3)
        maxncoeffs = 3;

    u = (double **) calloc((1 + maxncoeffs), sizeof(double *));
    memory += (float) (1 + maxncoeffs) * sizeof(double);
    for(i = 0; i <= maxncoeffs; i++) {
        memory += (float) (1 + maxncoeffs) * sizeof(double);
        u[i] = (double *) calloc((1 + maxncoeffs), sizeof(double));
    }

    v = (double **) calloc((1 + maxncoeffs), sizeof(double *));
    memory += (float) (1 + maxncoeffs) * sizeof(double);
    for(i = 0; i <= maxncoeffs; i++) {
        memory += (float) (1 + maxncoeffs) * sizeof(double);
        v[i] = (double *) calloc((1 + maxncoeffs), sizeof(double));
    }

    avals = (double *) calloc(1 + maxncoeffs, sizeof(double));
    w = (double *) calloc(1 + maxncoeffs, sizeof(double));
    memory += 2. * (1 + maxncoeffs) * sizeof(double);

    xcoeffs = (double *) calloc((maxncoeffs), sizeof(double));
    memory += (float) (maxncoeffs) * sizeof(double);
    if (xcoeffs == NULL) {
        printf("Memory allocation error for coeffs \n");
        exit(2);
    }
    ycoeffs = (double *) calloc((maxncoeffs), sizeof(double));
    memory += (float) (maxncoeffs) * sizeof(double);
    if (ycoeffs == NULL) {
        printf("Memory allocation error for coeffs \n");
        exit(2);
    }

    diffmax = -1.e10;
    diffmin = 1.e10;


    /*
     * got through each input image and accumulate the star postage stamps etc 
     */

    fflush(stdout);

    for(image = 0; image < nchip; image++) {

        /*
         * strip out the chip number and remember it 
         */
        bzero(firstname, 500);
        strcpy(firstname, image_file[image]);
	// work out the chip ID 
	if (nchip>1)
	  {
	    // if multiple chips, attempt to read the chip number from the filename
	    item = strtok(firstname, delims);
	    item = strtok(NULL, delims2);
	    i = atoi(item) - 1;
	  }
	else
	  {
	    // only one chip expected, so set chip ID to 0
	    i = 0;
	  }
       //ychip = i / nxchip;
       // xchip = i - ychip * nxchip;
        printf(" image %s chip %d %d \n", image_file[image], xchip[i] + 1,
               ychip[i] + 1);
        if (i < 0 || i >= nchip_expected) {
          fflush(stdout);
          fprintf(stderr," error reading chip number for image %d \n",image+1);
          fprintf(stderr," chip number %d out of %d \n",i+1,npsfchip);
          exit(2);
        }

        /*
         * strip off .fits extension on first argument if supplied 
         */
        len = strlen(image_file[image]);
        if (len <= 0) {
          fflush(stdout);
            fprintf(stderr, " error reading arguments from command line \n");
            exit(2);
        }
        if (strncmp(image_file[image] + len - 5, ".fits", 5) == 0
            || strncmp(image_file[image] + len - 5, ".FITS", 5) == 0) {
            len = len - 5;
        }

        bzero(firstname, 500);
        strncpy(firstname, image_file[image], len);

        bzero(imageweightname, 500);
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

        bzero(headerfile, 500);

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
         * routines instead but leave this as it is just now)
         */
        if (getimagesize
            (imagename, satname, dim, &gain, &satlev, &arcperpix, &angle,
             &badccd, &imageid, (float) minexposuretime)) {
            fflush(stdout);
            fprintf(stderr, " error checking input image %s \n", imagename);
            exit(2);
        }

	if (imageid <= 0)
      {
        fflush(stdout);
        fprintf(stderr," ERROR keyword IMAGEID not found with valid value in FITS header\n");
        exit(EXIT_FAILURE);
      }
    // imageid = image+1;

    /*
        if (imageid != (1 + image)) {
            fflush(stdout);
            fprintf(stderr, "error in identifying chip number, %s \n",
                    imagename);
            fflush(stderr);
            exit(2);
        }
    */

    ichip = imageid - 1;
    
    if (ichip<0 || ichip >= nchip)
      {
        fflush(stdout);
        fprintf(stderr," chip number %d out of range %d %d \n",ichip,0,nchip);
        exit(EXIT_FAILURE);
      }

    if (badccd == 1) {
      printf(" %s flagged as bad, skipping this image \n", imagename);
      continue;
    }

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
                exit(2);
            }
            fclose(filep);
        }

        imagesize = dim[0] * dim[1];

	if (VERBOSE == 1)
	  printf(" image dimensions = %d x %d \n", dim[0], dim[1]);

        //      printf (" CCD gain = %6.1f e/ADU \n",gain);

	if (VERBOSE == 1)
	  printf(" CCD saturation level = %8.0f ADU \n", satlev);

        /*
         * set maximum pixel value as 50 percent of CCD saturation level
         * in case of CCD non-linearity 
         */
        satlimit = 0.5 * satlev;

        /*
         * set maximum intensity used for object detection 
         */
        fmax_intensity = satlimit;

	if (VERBOSE == 1)
	  printf(" trying to allocate %d MB for images and arrays\n",
               (int) (15. * imagesize * 4. / 1024. / 1024.));

        region = (int *) calloc(imagesize, sizeof(int));
        if (region == NULL) {
            fprintf(stderr, " error allocating memory for region mask %d \n",
                    i);
            exit(1);
        }
        memloop += (float) imagesize *sizeof(int);

        /*
         * allocate memory for input image 
         */
        apix = (float *) calloc(imagesize, sizeof(float));
        if (apix == NULL) {
            fprintf(stderr, " error allocating memory for image \n");
            exit(1);
        }
        memloop += (float) imagesize *sizeof(float);
        badpix = (float *) calloc(imagesize, sizeof(float));
        if (badpix == NULL) {
            fprintf(stderr, " error allocating memory for image \n");
            exit(1);
        }
        memloop += (float) imagesize *sizeof(float);
	if (SUBTRACT_MEDIAN == 1)
	  sortarray = ml_calloc(imagesize, sizeof(float), &memloop,
				" error allocating memory for sortarray\n");
	if (FILTER_WEIGHTS == 1)
	  {
	    opix = ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for image \n");
	    changed = ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for image \n");
	    weightfilterarray = ml_calloc(1+2*fheight, sizeof(float), &memloop,
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
        getdata_badpix(imagename, imageweightname, dim, apix, badpix, &noise);

	if (FILTER_WEIGHTS == 1)
	  {
	    nchanged = weightfilter(badpix, opix, changed, weightfilterarray, dim, fheight);
	    dilatemask(opix,dim);
	  }
	else
	  {
	    nchanged = 0;
	    opix = badpix;
	  }

	// dilate pixel mask by one pixel in all directions
	dilatemask(badpix,dim);

	if (VERBOSE == 1)
	  printf(" estimated noise in image = %f ADU \n", noise);

        if (noise <= 0.) {
            fprintf(stderr, " Erroneous noise value \n");
            exit(1);
        }

        /*
         * subtract assumed constant median background 
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
        fintensity_limit = 3. * noise;

        /*
         * maximum object detection pixel value has already been set to half of
         * the sat level.  Now allocate arrays large enough to hold the 
         * digitised pixel values.  This must match with assumed array sizes
         * inside varylpthresholdf 
         */
        if (noise > 0.) {
            scale = 2. / noise;
            imax_intensity = (int) (fmax_intensity * scale);
        } else {
            printf(" noise value < 0 \n");
            exit(2);
        }

        nobjects = varylpthresholdf2(apix, region, noise, dim[0], dim[1],
                                     lower_area_limit, lower_merging_limit,
                                     arealimit, fintensity_limit,
                                     fmax_intensity, imax_intensity);
	if (VERBOSE == 1)
	  printf(" %d objects identified within image \n", nobjects);

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
            if (wcs_to_raw(wcs_raw, wcspos, rawpos) == RETURN_OK) {
                if (rawpos[0] > 0. && rawpos[0] < (double) rawfield->width) {
                  if (rawpos[1] > 0. && rawpos[1] < (double) rawfield->height) {
                    // convert back to wcs to check the lookup was correct
                    // (can fail for multivalued astrometric solutions)
                    raw_to_wcs(wcs_raw, rawpos, wtest);
                    if (*wtest == WCS_NOCOORD) {
                      fflush(stdout);
                      fprintf(stderr," error in inverse position lookup \n");
                      fprintf(stderr, " image pixel position %lf %lf \n", rawpos[0],rawpos[1]);
                      exit(EXIT_FAILURE);
                    }
                    wcstest =
                      pow(((wtest[0]-wcspos[0])/cos(wcspos[1]*pi/180.)),2) +
                      pow((wtest[1]-wcspos[1]),2);
                    if (wcstest < wcstesttol)
                      {
                        dobjx[ii] = rawpos[0];
                        dobjy[ii] = rawpos[1];
                        offscale[ii] = 0;
                      }
                    else
                      {
                        // give a warning but not an error if this test is failed
                        printf(" WARNING input object rejected by WCS test %lf %lf -> %lf %lf -> %lf %lf\n",
                               wcspos[0],wcspos[1],
                               rawpos[0],rawpos[1],
                               wtest[0],wtest[1]);
                        fflush(stdout);
                      }
                  }
                }
            }
            /*
            if (offscale[ii]==0)
              printf(" pos check %lf %lf %lf %lf %lf \n",radegs[ii],decdegs[ii],dobjx[ii],dobjy[ii],rmag[ii]);
            */
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
                //printf(" %lf %lf %f %f \n",radegs[i],decdegs[i],objx[nobj],objy[nobj]);
                nobj++;
            } else {
                // printf(" %lf %lf off-image \n",radegs[i],decdegs[i]);
            }
        }

        dbadpix = (double *) calloc((pwidth * pheight), sizeof(double));
        memloop += (float) (pwidth * pheight) * sizeof(double);
        /*
         * printf(" Allocated ppix for image %d, %d\n",i,pwidth*pheight);
         */
        if (dbadpix == NULL) {
            printf("Memory allocation error for sub-image \n");
            exit(2);
        }


        /*
         * set limits on portion of image that a PSF must occupy (excludes
         * a narrow band around the edges) 
         */
        xmin = pwidth / 2;
        ymin = pheight / 2;
        xmax = dim[0] - pwidth / 2;
        ymax = dim[1] - pheight / 2;

        /*
         * read, normalise, fft and store all the psfs 
         */

	if (VERBOSE == 1)
	  printf(" total memory allocated %d MB \n",
               (int) ((memory + memloop) / 1024. / 1024.));

        /*
         * loop through each star from the catalogue, find which psf it belongs to:  
         * then extract the data 
         */

        nstars = 0;
        num_snratio = 0;
        num_saturated = 0;
        num_off = 0;
        num_wrong_mag = 0;
        num_badpix = 0;
        ncross = 0;


        /*
         * start looping through list of PSF stars 
         */

        for(i = 0; i < nobj; i++) {
            x = (int) (objx[i] + 0.5) - 1;
            y = (int) (objy[i] + 0.5) - 1;

            if (x >= xmin && x < xmax && y >= ymin && y < ymax) {
                if (mag[i] >= brightmag && mag[i] <= faintmag) {
                    ximageoffset = 0.;
                    yimageoffset = 0.;
                    maxvalue = 0.;
                    maxlevel = satlev / 2.;
                    if (USE_SWARP == 0)
                      {
                        ngoodpix = extractpostagestamp
                          (apix, badpix, dim, dpix, dbadpix, temp, badtemp,
                           objx[i], objy[i], poserror,
                           ximageoffset, yimageoffset, pwidth, pheight,
                           noise, region, maxlevel);
                      }
                    else if (CORRECT_DISTORTION == 1) 
                      {
                        ngoodpix = swarpextract(wcs_raw, apix, badpix, dim, dpix,
						dbadpix, objx[i], objy[i], poserror,
						pwidth, pheight, ikernel, noise,
						region, maxlevel, scalefactor);
		      } 
		    else 
		      {
                flagbad = 1;
                ngoodpix = tbswarpextract_nodistortion
                  (wcs_raw, apix, badpix, opix, dim, dpix,
                   dbadpix, objx[i], objy[i],
                   poserror, pwidth, pheight,
                   noise, region, scalefactor, maxlevel,
                   rawpos, wcspos, wcsneg, wcscentre, flagbad);
		      }

                    fflush(stdout);

                    if (ngoodpix > 0) {
                        psfnorm = noisevalue = maxvalue = 0.;
                        for(iy = 0; iy < pheight; iy++) {
                            yy = iy > pheight / 2 ? iy - pheight : iy;
                            for(ix = 0; ix < pwidth; ix++) {
                                xx = ix > pwidth / 2 ? ix - pwidth : ix;
                                if (xx * xx + yy * yy <= poserror*poserror) {
                                    pixel = ix + iy * pwidth;
                                    psfnorm += dpix[pixel];
                                    noisevalue += noise * noise;
                                    if (dpix[pixel] > maxvalue)
                                        maxvalue = dpix[pixel];
                                }
                            }
                        }
                        if (noisevalue <= 0.) {
                            fflush(stdout);
                            fprintf(stderr, " negative noise value \n");
                            exit(2);
                        }
                        noisevalue = sqrt(noisevalue);
                        /*
                         * printf(" maxvalues %lf %lf %f %f %d %f %f \n",
                         * rdegs[i],ddegs[i],objx[i],objy[i],chipnumber+1,maxvalue/psfnorm,psfnorm/noisevalue);
                         */
                        if (psfnorm > (double) (snratio * noisevalue)) {
                            // check star isn't saturated nor a likely cosmic ray hit
                          if (TEST_MAXIMUM != 1 || (maxvalue <= satlimit && maxvalue < 0.15 * psfnorm)) {
                                /*
                                 * this star meets the selection criteria 
                                 */
                                nstars++;

                                /*
                                 * S/N 
                                 */
                                sn[i] = psfnorm / noisevalue;

                                /*
                                 * reconstruct the PSF here 
                                 */
                                /*
                                 * define global x,y position in range -1 < val < 1 
                                 */
                                if (imageid <= 0 || npsfchip <= 0) {
                                    fprintf(stderr,
                                            " error in identifying chip for image %d \n",
                                            image + 1);
                                    exit(2);
                                }
                                //ichip = imageid - 1;
                                //ychip = ichip / nxchip;
                                //xchip = ichip - ychip * nxchip;
                                xval =
                                    (objx[i] +
                                     xchip[ichip] * xchipsampling) / hxsize - 1.;
                                yval = objy[i] + ychip[ichip] * ychipsampling;
                                /*
                                 * assume megacam type of geometry 
                                 */
                                if (ychip[ichip] > 0)
                                    yval += big_gap;
                                if (ychip[ichip] > 2)
                                    yval += big_gap;
                                yval = yval / hysize - 1.;
                                for(pixel = 0; pixel < (pwidth * pheight);
                                    pixel++) {
                                    globalreconstruct(xval, yval, order,
                                                      crossterm,
                                                      chipvariation,
                                                      chiporder, ichip, nchip,
                                                      acoeffs[pixel],
                                                      &bpix[pixel]);
                                }
                                /*
                                 * find the centroid position 
                                 */
                                psfmoments(bpix, pheight, pwidth, psfe, centroid, moments);
                                /*
                                 * FFT the PSF 
                                 */
                                fftw_execute(psfplan);
                                /*
                                 * FFT the star 
                                 */
                                fftw_execute(starplan);
                                /*
                                 * measure the position shift 
                                 */
                                xcorr_measure(padinv, Bpix, Dpix,
                                              pheight, pwidth, padheight,
                                              padwidth, C, c, Cpad, shift);
                                /*
                                 * if the shift is a sensible amount, shift the PSF by the cross-correlation amount 
                                 */
                                //                      if (fabs(shift[0])<poserror && fabs(shift[1])<poserror)
                                //{
                                for(iy = 0; iy < pheight; iy++) {
                                    for(ix = 0; ix < halfpwidth; ix++) {
                                        pixel = ix + iy * halfpwidth;
                                        shiftpixft[pixel] = Bpix[pixel];
                                    }
                                }
                                /*
                                 * apply the shift 
                                 */
                                shiftft(pheight, pwidth, shiftpixft, shift);
                                /*
                                 * transform back the shifted PSF 
                                 */
                                fftw_execute(iqplan);
                                //}

                                cchwidth = poserror;
                                /*
                                 * calculate cross-correlation between psf and star only in central region 
                                 */
                                psfauto = 0.;
                                dauto = 0.;
                                cross = 0.;
                                norm = 0.;
                                for(iy = 0; iy < pheight; iy++) {
                                  yy = iy>pheight/2 ? iy-pheight : iy;
                                  for(ix = 0; ix < pwidth; ix++) {
                                    xx = ix>pwidth/2 ? ix-pwidth : ix;
                                    pixel = ix + iy * pwidth;
                                    if (xx*xx+yy*yy <= cchwidth*cchwidth)
                                      {
                                            psfauto +=
                                                shiftpix[pixel] *
                                                shiftpix[pixel];
                                            dauto +=
                                                dpix[pixel] * dpix[pixel];
                                      }
                                    }
                                }

                                pnorm = 0.;
                                dnorm = 0.;
                                for(iy = 0; iy < pheight; iy++) {
                                  yy = iy>pheight/2 ? iy-pheight : iy;
                                  for(ix = 0; ix < pwidth; ix++) {
                                    xx = ix>pwidth/2 ? ix-pwidth : ix;
                                    pixel = ix + iy * pwidth;
                                    if (xx*xx+yy*yy <= cchwidth*cchwidth)
                                      {
                                        cross +=
                                          dpix[pixel] *
                                          shiftpix[pixel] / sqrt(dauto *
                                                                 psfauto);
                                        norm +=
                                          dpix[pixel] *
                                          shiftpix[pixel] / psfauto;
                                        dnorm += dpix[pixel];
                                        pnorm += shiftpix[pixel];
                                      }
                                  }
                                }
                                sx = shift[0];
                                sy = shift[1];
                                /*
                                 * correct the shift for the extraction pixelisation 
                                 */
                                xs = objx[i] - (int) (objx[i] + 0.5);
                                ys = objy[i] - (int) (objy[i] + 0.5);
                                shift[0] -= xs;
                                shift[1] -= ys;
                                /*
                                 * correct for centroid shift 
                                 */
                                shift[0] += centroid[0];
                                shift[1] += centroid[1];
                                /*
                                 * record the corrected position shift for this star location 
                                 */
                                printf (" %f %f %f %f %f %f %f %f %f \n",
                                        rdegs[i],ddegs[i],objx[i],objy[i],shift[0],shift[1],
                                        centroid[0],centroid[1],cross);

                                /*
                                 * load shift values for good stars into arrays for fitting 
                                 */
                                if (cross > 0.95 && dnorm > 0. && pnorm > 0.) {
                                    chip[nn] = ichip;
                                    // modified S/N weighting
                                    wfits[nn] =
                                        1. / (1. +
                                              50. * 50. / (sn[i] * sn[i]));
                                    xfits[nn] = xval;
                                    yfits[nn] = yval;
                                    xshiftfits[nn] = shift[0];
                                    yshiftfits[nn] = shift[1];
                                    rx[nn] = rdegs[i];
                                    ry[nn] = ddegs[i];
                                    // printf(" %d %d %d %lf %lf \n",nn,i,chip[nn],rdegs[i],ddegs[i]);

				/* measure moments of shifted PSF and data and output to file */
                                    psfmoments(shiftpix, pheight, pwidth, psfe, centroid, moments);
                                    fprintf(ellfile,"%d %d %f %f %f %f %f %f %f %f %f %f ", 
                                            i+1, ichip+1, rdegs[i], ddegs[i], mag[i], 
                                            centroid[0], centroid[1], psfe[0], psfe[1], 
                                            moments[0], moments[1], moments[2]);
                                    psfmoments(dpix, pheight, pwidth, psfe, centroid, moments);
                                    fprintf(ellfile,"%f %f %f %f %f %f %f %f\n ", 
                                            centroid[0], centroid[1], psfe[0], psfe[1], 
                                            moments[0], moments[1], moments[2], wfits[nn]);
                                    
                                    /*
                                     * count total number of useful values for this pixel 
                                     */
                                    nn++;
                                    num_obj_used[image]++;

                                    // accumulate PSF-star differences
                                    for(iy = 0; iy < pheight; iy++) {
                                        yy = iy + pheight / 2;
                                        if (yy >= pheight)
                                            yy -= pheight;
                                        for(ix = 0; ix < pwidth; ix++) {
                                            xx = ix + pwidth / 2;
                                            if (xx >= pwidth)
                                                xx -= pwidth;
                                            pixel = ix + iy * pwidth;
                                            j = xx + yy * pwidth;
                                            diff[ichip][j] +=
                                                shiftpix[pixel] / pnorm -
                                                dpix[pixel] / dnorm;
                                        }
                                    }
                                    // printf("image %d %lf %lf %f %f %f %f \n",image,rdegs[i],ddegs[i],mag[i],sx,sy,shiftpix[0]/pnorm-dpix[0]/dnorm);
                                } else {
                                    ncross++;
                                }

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

        if (ichip>=0 && ichip<nchip) {
          for (iy = 0; iy < pheight; iy++) {
            for(ix = 0; ix < pwidth; ix++) {
              pixel = ix + iy * pwidth;
              if (num_obj_used[image] > 0) {
                diff[ichip][pixel] =
                  diff[ichip][pixel] / num_obj_used[image];
              } else {
                diff[ichip][pixel] = 0.;
              }
            }
          }
        }
        
        printf(" excluded stars: \n");
        printf(" %d stars outside image area \n", num_off);
        printf(" %d stars outside mag range \n", num_wrong_mag);
        printf(" %d stars with no good data \n", num_badpix);
        if (snratio > 0.)
            printf(" %d stars with peak S/N ratio < %8.1f \n", num_snratio,
                   snratio);

        if (TEST_MAXIMUM == 1 || num_saturated>0)
          printf(" %d stars with peak > %8.0f ADU \n", num_saturated, satlimit);
        else
          printf(" no filtering on star peak brightness applied\n");

        printf(" %d stars rejected after cross-correlation \n", ncross);

        printf(" %d stars matching selection criteria \n", nstars);

        fflush(stdout);

        /*
         * free memory ready for next image 
         */
        free(region);
        free(apix);
        free(badpix);
        if (FILTER_WEIGHTS == 1)
          {
            free(opix);
            free(changed);
            free(weightfilterarray);
          }
        if (SUBTRACT_MEDIAN == 1)
          free(sortarray);

        free(dbadpix);
        free_cat(&rawcat, 1);
        free(rawfield);
        
        //printf(" memory freed ready for next image\n");

        memloop = 0.;

        /*
         * next image 
         */

        fflush(stdout);
        
    }

    if (nn <= 0) {
        fflush(stdout);
        fprintf(stderr, " no stars selected for fitting \n");
        exit(2);
    }

    /*
     * write diff image to fits file 
     */
    int anaxis = 3;
    long fpixel[3] = { 1, 1, 1 };
    long anaxes[3];
    anaxes[0] = pwidth;
    anaxes[1] = pheight;
    anaxes[2] = nchip;
    int bitpix = -64;
    int testsize = pwidth * pheight;
    fitsfile *tfptr;

    // if this file already exists make a new name
    while (access(stampname, F_OK) == 0) {
        remove(stampname);
    }
    printf(" writing PSF-star difference image to file %s \n", stampname);

    // create the file and write the postage stamps array 
    fits_create_file(&tfptr, stampname, &status);
    fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
    if (status) {
        printf(" error creating difference image  %s \n", stampname);
        fits_report_error(stderr, status);
    }
    for(i = 0; i < nchip; i++) {
        fpixel[2] = i + 1;
        if (fits_write_pix
            (tfptr, TDOUBLE, fpixel, testsize, diff[i], &status)) {
            printf(" error writing difference image pixel data %s %d \n",
                   stampname, i);
            fits_report_error(stderr, status);
        }
    }
    fits_close_file(tfptr, &status);
    if (status) {
        printf(" error closing difference image file %s \n", stampname);
        fits_report_error(stderr, status);
    }


    /*
     * store original weights 
     */
    for(i = 0; i < nn; i++)
        ow[i] = wfits[i];

    xrms = yrms = 1.;

    for(iter = 0; iter < 3; iter++) {
        printf(" star shift fitting, iteration %d \n", iter + 1);
        /*
         * fit polynomials to x & y shifts across mosaic (use same routine as for PSF pixel values 
         */
        globalsvdfit(xfits, yfits, xshiftfits, wfits, chip, nn,
                     nchip, sorder, scrossterm, schipvariation, schiporder,
                     sncoeffs, avals, u, v, w);
        /*
         * put fit coefficients into global pixel array 
         */
        for(j = 0; j < sncoeffs; j++) {
            xcoeffs[j] = avals[j + 1];  /* shift by 1 from NR routine */
        }
        globalsvdfit(xfits, yfits, yshiftfits, wfits, chip, nn,
                     nchip, sorder, scrossterm, schipvariation, schiporder,
                     sncoeffs, avals, u, v, w);
        /*
         * put fit coefficients into global pixel array 
         */
        for(j = 0; j < sncoeffs; j++) {
            ycoeffs[j] = avals[j + 1];  /* shift by 1 from NR routine */
        }

        sum = sumx = sumy = sumxsq = sumysq = 0.;

        /*
         * test if any stars look way off, set their weights to zero and refit 
         */
        for(i = 0; i < nn; i++) {
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chip[i], nchip,
                              xcoeffs, &xs);
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chip[i], nchip,
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
                           rx[i], ry[i], chip[i], wfits[i], xshiftfits[i], xs,
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

    fprintf(shiftsfile,
     "# good stars ra, dec, chip, global x, global y, fitted xshift, fitted yshift, measured xshift, measured yshift\n");
    fprintf(shiftsfile, "# hxsize %f\n",hxsize);
    fprintf(shiftsfile, "# hysize %f\n",hysize);
    fprintf(shiftsfile, "# nchip %d\n",nchip);

    for(i = 0; i < nn; i++) {
        if (wfits[i] > 0.) {
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chip[i], nchip,
                              xcoeffs, &xs);
            globalreconstruct(xfits[i], yfits[i], sorder, scrossterm,
                              schipvariation, schiporder, chip[i], nchip,
                              ycoeffs, &ys);
            // print to file with conversion of pixel coordinates back to pixel values
            // in the global mosaic frame
            fprintf(shiftsfile," %lf %lf %d %f %f %f %f %f %f \n", rx[i], ry[i], 1+chip[i],
                    (1.+xfits[i])*hxsize, (1.+yfits[i])*hysize,
                    xs, ys, xshiftfits[i], yshiftfits[i]);
        }
    }

    updateglobalcoeffs(coeffname, sncoeffs, sorder, schipvariation,
                       schiporder, xcoeffs, ycoeffs);

    printf(" coefficients FITS file %s updated \n", coeffname);

    fclose(ellfile);
    fclose(shiftsfile);
            
    return 0;
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
        exit(2);
    }

/* read dimensions */
    fits_get_img_dim(afptr, &anaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(2);
    }

    if (anaxis != 2) {
        printf
            ("Error: images with other than 2 dimensions are not supported\n");
        exit(2);
    }

    if (dim[0] != anaxes[0]) {
        fprintf(stderr, " error reading image dimensions \n");
        exit(1);
    }
    if (dim[1] != anaxes[1]) {
        fprintf(stderr, " error reading image dimensions \n");
        exit(1);
    }

    size = dim[0] * dim[1];
    bsize = dim[0] * dim[1];

/* read input data into image array */

    if (fits_read_pix(afptr, TFLOAT, fpixel, size, NULL, apix, NULL, &status)) {
        printf(" error reading pixel data \n");
        exit(2);
    }

/* close main image file */

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(2);
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
            exit(2);
        }

        if (bnaxis != 2) {
            printf
                ("Error: images with other than 2 dimensions are not supported\n");
            exit(2);
        }

        if (dim[0] != bnaxes[0]) {
            fprintf(stderr, " error reading bad pixel image dimensions \n");
            exit(1);
        }
        if (dim[1] != bnaxes[1]) {
            fprintf(stderr, " error reading bad pixel image dimensions \n");
            exit(1);
        }

/* read input data into image array */

        if (fits_read_pix(bfptr, TFLOAT, fpixel, bsize, NULL, badpix,
                          NULL, &status)) {
            printf(" error reading bad pixel mask \n");
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }

        fits_close_file(bfptr, &status);

        if (status) {
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }

    } else {
        if (WEIGHTS_REQUIRED > 0) {
            fprintf(stderr, " bad pixel/weight mask not found \n");
            printf(" bad pixel/weight mask not found %s \n", badpixelname);
            exit(2);
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

    if (medbin > 0) {
        qsort(noiseval, medbin, sizeof(float), compare);
        *mednoiseval = noiseval[(medbin / 2)];
    } else {
        *mednoiseval = 0.;
    }

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
