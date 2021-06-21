/*-----------------------------------------------------------------

createPSFcube to create PSF models at specified locations

Read catalogue of RA, dec.
Reads FITS images for one exposure and uses the astrometry info to calculate
x,y position of objects on each image
Reads PSF model coefficients and creates PSF model at each location
Outputs PSF images to FITS cube and numerical information to ascii file

Revision $Rev: 238 $ 
last changed by $Author: miller $
on $LastChangedDate: 2014-02-22 07:18:57 +0000 (Sat, 22 Feb 2014) $

---------------------------------------------------------------------------*/

#define DATE_STRING "$LastChangedDate: 2014-02-22 07:18:57 +0000 (Sat, 22 Feb 2014) $"

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


int readpsfcatsize(char *);     // read size of catalogue
int readpsfcat_simple_wcs(char *, double *, double *, float *);
int getimagesize(char *, char *, int *, float *, float *, float *, float *,
                 int *, int *, float);
int getglobalcoeffs(char *, double **, int, int, int, int, int, int, double *,
                    double *);
int getglobalpsfsizes(char *, int *, int *, int *, int *, int *, int *, int *, int *,
                      int *, int *, int *, int *, int *, double *, double *, float *,
                      double *, int *, int *, int *, int *);
void globalreconstruct(double, double, int, int, int, int, int, int, double *,
                       double *);
int wcs_to_raw(wcsstruct *, double *, double *);
void load_dfield(fieldstruct *, catstruct *, char *, int, int);
ikernelstruct *init_ikernel(interpenum *, int);
char *create_psf_name(const char *, char *, const char *);
char *create_ospsf_name(const char *, char *, const char *);
void get_poly_order(int, int, int, int, int *, int *);
void psfmoments(double *, int, int, double *, double *, double *);
int open3Dfits(char*, fitsfile**, int);
int write3Dfits(fitsfile*, int, double*);


/* define whether or not to use swarp astrometric corrections
   requires WCS  */
#define USE_SWARP 1

/* Set the instrument here by defining the flag that is checked later
   in the code KIDS, CFHT, SUPRIME or SINGLE */
#define KIDS (1)
/*---------------------------------------------------------------------------*/

/* global variables */
float poserror;
int pwidth, pheight; // sizes of  arrays
int nchip, nstar, order, ncoeffs;
int chipvariation, chiporder;
double xval, yval;
double **shift;
double *sf;
float gain;
float minexposuretime = 0.;

/* main program */
int main(int argc, char *argv[])
{
    int dim[2], badccd, imageid;
    int chipnumber, correctd;
    int i, ii; //, p, q;
    int len, nobjt, nobj2;
    int *numstar;
    int pixel;
    int naxis = 2;
    int narg = 0;
    int image, xchipsampling, ychipsampling, big_gap, nchip_expected,
        nxchip, nychip;
    int xchip, ychip;
    int sorder, schipvariation, schiporder, scrossterm, sncoeffs;

    float memory;

    double **psfevals, centroid[2], **pmoments;

    double xsize, ysize, hxsize, hysize;
    double *radegs, *decdegs;   // ra,dec in input catalogue
    double scalefactor[1];
    double **psf;
    double *xcoeffs, *ycoeffs;

    float psfposerror;
    float *objx, *objy, *mag, *rmag;
    float satlev, arcperpix, angle;
    int crossterm, maxncoeffs;
    int *chipid;

    double **acoeffs;
    double rawpos[2], wcspos[2], *kern[2];

    char *evariable, *headevariable, *psfdir;
    char *prefsname;
    char *catname;
    char **image_file;
    char *imagename, *headerfile;
    char *psfname;
    char *firstname, *shiftsname, *psfcubefile;
    char **argkey, **argval;
    char *pstr, delimiter[5], headdelimiter[5], *satname;

    fitsfile *psfcubepointer; 
    FILE *shiftsfile, *filep;

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


    // variables for output fits file of stars used 
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

    printf("\n version %s ",argv[0]);
    char datedelims[]="()";
    char *dstring = NULL;
    char date_string[200];
    bzero(date_string,200);
    strcpy(date_string, DATE_STRING);
    dstring  = strtok(date_string, datedelims);
    if (dstring != NULL) dstring = strtok(NULL, datedelims);
    if (dstring != NULL) printf(" %s ",dstring);
    printf("\n");

    if (argc != 2) {
        printf
            (" %s <exposurename list of images>\n",
             argv[0]);
        exit(EXIT_FAILURE);
    }

    if (USE_SWARP == 1) {
      //printf(" swarp astrometric corrections will be applied \n");
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
     * initialise memory counters 
     */
    memory = 0.;

    evariable = ml_calloc(PATH_MAX, sizeof(char), &memory, "evariable");
    headevariable = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "headevariable");
    psfdir = ml_calloc(PATH_MAX, sizeof(char), &memory, "psfdir");
    prefsname = ml_calloc(PATH_MAX, sizeof(char), &memory, "prefsname");
    catname = ml_calloc(PATH_MAX, sizeof(char), &memory, "catname");
    imagename = ml_calloc(PATH_MAX, sizeof(char), &memory, "imagename");
    headerfile = ml_calloc(PATH_MAX, sizeof(char), &memory, "headerfile");
    firstname = ml_calloc(PATH_MAX, sizeof(char), &memory, "firstname");
    satname = ml_calloc(PATH_MAX, sizeof(char), &memory, "satname");

    /*
     * create output PSF filenames
     */

    // copy PSF_DIR into filename string
    psfdir = getenv("PSF_DIR");
    strcpy(firstname, argv[1]);

    psfname = create_psf_name(psfdir, firstname, dotdelimiter);

    // then create and append the coefficients file name, depending on 
    // whether we want to use native sampling or oversampling

    /*
    if (strncmp(argv[2],"n",1) == 0)
      {
	// PSF coefficients file at native sampling
	psfname = create_psf_name(psfdir, firstname, dotdelimiter);
      }
    else if (strncmp(argv[2],"o",1) == 0)
      {
	// PSF coefficients file with oversampling
	psfname = create_ospsf_name(psfdir, firstname, dotdelimiter);
      }
    else
      {
	fflush(stdout);
	fprintf(stderr," error in command-line option 2 %s\n",argv[2]);
        fprintf
	  (stderr," %s <exposurename list of images <n(ative)/o(versampled) >\n",
	   argv[0]);
	exit(EXIT_FAILURE);
      }
    */

    if (access(psfname, F_OK) != 0) {
      fflush(stdout);
      fprintf(stderr," Can't read PSF coefficients file %s \n", psfname);
      exit(EXIT_FAILURE);
    } else {
      printf(" using PSF file %s \n",psfname);
    }

    strcpy(firstname, argv[1]);

    shiftsname = ml_calloc(PATH_MAX, sizeof(char), &memory, "shiftsname");
    bzero(shiftsname, PATH_MAX);
    strcpy(shiftsname, firstname);
    strcat(shiftsname, "_coords_shifts.dat");
    
    psfcubefile = ml_calloc(PATH_MAX, sizeof(char), &memory, "psfcubefile");
    bzero(psfcubefile, PATH_MAX);
    strcpy(psfcubefile, firstname);
    strcat(psfcubefile, "_psfcube.fits");

    numstar = calloc(1000, sizeof(int));
    if (numstar == NULL) {
      fflush(stdout);
      fprintf(stderr, " error allocating memory for numstar \n");
      exit(EXIT_FAILURE);
    }

    // get the PSF file parameters from the header
    getglobalpsfsizes(psfname, &correctd, &pwidth, &pheight,
		      &order, &chipvariation, &chiporder,
		      &ncoeffs, numstar, &nxchip, &nychip,
		      &xchipsampling, &ychipsampling, &big_gap,
		      &hxsize, &hysize, &psfposerror, scalefactor,
		      &sorder, &schipvariation, &schiporder,
		      &sncoeffs);

    printf(" creating PSF postage stamps of size %d x %d \n", pwidth, pheight);

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
	strcpy(delimiter,"S");
	strcpy(headdelimiter,"O");
#endif
#ifdef CFHT
	strcpy(delimiter,"C");
	strcpy(headdelimiter,"C");
#endif
#ifdef SINGLE
	strcpy(delimiter,"C");
	strcpy(headdelimiter,"C");
#endif
#ifdef SUPRIME
	strcpy(delimiter,"O");
	strcpy(headdelimiter,"O");
#endif
      }
    else
      {
	strcpy(delimiter,".");
      }

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

    get_poly_order(sorder, schiporder, schipvariation, nchip, &sncoeffs,
                   &scrossterm);

    maxncoeffs = ncoeffs;
    if (sncoeffs > maxncoeffs) maxncoeffs = sncoeffs;
    if (maxncoeffs < 6)
        maxncoeffs = 6;

    xcoeffs = (double *) calloc(maxncoeffs, sizeof(double));
    ycoeffs = (double *) calloc(maxncoeffs, sizeof(double));
    memory += 2. *maxncoeffs * sizeof(double);

    acoeffs = (double **) calloc((pwidth * pheight), sizeof(double *));
    memory += (float) (pwidth * pheight) * sizeof(double *);
    if (acoeffs == NULL) {
        printf("Memory allocation error for sub-image pointers\n");
        exit(EXIT_FAILURE);

    }

    sf = (double*)calloc((pwidth*pheight), sizeof(double));

    for(i = 0; i < (pwidth * pheight); i++) {
        acoeffs[i] = (double *) calloc((ncoeffs), sizeof(double));
        memory += (float) (ncoeffs) * sizeof(double);
        if (acoeffs[i] == NULL) {
            printf("Memory allocation error for sub-image \n");
            exit(EXIT_FAILURE);
        }
    }

    // get the coefficients from the PSF coeffs file
    getglobalcoeffs(psfname, acoeffs, pwidth, pheight,
		    pwidth, pheight, ncoeffs,
		    sncoeffs, xcoeffs, ycoeffs);

    /*
     * get environment variables holding directory names 
     */

    evariable = getenv("DATA_DIR");
    headevariable = getenv("HEAD_DIR");

    /*
     * create the catalogue filename 
     */
    catname = getenv("CATALOGUE");
    if (catname == NULL) {
        fprintf(stderr,
                " catalogue name not set, use environment variable CATALOGUE \n");
        exit(EXIT_FAILURE);
    }

    if (access(catname, F_OK) != 0) {
        printf(" Can't read catalogue %s \n", catname);
        exit(EXIT_FAILURE);
    }


    /*
     * read catalogue of objects 
     */

    nobjt = readpsfcatsize(catname);

    objx = ml_calloc(nobjt, sizeof(float), &memory, "objx");
    objy = ml_calloc(nobjt, sizeof(float), &memory, "objy");
    chipid = ml_calloc(nobjt, sizeof(int), &memory, "chipid");
    mag = ml_calloc(nobjt, sizeof(float), &memory, "mag");
    radegs = ml_calloc(nobjt, sizeof(double), &memory, "radegs");
    decdegs = ml_calloc(nobjt, sizeof(double), &memory, "decdegs");
    rmag = ml_calloc(nobjt, sizeof(float), &memory, "rmag");
    shift = ml_calloc(nobjt, sizeof(double*), &memory, "shift");
    psf = ml_calloc(nobjt, sizeof(double*), &memory, "psf");    
    pmoments =  ml_calloc(nobjt, sizeof(double*), &memory, "pmoments");  
    psfevals =  ml_calloc(nobjt, sizeof(double*), &memory, "psfevals");  

    for (ii=0; ii<nobjt; ii++)
      {
	psf[ii] = ml_calloc((pwidth*pheight), sizeof(double), &memory, "psf[i]"); 
	shift[ii] = ml_calloc(2, sizeof(double), &memory, "shift[i]"); 
	pmoments[ii] = ml_calloc(3, sizeof(double), &memory, "pmoments[i]"); 
	psfevals[ii] = ml_calloc(2, sizeof(double), &memory, "psfevals[i]"); 
      }	      		      
    
    // read in coordinates assuming in decimal WCS 
    nobj2 = readpsfcat_simple_wcs(catname, radegs, decdegs, rmag);
    if (nobj2 != nobjt) {
      printf(" error reading catalogue \n");
      printf(" nobjt = %d, nobj2 = %d \n", nobjt, nobj2);
      exit(EXIT_FAILURE);
    }

    if (nobjt > 0) {
        printf(" %d objects read from PSF catalogue %s \n", nobjt, catname);
    } else {
      fflush(stdout);
      fprintf(stderr," no objects read from PSF catalogue %s \n", catname);
      exit(EXIT_FAILURE);
    }

    // check coords are in sensible range 
    for(i = 0; i < nobjt; i++) {
      if (radegs[i] < 0. || radegs[i] > 360. || decdegs[i] < -90.
	  || decdegs[i] > 90.) {
	fprintf(stderr,
		" these don't look like world coordinates \n");
	fprintf(stderr, " %lf %lf \n", radegs[i], decdegs[i]);
	exit(EXIT_FAILURE);
      }
    }

    nstar = 0;

    for (image = 0; image < nchip; image++) {

      // strip out the chip number and remember it 
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
            (imagename, satname, dim, &gain, &satlev, &arcperpix, &angle,
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

	if (chipnumber+1 != imageid)
	  {
	    fflush(stdout);
	    fprintf(stderr," error matching up chip number identifiers, %d %d \n",imageid,chipnumber+1);
	    exit(EXIT_FAILURE);
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
                exit(EXIT_FAILURE);
            }
            fclose(filep);
        }

        /*
         * use swarp code to get the full wcs info 
         */
        rawcat = read_cat(imagename);
        QCALLOC(rawfield, fieldstruct, 1);
        load_dfield(rawfield, rawcat, headerfile, 0, 0);        // modified version of load_field
        wcs_raw = copy_wcs(rawfield->wcs);      // get the wcs info

        xsize = (double) dim[0];
        ysize = (double) dim[1];

        /*
         * convert WCS to pixel positions using swarp WCS library 
         */
        for(ii = 0; ii < nobjt; ii++) {
            wcspos[0] = radegs[ii];
            wcspos[1] = decdegs[ii];
            if (wcs_to_raw(wcs_raw, wcspos, rawpos) == RETURN_OK) {
                if (rawpos[0] > 0. && rawpos[0] < (double) rawfield->width) {
                    if (rawpos[1] > 0.
                        && rawpos[1] < (double) rawfield->height) {

		      // object has been found on this chip, store the local position info
		      objx[ii] = rawpos[0];
		      objy[ii] = rawpos[1];
		      chipid[ii] = chipnumber+1;

		      // evaluate its global coordinates
		      // in the PSF coordinate system
		      xval = (objx[ii] + xchip * xchipsampling) / hxsize - 1.;
		      yval = objy[ii] + ychip * ychipsampling;
		      if (ychip > 0)
			yval += big_gap;
		      if (ychip > 2)
			yval += big_gap;
		      yval = yval/ hysize - 1.;

		      // calculate the PSF image given the above info
		      for(pixel = 0; pixel < (pwidth * pheight); pixel++) {
			globalreconstruct(xval, yval, order,
					  crossterm, chipvariation,
					  chiporder, chipnumber, nchip,
					  acoeffs[pixel],
					  &psf[ii][pixel]);
		      }

		      // recover also the fitted position shift at this location
		      globalreconstruct(xval, yval, sorder,
					scrossterm, schipvariation,
					schiporder, chipnumber, nchip,
					xcoeffs, &shift[ii][0]);
		      globalreconstruct(xval, yval, sorder,
					scrossterm, schipvariation,
					schiporder, chipnumber, nchip,
					ycoeffs, &shift[ii][1]);

		      // measure moments and centroid
		      psfmoments(psf[ii], pheight, pwidth,
				 psfevals[ii], centroid, pmoments[ii]);

		      // apply pixelisation correction
		      // (commented out here as this is only relevant to
		      // correct for the extraction of a postage stamp
		      // as a discrete set of non-interpolated pixels)
		      //shift[ii][0] += objx[ii] - (int) (objx[ii] + 0.5) 
		      //shift[ii][1] += objy[ii] - (int) (objy[ii] + 0.5) 

		      // apply centroid shift (this was removed before the shifts
		      // were calculated so must now be added back in, using the
		      // same centroid measure as was initially subtracted
		      shift[ii][0] -= centroid[0];
		      shift[ii][1] -= centroid[1];

		      // count number of objects found
		      nstar++;
                    }
                }
            }
        }

	// free swarp field structure
        free_cat(&rawcat, 1);
        free(rawfield);

	// next image 
    }
    
    printf(" %d stars found on images out of total input number %d \n",nstar,nobjt);

    // write out the PSF cube and the shifts information
    remove(psfcubefile);
    open3Dfits(psfcubefile, &psfcubepointer, nobjt);
    if ((shiftsfile = fopen(shiftsname, "w")) == NULL) {
      fprintf(stderr, " failed to make file %s \n", shiftsname);
    }
    fprintf(shiftsfile,
       "#     RA        dec      chipID     X     Y     x-shift    y-shift   e1    e2   m11   m22   m12 \n");
    for (i=0; i<nobjt; i++)
      {
	write3Dfits(psfcubepointer, i, psf[i]);
	fprintf(shiftsfile," %12.5f %12.5f %3d %8.2f %8.2f %8.5f %8.5f %8.5f %8.5f %8.2f %8.2f %8.2f\n",
		radegs[i],decdegs[i],chipid[i],objx[i],objy[i],shift[i][0],shift[i][1],
		psfevals[i][0], psfevals[i][1], pmoments[i][0], pmoments[i][1], pmoments[i][2]);
      }
    fits_close_file(psfcubepointer, &status);
    fclose(shiftsfile);

    return 0;
}  



char *
create_psf_name(const char *psfdir, char *firstname, const char *dotdelimiter)
{
    char *pstr;
    int len;
    char *coeffname;

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

    /*
    if (access(coeffname, F_OK) == 0) {
        printf(" *** PSF file exists: removing \n");
        remove(coeffname);
    }
    */

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

    /*
    if (access(oscoeffname, F_OK) == 0) {
        printf(" *** oversampled PSF file exists: removing \n");
        remove(oscoeffname);
    }
    */

    return oscoeffname;
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

  printf(" creating 3D FITS image %s size %d %d %d \n",
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
  int x, y, xx, yy, ipix, opix;

  psfsize = pwidth*pheight;

  // swap quadrants
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

  return 0;

}


