#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

#define VERBOSE 0

int
getimagesize(char *imagename, char *satlevkey, int *dim, float *gain,
             float *satlev, float *arcperpix, float *angle, int *badccd,
             int *imageid, float minexposuretime)
{
/* opens main FITS image and gets basic info */

    fitsfile *afptr;            /* FITS file pointer */
    int status = 0;             /* CFITSIO status value must be initialized to zero */
    int anaxis;
    long anaxes[2] = { 1, 1 };
    int size;
    float exptime;
    char *anglekey = "ORIENTAT";
    char *gainkey = "GAIN";
    char *satlevkey2 = "SATLEVEL";
    char *arcperpixkey = "PIXSCAL1";
    char comment[FLEN_COMMENT];

/* open input images */
    fits_open_file(&afptr, imagename, READONLY, &status);

/* read dimensions */
    fits_get_img_dim(afptr, &anaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis != 2) {
        fflush(stdout);
        fprintf(stderr,
                "Error: images with other than 2 dimensions are not supported\n");
        exit(2);
    }

    dim[0] = (int) anaxes[0];
    dim[1] = (int) anaxes[1];

    size = dim[0] * dim[1];


    /*
     * get the PA value from the fits header 
     */
    ffgky(afptr, TFLOAT, anglekey, angle, comment, &status);
    if (status) {
        angle[0] = 0.;
        status = 0;
    } else {
        if (VERBOSE == 1)
            printf(" image position angle = %f \n", angle[0]);
    }

    /*
     * get the gain value from the fits header 
     */
    ffgky(afptr, TFLOAT, gainkey, gain, comment, &status);
    if (status) {
        // printf (" GAIN keyword not found, setting gain=1 \n");
        gain[0] = 1.;
        status = 0;
    }

    /*
     * try two alternate keyword names for the saturation level 
     */
    ffgky(afptr, TFLOAT, satlevkey, satlev, comment, &status);
    if (status) {
        status = 0;
        ffgky(afptr, TFLOAT, satlevkey2, satlev, comment, &status);
        if (status) {
            printf(" saturation keyword %s or %s not found \n", satlevkey,
                   satlevkey2);
            printf
                (" saturation keyword must be set in header for image %s \n",
                 imagename);
            printf(" setting value 60000 asuming raw pixel values \n");
            *satlev = 60000.;
        }
    }

    /*
     * get the arcsec per pixel value from the fits header 
     */
    ffgky(afptr, TFLOAT, arcperpixkey, arcperpix, comment, &status);
    if (status) {
        //      printf (" pixel scale keyword not found, setting pixel scale = 0.185 arcsec \n");
        arcperpix[0] = 0.185;
        status = 0;
    }

    /*
     * get the badccd flag from the fits header 
     */
    ffgky(afptr, TINT, "BADCCD", badccd, comment, &status);
    if (status) {
        *badccd = 0;
        status = 0;
    }

    /*
     * get the chip number from the fits header 
     */
    ffgky(afptr, TINT, "IMAGEID", imageid, comment, &status);
    if (status) {
        *imageid = 0;
        status = 0;
    }

    /*
     * get the exposure time value from the fits header - if this is less than minexposuretime secs
     * set badccd=1 
     */
    ffgky(afptr, TFLOAT, "EXPTIME", &exptime, comment, &status);
    if (status) {
        status = 0;
    } else {
        if (VERBOSE == 1) {
            printf(" exposure time %f \n", exptime);
            fflush(stdout);
        }
        if (exptime < minexposuretime) {
            *badccd = 1;
            if (VERBOSE == 1) {
                printf
                    (" *** exposure time < %f secs: this image flagged as bad \n",
                     minexposuretime);
                fflush(stdout);
            }
        }
    }

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    return status;

}
