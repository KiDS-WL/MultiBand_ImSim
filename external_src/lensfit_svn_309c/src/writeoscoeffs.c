#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int writeoscoeffs(char *psfname, int correct_distortion, int *pindex, double *soln, 
		  int nfit, int ncoeffs, int order,
		  int chipvariation, int chiporder, double *scalefactor,
		  int fwidth, int fheight, int oversampling,
		  int nxchip, int nychip,
		  int *numstar, int xchipsampling, int ychipsampling,
		  int big_gap, double hxsize, double hysize, float poserror)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[2], fpixel[2] = { 1, 1 };
    int i, bitpix;
    int psfsize, row;

    // table definition for soln vector
    int sfields = 1;            /* number of table columns in table */
    char sextname[] = "solution vector";        /* extension name */
    char *stype[] = { "SOLN"};        /* names of columns */
    char *sform[] = { "1D"};     /* data types of columns */
    char *sunit[] = { ""}; /* units */

    // table definition for numbers of stars 
    int tfields = 2;            /* number of table columns in table */
    char extname[] = "PSF star numbers";        /* extension name */
    char *ttype[] = { "CHIP", "NSTAR" };        /* names of columns */
    char *tform[] = { "1J", "1J" };     /* data types of columns */
    char *tunit[] = { "", "" }; /* units */

    // write out pindex data as a fits image

    if (fwidth != fheight) {
        fprintf(stderr, " PSF subimages must be square! \n");
        return 1;
    }
    if (fwidth <= 0 || fheight <= 0) {
        fprintf(stderr, " error in subimage dimensions \n");
        return 1;
    }

    bitpix = 32; // 32-bit integers
    anaxis = 2;
    anaxes[0] = fwidth;
    anaxes[1] = fheight;
    psfsize = anaxes[0] * anaxes[1];

    fits_create_file(&afptr, psfname, &status);
    fits_create_img(afptr, bitpix, anaxis, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    // write pindex data into image array
    if (fits_write_pix(afptr, TINT, fpixel, psfsize, pindex, &status)) {
        printf(" error writing psf pindex data \n");
        exit(2);
    }

    // write keywords into header 
    if (fits_update_key(afptr, TINT, "NCOEFFS", &ncoeffs,
                        "number of coefficients", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NFIT", &nfit,
                        "number of fittable pixels", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "CORRECTD", &correct_distortion,
                        "distortion correction", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TDOUBLE, "PIXSCALE", scalefactor,
                        "pixel scalefactor", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "ORDER", &order,
                        "order of global fit", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "CHIPVAR", &chipvariation,
                        "chip-dependency allowed?", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "CHORDER", &chiporder,
                        "chip fit order", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "FWIDTH", &fwidth,
                        "postage stamp width", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "FHEIGHT", &fheight,
                        "postage stamp height", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "OVERSAMP", &oversampling,
                        "oversampling", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NXCHIP", &nxchip,
                        "number of chips along x axis", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NYCHIP", &nychip,
                        "number of chips along y axis", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "XSAMPL", &xchipsampling,
                        "x chip spacing", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "YSAMPL", &ychipsampling,
                        "y chip spacing", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "BIG_GAP", &big_gap,
                        "size of extra gap", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TDOUBLE, "HXSIZE", &hxsize,
                        "x size scaling", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TDOUBLE, "HYSIZE", &hysize,
                        "y size scaling", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TFLOAT, "TOL", &poserror,
                        "tolerance radius", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    // create table for solution vector
    fits_create_tbl(afptr, BINARY_TBL, nfit*ncoeffs, sfields, stype, sform,
                    sunit, sextname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < nfit*ncoeffs; i++) {
        row = i + 1;
        fits_write_col(afptr, TDOUBLE, 1, row, 1, 1, &soln[i], &status);
	if (status) {
	  fflush(stdout);
	  fits_report_error(stderr, status);      /* print error message */
	  fprintf(stderr," error at i=%d \n",i);
	  break;
	}
    }
    printf(" writeoscoeffs: %d x %d = %d solution coefficients\n",
	   nfit,ncoeffs,nfit*ncoeffs);
    fflush(stdout);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    // create table for number of stars
    fits_create_tbl(afptr, BINARY_TBL, nxchip * nychip, tfields, ttype, tform,
                    tunit, extname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < nxchip * nychip; i++) {
        row = i + 1;
        printf(" chip %d numstars %d \n", row, numstar[i]);
        fits_write_col(afptr, TINT, 1, row, 1, 1, &row, &status);
        fits_write_col(afptr, TINT, 2, row, 1, 1, &numstar[i], &status);
    }
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    // close the fits file 
    fits_close_file(afptr, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    return (0);
}
