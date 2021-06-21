#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int readoscoeffs(char *psfname, int *pindex, double *soln, 
		 int nfit, int ncoeffs, int sncoeffs,
		 double *xc, double *yc)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    int hdunum, hdutype;
    long anaxes[2], fpixel[2] = { 1, 1 };
    long nrows;
    int bitpix;
    int pindexsize;
    double dnull;
    int anynull, ncol;
    long frow, felem;

    bitpix = 32; // 32-bit integers

    fits_open_file(&afptr, psfname, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * read dimensions 
     */
    fits_get_img_dim(afptr, &anaxis, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis != 2) {
        fprintf(stderr, " error in dimensions of PSF fits file, pindex array\n");
        exit(2);
    }

    fits_get_img_size(afptr, anaxis, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    pindexsize = anaxes[0] * anaxes[1];

    // read pindex as an image
    if (fits_read_pix(afptr, TINT, fpixel, pindexsize, NULL, pindex,
                      NULL, &status)) {
        printf(" error reading pindex data \n");
        exit(EXIT_FAILURE);
    }

    // skip to solution vector and read it as a column of numbers
    hdunum = 2;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        fflush(stdout);
        fprintf(stderr, " error reading table of position shifts \n");
        fits_report_error(stderr, status);      /* print error message */
        fprintf(stderr, " run program globalshifts first \n");
        exit(2);
    }
    if (hdutype != BINARY_TBL) {
        fprintf(stderr, " PSF solution vector is not a binary table \n");
        exit(2);
    }
    if (fits_get_num_cols(afptr, &ncol, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (ncol != 1) {
        fprintf(stderr, 
		" error in table format for input PSF solution vector, no. columns = %d \n", 
		ncol);
        exit(2);
    }
    if (fits_get_num_rows(afptr, &nrows, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (nrows != (long)(nfit*ncoeffs))
      {
	fflush(stdout);
	fprintf(stderr," error in number of coefficients for solution vector\n");
	exit(EXIT_FAILURE);
      }

    frow = 1;
    felem = 1;
    dnull = 0.;

    //    printf(" reading %d\n",(int)nrows);

    if (fits_read_col(afptr, TDOUBLE, 1, frow, felem, nrows, &dnull, soln,
                      &anynull, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * skip over the table of star numbers and try to read the table of x,y shifts 
     */

    hdunum = 4;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        fflush(stdout);
        fprintf(stderr, " error reading table of position shifts \n");
        fits_report_error(stderr, status);      /* print error message */
        exit(2);
    }

    if (hdutype != BINARY_TBL) {
      fprintf(stderr, " PSF file table %d not a binary table \n", hdunum);
        exit(2);
    }

    if (fits_get_num_cols(afptr, &ncol, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_get_num_rows(afptr, &nrows, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    printf(" shift table size %d %d \n",(int)nrows, (int)ncol); fflush(stdout);

    if (ncol != 2) {
        fprintf(stderr, " error in table format, no. columns = %d \n", ncol);
        exit(2);
    }

    if (nrows != sncoeffs) {
        fflush(stdout);
        fprintf(stderr, " error in number of shifts coefficients \n");
        exit(2);
    }

    frow = 1;
    felem = 1;
    dnull = 0.;

    if (fits_read_col(afptr, TDOUBLE, 1, frow, felem, nrows, &dnull, xc,
                      &anynull, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_read_col(afptr, TDOUBLE, 2, frow, felem, nrows, &dnull, yc,
                      &anynull, &status)) {
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
