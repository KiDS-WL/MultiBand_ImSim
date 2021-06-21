#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
getpsfsize_v5(char *psfname, int *pwidth, int *pheight, int *order,
              int *ncoeffs, int *boxsize, int *xboxnum, int *yboxnum)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis, hdunum, hdutype;
    long anaxes[3] = { 1, 1, 1 };
    char comment[64];

    /*
     * open input images 
     */
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

    if (anaxis < 2 || anaxis > 3) {
        printf("Error: expected 2 or 3 dimensions in PSF fits file\n");
        exit(2);
    }

    fits_get_img_size(afptr, anaxis, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis == 2) {
        anaxes[2] = 1;
        order[0] = 0;
    } else {
        if (fits_read_key(afptr, TINT, "ORDER", order, comment, &status)) {
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword ORDER\n");
            return (status);
        }
    }

    /*
     * move to the HDU 
     */
    hdunum = 2;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        status = 0;
        return (status);
    }
    if (hdutype != BINARY_TBL) {
        fprintf(stderr, " PSF file star numbers not a binary table \n");
        exit(2);
    }


    /*
     * get PSF box keywords, but only give a warning if these are not found,
     * for backwards compatibility with older PSF files 
     */
    if (fits_read_key(afptr, TINT, "PSFBOX", boxsize, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword PSFBOX\n");
        *boxsize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NXPSF", xboxnum, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword NXPSF\n");
        *xboxnum = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NYPSF", yboxnum, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword NYPSF\n");
        *yboxnum = 0;
        status = 0;
    }

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    pwidth[0] = (int) anaxes[0];
    pheight[0] = (int) anaxes[1];
    ncoeffs[0] = (int) anaxes[2];

    return (0);

}
