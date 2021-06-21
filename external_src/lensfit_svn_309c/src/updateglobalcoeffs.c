#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

/* opens a global PSF coefficients file and appends a table of shift coefficients */

int
updateglobalcoeffs(char *psfname, int ncoeffs, int order,
                   int chipvariation, int chiporder, double *xc, double *yc)
{
    fitsfile *afptr;
    int status = 0;
    int i;
    int row;
    int hdunum, hdutype;

    int tfields = 2;            /* number of table columns in table */
    char extname[] = "star shift coefficients"; /* extension name */
    char *ttype[] = { "XC", "YC" };     /* names of columns */
    char *tform[] = { "1D", "1D" };     /* data types of columns */
    char *tunit[] = { "", "" }; /* units */

    fits_open_file(&afptr, psfname, READWRITE, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * write keywords into header 
     */

    if (fits_update_key(afptr, TINT, "SORDER", &order,
                        "order of shifts global fit", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "SCHIPVAR", &chipvariation,
                        "shifts chip-dependency allowed?", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_update_key(afptr, TINT, "SCHORDER", &chiporder,
                        "shifts chip fit order", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    hdunum = 2;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        status = 0;
        fprintf(stderr, " error moving to table of star information \n");
        fits_report_error(stderr, status);      /* print error message */
        exit(2);
    }


    /*
     * create new table 
     */
    fits_create_tbl(afptr, BINARY_TBL, ncoeffs, tfields, ttype, tform, tunit,
                    extname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < ncoeffs; i++) {
        row = i + 1;
        fits_write_col(afptr, TDOUBLE, 1, row, 1, 1, &xc[i], &status);
        fits_write_col(afptr, TDOUBLE, 2, row, 1, 1, &yc[i], &status);
    }
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * close the fits file 
     */

    fits_close_file(afptr, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    return (0);

}
