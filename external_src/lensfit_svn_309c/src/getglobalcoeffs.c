#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
getglobalcoeffs(char *psfname, double **acoeffs, int newpwidth,
                int newpheight, int pwidth, int pheight, int ncoeffs,
                int sncoeffs, double *xc, double *yc)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[3] = { 1, 1, 1 }, fpixel[3] = {
    1, 1, 1};
    double *f, factor;
    int psfsize, psfcubesize;
    int i, ii, oldj, newj, halfpheight, halfpwidth;
    int x, y, x1, y1;
    int hdunum, anynull, hdutype, ncol;
    long frow, felem, nrows;
    double inner, outer, radius, radius2, pi, dnull;

    pi = 3.1415927;

    /*
     * read in psf data as a fits file 
     */

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
        fprintf(stderr, " error in dimensions of PSF fits file \n");
        exit(2);
    }

    fits_get_img_size(afptr, anaxis, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis == 2) {
        anaxes[2] = 1;
    }

    if (anaxes[0] != pwidth) {
        printf(" error in PSF image dimensions in getpsf_fits\n");
        printf(" expected %d, got %ld, for dimension 1\n", pwidth, anaxes[0]);
        exit(2);
    }
    if (anaxes[1] != pheight) {
        printf(" error in PSF image dimensions in getpsf_fits\n");
        printf(" expected %d, got %ld, for dimension 2\n", pheight,
               anaxes[1]);
        exit(2);
    }
    if (anaxes[2] != ncoeffs) {
        printf(" error in PSF image dimensions in getpsf_fits\n");
        printf(" expected %d, got %ld, for dimension 3\n", ncoeffs,
               anaxes[2]);
        exit(2);
    }

    psfsize = pwidth * pheight;
    psfcubesize = anaxes[0] * anaxes[1] * anaxes[2];

    f = (double *) calloc(psfcubesize, sizeof(double));
    if (f == NULL) {
        printf(" memory allocation error for array f in getpsf_fits \n");
        printf(" requested %d doubles \n", psfcubesize);
        exit(2);
    }


    /*
     * read input data into image array 
     */

    if (fits_read_pix(afptr, TDOUBLE, fpixel, psfcubesize, NULL, f,
                      NULL, &status)) {
        printf(" error reading pixel data \n");
        exit(2);
    }

    /*
     * skip over the table of star numbers and try to read the table of x,y shifts 
     */

    // try to read the last table in the file
    if (fits_get_num_hdus(afptr, &hdunum, &status)) {
      fflush(stdout);
      fprintf(stderr, " error reading number of HDUs in FITS header \n");
      fits_report_error(stderr, status);      /* print error message */
      exit(2);
    }

    if (hdunum<3)
      {
        fflush(stdout);
        fprintf(stderr, " hdunum = %d, too few HDUs in FITS file\n",hdunum);
        fprintf(stderr, " run program globalshifts first \n");
        exit(2);
      }        
    
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        fflush(stdout);
        fprintf(stderr, " error reading table of position shifts from HDU %d\n",hdunum);
        fits_report_error(stderr, status);      /* print error message */
        exit(2);
    }

    if (hdutype != BINARY_TBL) {
        fprintf(stderr, " PSF file star numbers not a binary table \n");
        exit(2);
    }

    if (fits_get_num_cols(afptr, &ncol, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (ncol != 2) {
        fprintf(stderr, " error in table format, no. columns = %d \n", ncol);
        exit(2);
    }

    if (fits_get_num_rows(afptr, &nrows, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
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


    /*
     * close the fits file 
     */

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * change the PSF size to the new value, with an apodising function applied  
     */

    if (newpheight > pheight) {
        halfpheight = pheight / 2;
    } else {
        halfpheight = newpheight / 2;
    }
    if (newpwidth > pwidth) {
        halfpwidth = pwidth / 2;
    } else {
        halfpwidth = newpwidth / 2;
    }

    /*
     * set range of radius from original outside which PSF will be forced to zero 
     */
    outer = (double) pwidth / 2;

    /*
     * set inner edge of radius at which the tapering starts (inner<outer) 
     */
    inner = 3. * (double) pwidth / 8.;

    for(ii = 0; ii < (newpwidth * newpheight); ii++)
        for(i = 0; i < ncoeffs; i++)
            acoeffs[ii][i] = 0.;

    for(i = 0; i < ncoeffs; i++) {
        for(y = 0; y < halfpheight; y++) {
            for(x = 0; x < halfpwidth; x++) {
                oldj = y * pwidth + x;
                ii = i * psfsize + oldj;
                newj = y * newpwidth + x;
                acoeffs[newj][i] = f[ii];
            }
            for(x = halfpwidth; x < pwidth; x++) {
                oldj = y * pwidth + x;
                ii = i * psfsize + oldj;
                newj = y * newpwidth + x + (newpwidth - pwidth);
                acoeffs[newj][i] = f[ii];
            }
        }
        for(y = halfpheight; y < pheight; y++) {
            for(x = 0; x < halfpwidth; x++) {
                oldj = y * pwidth + x;
                ii = i * psfsize + oldj;
                newj = (y + newpheight - pheight) * newpwidth + x;
                acoeffs[newj][i] = f[ii];
            }
            for(x = halfpwidth; x < pwidth; x++) {
                oldj = y * pwidth + x;
                ii = i * psfsize + oldj;
                newj =
                    (y + newpheight - pheight) * newpwidth + x + (newpwidth -
                                                                  pwidth);
                acoeffs[newj][i] = f[ii];
            }
        }

        for(y = 0; y < newpheight; y++) {
            for(x = 0; x < newpwidth; x++) {
                x1 = newpwidth - x;
                y1 = newpheight - y;
                radius = sqrt((double) (x * x + y * y));
                radius2 = sqrt((double) (x1 * x1 + y * y));
                if (radius2 < radius)
                    radius = radius2;
                radius2 = sqrt((double) (x * x + y1 * y1));
                if (radius2 < radius)
                    radius = radius2;
                radius2 = sqrt((double) (x1 * x1 + y1 * y1));
                if (radius2 < radius)
                    radius = radius2;
                newj = y * newpwidth + x;
                if (radius > inner) {
                    factor =
                        cos((radius - inner) * pi / 2. / (outer - inner));
                    acoeffs[newj][i] = acoeffs[newj][i] * factor * factor;
                }
                if (radius > outer)
                    acoeffs[newj][i] = 0.;
            }
        }
    }

    return (0);

}
