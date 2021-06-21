#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
getcoeffs_fits_expand_v5(char *psfname, double **acoeffs, int newpwidth,
                         int newpheight, int pwidth, int pheight, int ncoeffs,
                         int *listnum)
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
    int hdunum, intnull, anynull, hdutype, ncol;
    long frow, felem, numpsfs;
    double inner, outer, radius, radius2, pi;
    char comment[FLEN_COMMENT];

    pi = M_PI;

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


    // test if the NPSFCHIP keyword is set - if it is, then then is a new-style PSF
    // coeffs file, and we need to skip over the table of xchip, ychip values (these
    // have already been read by function getglobalpsfsize)
    int npsfchip = 0;
    if (fits_read_key(afptr, TINT, "NPSFCHIP", &npsfchip, comment, &status)) {
      //fits_report_error(stderr, status);      /* print error message */
      //printf(" keyword NPSFCHIP\n");
        npsfchip = 0;
        status = 0;
    }        
    
    // get the table of star numbers if it exists and if we want it 
    
    // if the table array is allocated, try to move to the HDU, don't worry if can't 
    hdunum = 2;
    if (npsfchip > 0) hdunum++;

    if (listnum != NULL) {
        if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
            status = 0;
            //  fits_report_error(stderr, status); /* print error message */
            //  return(status);
        } else {
            if (hdutype != BINARY_TBL) {
                fprintf(stderr,
                        " PSF file star numbers not a binary table \n");
                exit(2);
            }

            if (fits_get_num_cols(afptr, &ncol, &status)) {
                fits_report_error(stderr, status);      /* print error message */
                return (status);
            }
            if (ncol != 3) {
                fprintf(stderr, " error in table format, no. columns = %d \n",
                        ncol);
                exit(2);
            }

            if (fits_get_num_rows(afptr, &numpsfs, &status)) {
                fits_report_error(stderr, status);      /* print error message */
                return (status);
            }

            frow = 1;
            felem = 1;
            intnull = 0;

            /*
             * read the number of stars column 
             */
            if (fits_read_col
                (afptr, TINT, 3, frow, felem, numpsfs, &intnull, listnum,
                 &anynull, &status)) {
                fits_report_error(stderr, status);      /* print error message */
                return (status);
            }
        }
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
