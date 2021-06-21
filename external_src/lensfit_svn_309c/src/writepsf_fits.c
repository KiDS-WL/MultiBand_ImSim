#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
writepsf_fits(char *psfname, double **f, int xboxnum, int yboxnum, int pwidth,
              int pheight, int *listnum, float **listx, float **listy,
              int boxsize, int boxsmooth)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[3], fpixel[3] = { 1, 1, 1 };
    int i, bitpix, row, psfnum, psfnumx, psfnumy;
    float *psf, *psf3d;
    int psfsize, psfsize3d, dwidth, dfactor, maxboxno;
    int ix, iy, x, y, j, ij, npixels, maxlistnum;

    /*
     * table of which stars have been used to make each PSF 
     */
    int tfields = 5;            /* number of table columns in table */
    char extname[] = "PSF stars";       /* extension name */
    char *ttype[] = { "NX", "NY", "NSTAR", "STARX", "STARY" };  /* names of columns */
    char *tform[5];             /* data types of columns */
    char *tunit[] = { "", "", "", "pixels", "pixels" }; /* units */

    /*
     * write out psf data as a fits file 
     */

    if (pwidth != pheight) {
        fprintf(stderr, " PSF subimages must be square! \n");
        return 1;
    }

    if (pwidth <= 0 || pheight <= 0 || xboxnum <= 0 || yboxnum <= 0) {
        fprintf(stderr, " error in subimage dimensions or number \n");
        return 1;
    }

    bitpix = -32;
    anaxis = 3;
    anaxes[0] = pwidth;
    anaxes[1] = pheight;
    maxboxno = xboxnum * yboxnum;
    anaxes[2] = maxboxno;

    psfsize = anaxes[0] * anaxes[1];
    psf = (float *) calloc(psfsize, sizeof(float));
    if (psf == NULL) {
        fprintf(stderr,
                " error allocating memory for psf in writepsf_fits \n");
        fprintf(stderr, " memory requested = %d floats \n", psfsize);
        exit(2);
    }


    psfsize3d = maxboxno * psfsize;
    psf3d = (float *) calloc(psfsize3d, sizeof(float));
    if (psf3d == NULL) {
        fprintf(stderr,
                " error allocating memory for 3D psf in writepsf_fits \n");
        fprintf(stderr, " memory requested = %d floats \n", psfsize3d);
        exit(2);
    }

    fits_create_file(&afptr, psfname, &status);
    fits_create_img(afptr, bitpix, anaxis, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * normalise the each psf by the sum of its flux in order to
     * conserve total flux 
     */

    for(psfnum = 0; psfnum < maxboxno; psfnum++) {

        /*
         * downsample by a factor 1  and fill up array of size pwidth 
         */

        dfactor = 1;
        dwidth = pwidth;        /* pwidth/dfactor; */

        for(iy = 0; iy < pheight; iy++) {
            for(ix = 0; ix < pwidth; ix++) {
                ij = iy * pwidth + ix;
                psf[ij] = 0.;
            }
        }

        /*
         * swap quadrants in the data so that the psf is centred 
         * on pixel (0,0) 
         */

        npixels = 0;

        /*
         * fill up quadrant 1 
         */
        y = pwidth / 2;
        iy = 0;
        while (y < pwidth) {
            x = pwidth / 2;
            ix = 0;
            while (x < pwidth) {
                j = x + y * pwidth;
                ij = iy * pwidth + ix;
                psf[ij] = (float) f[psfnum][j];
                x += dfactor;
                ix++;
                npixels++;
            }
            y += dfactor;
            iy++;
        }

        /*
         * fill up quadrant 2 
         */
        y = pwidth / 2;
        iy = 0;
        while (y < pwidth) {
            x = 0;
            ix = pwidth - dwidth / 2;
            while (x < (pwidth / 2)) {
                j = x + y * pwidth;
                ij = iy * pwidth + ix;
                psf[ij] = (float) f[psfnum][j];
                x += dfactor;
                ix++;
                npixels++;
            }
            y += dfactor;
            iy++;
        }

        /*
         * fill up quadrant 3 
         */
        y = 0;
        iy = pheight - dwidth / 2;
        while (y < (pwidth / 2)) {
            x = 0;
            ix = pwidth - dwidth / 2;
            while (x < (pwidth / 2)) {
                j = x + y * pwidth;
                ij = iy * pwidth + ix;
                psf[ij] = (float) f[psfnum][j];
                x += dfactor;
                ix++;
                npixels++;
            }
            y += dfactor;
            iy++;
        }

        /*
         * fill up quadrant 4 
         */
        y = 0;
        iy = pheight - dwidth / 2;
        while (y < (pwidth / 2)) {
            x = pwidth / 2;
            ix = 0;
            while (x < pwidth) {
                j = x + y * pwidth;
                ij = iy * pwidth + ix;
                psf[ij] = (float) f[psfnum][j];
                x += dfactor;
                ix++;
                npixels++;
            }
            y += dfactor;
            iy++;
        }

        if (npixels != (dwidth * dwidth)) {
            printf(" npixels does not match \n ");
            exit(2);
        }


        /*
         * write psf into 3D array 
         */
        for(i = 0; i < npixels; i++) {
            j = psfnum * npixels + i;
            psf3d[j] = psf[i];
        }

    }

    /*
     * write all PSF data into image array 
     */

    if (fits_write_pix(afptr, TFLOAT, fpixel, psfsize3d, psf3d, &status)) {
        printf(" error reading pixel data \n");
        exit(2);
    }

    /*
     * write which stars have been used for each PSF into a FITS table 
     */

    maxlistnum = 0;
    for(psfnum = 0; psfnum < maxboxno; psfnum++) {
        if (listnum[psfnum] > maxlistnum)
            maxlistnum = listnum[psfnum];
    }

    for(i = 0; i < 5; i++) {
        tform[i] = (char *) calloc(8, sizeof(char));
    }

    sprintf(tform[0], "1J");
    sprintf(tform[1], "1J");
    sprintf(tform[2], "1J");
    sprintf(tform[3], "%dE", maxlistnum);
    sprintf(tform[4], "%dE", maxlistnum);

    fits_create_tbl(afptr, BINARY_TBL, maxboxno, tfields, ttype, tform, tunit,
                    extname, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    for(i = 0; i < maxboxno; i++) {
        psfnumy = i / xboxnum;
        psfnumx = i - psfnumy * xboxnum;
        psfnumx += 1;
        psfnumy += 1;
        row = i + 1;
        fits_write_col(afptr, TINT, 1, row, 1, 1, &psfnumx, &status);
        fits_write_col(afptr, TINT, 2, row, 1, 1, &psfnumy, &status);
        fits_write_col(afptr, TINT, 3, row, 1, 1, &listnum[i], &status);
        if (listnum[i] > 0) {
            fits_write_col(afptr, TFLOAT, 4, row, 1, listnum[i], listx[i],
                           &status);
            fits_write_col(afptr, TFLOAT, 5, row, 1, listnum[i], listy[i],
                           &status);
        }
    }

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * write keywords into header 
     */

    if (fits_update_key(afptr, TINT, "NXPSF", &xboxnum,
                        "number of PSFs along x axis", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NYPSF", &yboxnum,
                        "number of PSFs along y axis", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "PSFBOX", &boxsize,
                        "PSF box size (pixels)", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_update_key(afptr, TINT, "PSMOOTH", &boxsmooth,
                        "PSF smoothing length (pixels)", &status)) {
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
