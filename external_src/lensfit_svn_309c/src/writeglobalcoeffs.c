#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
writeglobalcoeffs(char *psfname, int correct_distortion, double **f, int ncoeffs, int order,
                  int chipvariation, int chiporder, double *scalefactor,
                  int pwidth, int pheight, int npsfchip, int *xchip, int *ychip,
                  int *numstar, int xchipsampling, int ychipsampling,
                  int big_gap, double hxsize, double hysize, float poserror)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[3], fpixel[3] = { 1, 1, 1 };
    int i, j, ii, bitpix;
    double *psf3d;
    int psfsize, psfsize3d, row;

    int posfields = 2;            /* number of table columns in table */
    char posextname[] = "chip position identifiers";        /* extension name */
    char *postype[] = { "XCHIP", "YCHIP" };        /* names of columns */
    char *posform[] = { "1J", "1J" };     /* data types of columns */
    char *posunit[] = { "", "" }; /* units */
    
    int tfields = 2;            /* number of table columns in table */
    char extname[] = "PSF star numbers";        /* extension name */
    char *ttype[] = { "CHIP", "NSTAR" };        /* names of columns */
    char *tform[] = { "1J", "1J" };     /* data types of columns */
    char *tunit[] = { "", "" }; /* units */

    /*
     * write out psf data as a fits file 
     */

    if (pwidth != pheight) {
        fprintf(stderr, " PSF subimages must be square! \n");
        return 1;
    }

    if (pwidth <= 0 || pheight <= 0) {
        fprintf(stderr, " error in subimage dimensions \n");
        return 1;
    }

    bitpix = -64;
    anaxis = 3;
    anaxes[0] = pwidth;
    anaxes[1] = pheight;
    anaxes[2] = ncoeffs;

    psfsize = anaxes[0] * anaxes[1];
    psfsize3d = ncoeffs * psfsize;
    psf3d = (double *) calloc(psfsize3d, sizeof(double));
    if (psf3d == NULL) {
        fprintf(stderr,
                " error allocating memory for 3D psf coeffs in writecoeffs_fits \n");
        fprintf(stderr, " memory requested = %d doubles \n", psfsize3d);
        exit(2);
    }

    fits_create_file(&afptr, psfname, &status);
    fits_create_img(afptr, bitpix, anaxis, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    for(i = 0; i < ncoeffs; i++) {
        for(j = 0; j < psfsize; j++) {
            ii = i * psfsize + j;
            psf3d[ii] = f[j][i];
        }
    }

    /*
     * write all PSF data into image array 
     */

    if (fits_write_pix(afptr, TDOUBLE, fpixel, psfsize3d, psf3d, &status)) {
        printf(" error writing psf coefficient data \n");
        exit(2);
    }

    /*
     * write keywords into header 
     */

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

    if (fits_update_key(afptr, TINT, "PWIDTH", &pwidth,
                        "postage stamp width", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "PHEIGHT", &pheight,
                        "postage stamp height", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NPSFCHIP", &npsfchip,
                        "number of CCD chips", &status)) {
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

    // chip positions table
    fits_create_tbl(afptr, BINARY_TBL, npsfchip, posfields, postype, posform,
                    posunit, posextname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < npsfchip; i++) {
        row = i + 1;
        printf(" xchip %d ychip %d \n", xchip[i], ychip[i]);
        fits_write_col(afptr, TINT, 1, row, 1, 1, &xchip[i], &status);
        fits_write_col(afptr, TINT, 2, row, 1, 1, &ychip[i], &status);
    }
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    
    // star numbers table
    fits_create_tbl(afptr, BINARY_TBL, npsfchip, tfields, ttype, tform,
                    tunit, extname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < npsfchip; i++) {
        row = i + 1;
        printf(" chip %d numstars %d \n", row, numstar[i]);
        fits_write_col(afptr, TINT, 1, row, 1, 1, &row, &status);
        fits_write_col(afptr, TINT, 2, row, 1, 1, &numstar[i], &status);
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

    free(psf3d);

    return (0);

}
