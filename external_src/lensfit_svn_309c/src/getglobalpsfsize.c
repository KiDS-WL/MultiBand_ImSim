#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
getglobalpsfsize(char *psfname, int *correct_distortion, int *pwidth, int *pheight, int *order,
                 int *chipvariation, int *chiporder, int *ncoeffs,
                 int *npsfchip, int *xchip, int *ychip, int *xchipsampling,
                 int *ychipsampling, int *big_gap, double *hxsize,
                 double *hysize, float *poserror)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis, hdunum, hdutype;
    int nxchip[1], nychip[1], x, y;
    long anaxes[3] = { 1, 1, 1 };
    char comment[FLEN_COMMENT];
    long nrows;
    double dnull;
    float fnull;
    int anynull, ncol;
    long frow, felem, nelem;

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
     * get remaining keywords 
     */
    if (fits_read_key(afptr, TINT, "CORRECTD", correct_distortion, comment, &status)) {
      fflush(stdout);
      fits_report_error(stderr, status);  /* print error message */
      fprintf(stderr," keyword CORRECTD\n");
      fprintf(stderr," distortion correction keyword\n");
      return (status);
    }
    if (fits_read_key
        (afptr, TINT, "CHIPVAR", chipvariation, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHIPVAR\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "CHORDER", chiporder, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHORDER\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "NXCHIP", nxchip, comment, &status)) {
      //fits_report_error(stderr, status);      /* print error message */
      //printf(" keyword NXCHIP\n");
        *nxchip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NYCHIP", nychip, comment, &status)) {
      //fits_report_error(stderr, status);      /* print error message */
      //printf(" keyword NYCHIP\n");
        *nychip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NPSFCHIP", npsfchip, comment, &status)) {
      //fits_report_error(stderr, status);      /* print error message */
      //printf(" keyword NPSFCHIP\n");
        *npsfchip = 0;
        status = 0;
    }    
    if (fits_read_key(afptr, TINT, "XSAMPL", xchipsampling, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword XSAMPL\n");
        *xchipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "YSAMPL", ychipsampling, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword YSAMPL\n");
        *ychipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "BIG_GAP", big_gap, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword BIG_GAP\n");
        *big_gap = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HXSIZE", hxsize, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword HXSIZE\n");
        *hxsize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HYSIZE", hysize, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword HYSIZE\n");
        *hysize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TFLOAT, "TOL", poserror, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword TOL\n");
        *poserror = 0;
        status = 0;
    }

    // set current FITS table number
    hdunum = 1;
    
    // if old style PSF file with rectilinear chip arrangement, fill up the appropriate arrays
    if (npsfchip[0]==0)
      {
        if (nxchip[0] != 0 && nychip[0] != 0)
          {
            npsfchip[0] = 0;
            for (y=0; y<nychip[0]; y++)
              {
                for (x=0; x<nxchip[0]; x++)
                  {
                    xchip[npsfchip[0]] = x;
                    ychip[npsfchip[0]] = y;
                    npsfchip[0]++;
                  }
              }
          }
        else
          {
            fflush(stdout);
            fprintf(stderr," error in chip format in PSF file header \n");
            exit(EXIT_FAILURE);
          }
      }
    else
      {
        // read table of chip position identifiers
        hdunum++;
        if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
          status = 0;
          fflush(stdout);
          fprintf(stderr, " error moving to table of chip position information \n");
          fits_report_error(stderr, status);      /* print error message */
          status = 0;
          return(0);
        } else {
          if (hdutype != BINARY_TBL) {
            fflush(stdout);
            fprintf(stderr,
                    " Error file format not as expected - Table of chip position info is not a binary table \n");
            exit(2);
          }
          if (fits_get_num_cols(afptr, &ncol, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            return (status);
          }
          if (ncol != 2) {
            fflush(stdout);
            fprintf(stderr, " error in table format, no. columns = %d \n",
                    ncol);
            exit(2);
          }
          if (fits_get_num_rows(afptr, &nrows, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
          }
          if (nrows != npsfchip[0]) {
            fflush(stdout);
            fprintf(stderr,
                    " error, number of chip position entries does not match array size \n");
            exit(2);
          }
        frow = 1;
        felem = 1;
        nelem = npsfchip[0];
        dnull = 0.;
        fnull = 0.;
        anynull = 0;
        if (fits_read_col
            (afptr, TINT, 1, frow, felem, nelem, &dnull, xchip, &anynull,
             &status)) {
            fflush(stdout);
            fprintf(stderr, " error reading chip x position column in table \n");
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }
        if (fits_read_col
            (afptr, TINT, 2, frow, felem, nelem, &dnull, ychip, &anynull,
             &status)) {
            fflush(stdout);
            fprintf(stderr, " error reading chip y position column in table \n");
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }
        
        }
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
