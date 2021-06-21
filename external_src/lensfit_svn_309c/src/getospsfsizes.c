#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

/* this version also gets keywords for shifts fits */

int
getospsfsizes(char *psfname, int *correct_distortion, int *fwidth, int *fheight, 
	      int *oversampling, int *order,
	      int *chipvariation, int *chiporder, int *ncoeffs, int *nfit,
          int *numstar, int *npsfchip, int *xchip, int *ychip, int *xchipsampling,
	      int *ychipsampling, int *big_gap, double *hxsize,
	      double *hysize, float *poserror, double *scalefactor, int *sorder,
	      int *schipvariation, int *schiporder, int *sncoeffs)
{
    fitsfile *afptr;
    int nxchip[1], nychip[1], x, y;
    int status = 0;
    int anaxis, hdunum, hdutype;
    long anaxes[2] = { 1, 1 };
    char comment[64];
    long nrows;
    double dnull;
    float fnull;
    int anynull, ncol;
    int noglobalshift = { 0 };
    long frow, felem, nelem;

    /*
     * open input images 
     */
    fits_open_file(&afptr, psfname, READONLY, &status);
    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * read dimensions 
     */
    fits_get_img_dim(afptr, &anaxis, &status);
    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis != 2) {
        fflush(stdout);
        fprintf(stderr,
                "Error: expected 2 dimensions in PSF fits file\n");
        exit(2);
    }

    fits_get_img_size(afptr, anaxis, anaxes, &status);
    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_read_key(afptr, TINT, "ORDER", order, comment, &status)) {
      fflush(stdout);
      fits_report_error(stderr, status);  /* print error message */
      printf(" keyword ORDER\n");
      return (status);
    }

    if (fits_read_key(afptr, TINT, "NCOEFFS", ncoeffs, comment, &status)) {
      fflush(stdout);
      fits_report_error(stderr, status);  /* print error message */
      printf(" keyword NCOEFFS\n");
      return (status);
    }

    if (fits_read_key(afptr, TINT, "NFIT", nfit, comment, &status)) {
      fflush(stdout);
      fits_report_error(stderr, status);  /* print error message */
      printf(" keyword NFIT\n");
      return (status);
    }

    // get oversampled PSF keywords if they are available
    if (fits_read_key(afptr, TINT, "FWIDTH", fwidth,
                        comment, &status)) {
      *fwidth = 0;
      status = 0;
    }
    if (fits_read_key(afptr, TINT, "FHEIGHT", fheight,
                        comment, &status)) {
      *fheight = 0;
      status = 0;
    }
    if (fits_read_key(afptr, TINT, "OVERSAMP", oversampling,
                        comment, &status)) {
      *oversampling = 0;
      status = 0;
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
        (afptr, TDOUBLE, "PIXSCALE", scalefactor, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        fprintf(stderr, " keyword PIXSCALE\n");
        exit(2);
    }
    if (fits_read_key
        (afptr, TINT, "CHIPVAR", chipvariation, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHIPVAR\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "CHORDER", chiporder, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHORDER\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "NXCHIP", nxchip, comment, &status)) {
        fflush(stdout);
        //fits_report_error(stderr, status);      /* print error message */
        //printf(" keyword NXCHIP\n");
        *nxchip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NYCHIP", nychip, comment, &status)) {
        fflush(stdout);
        //fits_report_error(stderr, status);      /* print error message */
        //printf(" keyword NYCHIP\n");
        *nychip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NPSFCHIP", npsfchip, comment, &status)) {
        fflush(stdout);
        //fits_report_error(stderr, status);      /* print error message */
        //printf(" keyword NYCHIP\n");
        *npsfchip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "XSAMPL", xchipsampling, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword XSAMPL\n");
        *xchipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "YSAMPL", ychipsampling, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword YSAMPL\n");
        *ychipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "BIG_GAP", big_gap, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword BIG_GAP\n");
        *big_gap = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HXSIZE", hxsize, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword HXSIZE\n");
        *hxsize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HYSIZE", hysize, comment, &status)) {
        fflush(stdout);
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

    if (fits_read_key(afptr, TINT, "SORDER", sorder, comment, &status)) {
        fflush(stdout);
        printf
            (" warning from getospsfsizes: globalshifts keyword SORDER has not been found \n");
        //fits_report_error(stderr, status); /* print error message */
        //printf(" keyword SORDER\n");
        status=0;
        noglobalshift = 1;
    }

    if (noglobalshift == 0) {
        if (fits_read_key
            (afptr, TINT, "SCHIPVAR", schipvariation, comment, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword SCHIPVAR\n");
            exit(2);
        }
        if (fits_read_key
            (afptr, TINT, "SCHORDER", schiporder, comment, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword SCHORDER\n");
            exit(2);
        }
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
            fprintf(stderr," error in chip format in PSF file \n");
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
      

    /*
     * read table of star numbers in each chip 
     */
    /*
     * get the table of star info 
     */
    hdunum++;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        status = 0;
        fflush(stdout);
        fprintf(stderr, " error moving to table of star information \n");
        fits_report_error(stderr, status);      /* print error message */
        status=0;
	return(0);
    } else {
        if (hdutype != BINARY_TBL) {
            fflush(stdout);
            fprintf(stderr,
                    " Error file format not as expected - Table of star info is not a binary table \n");
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
                    " error, number of table entries does not match array size \n");
            exit(2);
        }

        frow = 1;
        felem = 1;
        nelem = npsfchip[0];
        dnull = 0.;
        fnull = 0.;
        anynull = 0;
        if (fits_read_col
            (afptr, TINT, 2, frow, felem, nelem, &dnull, numstar, &anynull,
             &status)) {
            fflush(stdout);
            fprintf(stderr, " error reading dec column in table \n");
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }
    }
      

    /*
     * repeat for x,y shifts coefficients 
     */
    /*
     * read table of star numbers in each chip 
     */
    /*
     * get the table of star info 
     */

    if (noglobalshift == 0) {
        hdunum++;
        if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
            status = 0;
            fflush(stdout);
            fprintf(stderr,
                    " error moving to table of PSF shifts information \n");
            fits_report_error(stderr, status);  /* print error message */
            fprintf(stderr, " run program globalshifts first \n");
            exit(2);
        } else {
            if (hdutype != BINARY_TBL) {
                fflush(stdout);
                fprintf(stderr,
                        " Error file format not as expected - Table of star info is not a binary table \n");
                exit(2);
            }
            if (fits_get_num_cols(afptr, &ncol, &status)) {
                fflush(stdout);
                fits_report_error(stderr, status);      /* print error message */
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
                fits_report_error(stderr, status);      /* print error message */
                exit(2);
            }

            *sncoeffs = nrows;
        }
    }
       
    
    fits_close_file(afptr, &status);

    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    return (0);

}
