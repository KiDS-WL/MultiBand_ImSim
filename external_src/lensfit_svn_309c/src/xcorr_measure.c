#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

#include <complex.h>
#include <fftw3.h>

int
xcorr_measure(fftw_plan padinv, fftw_complex * B, fftw_complex * D,
              int pheight, int pwidth, int padheight, int padwidth,
              fftw_complex * C, double *c, fftw_complex * Cpad, double *shift)
{
    double maxval, xx, yy;
    double realpart, impart;
    double rsq, xd, yd;
    double rangemax, rangemaxsq;
    double fsize;
    int i, x, y, xy, yx;
    int xmax, ymax;
    int halfpwidth, halfpheight, halfpadwidth, newxy;
    int dsfactor;

    /*
     * measures the shift between two stars/galaxies based on routine xcorr.  Inputs
     * are the FTs of the two stars, returned output is the shift in x,y 
     */

    /*
     * calculate the expansion factor 
     */
    dsfactor = padheight / pheight;

    /*
     * the maximum of the cross-correlation is searched within a radius
     * rangemax of the nominal position at (0,0) in units of the expanded grid
     */
    rangemax = dsfactor * 5.;
    rangemaxsq = rangemax * rangemax;

    fsize = (double) pheight *pwidth;

    halfpwidth = 1 + pwidth / 2;
    halfpheight = pheight / 2;

    i = 0;

    for(y = 0; y < pheight; y++) {
        for(x = 0; x < halfpwidth; x++) {
            xy = y * halfpwidth + x;
            realpart = creal(B[xy]) / fsize;
            impart = cimag(B[xy]) / fsize;
            C[xy] = (realpart * creal(D[xy]) + impart * cimag(D[xy])) +
                (realpart * cimag(D[xy]) - impart * creal(D[xy])) * I;
        }
    }

    /*
     * pad to make a higher-resolution image 
     */

    halfpadwidth = 1 + padwidth / 2;

    for(i = 0; i < (halfpadwidth * padheight); i++)
        Cpad[i] = 0. + 0. * I;

    for(y = 0; y < halfpheight; y++) {
        for(x = 0; x < halfpwidth; x++) {
            xy = y * halfpwidth + x;
            newxy = y * halfpadwidth + x;
            Cpad[newxy] = C[xy];
        }
    }

    /*
     * deal with Nyquist frequency 
     */
    y = halfpheight;
    for(x = 0; x < halfpwidth; x++) {
        xy = y * halfpwidth + x;
        newxy = y * halfpadwidth + x;
        Cpad[newxy] = C[xy] / 2.;
    }
    /*
     * top quadrant 
     */
    for(y = (halfpheight + 1); y < pheight; y++) {
        for(x = 0; x < halfpwidth; x++) {
            xy = y * halfpwidth + x;
            newxy = (y + padheight - pheight) * halfpadwidth + x;
            Cpad[newxy] = C[xy];
        }
    }

    /*
     * deal with Nyquist frequency 
     */
    y = halfpheight;
    for(x = 0; x < halfpwidth; x++) {
        xy = y * halfpwidth + x;
        newxy = (y + padheight - pheight) * halfpadwidth + x;
        Cpad[newxy] = C[xy] / 2.;
    }

    fftw_execute(padinv);

    /*
     * now locate the maximum. 
     */

    maxval = c[0];
    xmax = 0;
    ymax = 0;

    for(y = 0; y < rangemax; y++) {
        yd = (double) y;
        yy = yd * yd;
        for(x = 0; x < rangemax; x++) {
            xd = (double) x;
            xx = xd * xd;
            rsq = xx + yy;
            if (rsq < rangemaxsq) {
                yx = y * padwidth + x;
                if (c[yx] > maxval) {
                    xmax = x;
                    ymax = y;
                    maxval = c[yx];
                }
            }
        }
        for(x = (padwidth - rangemax); x < padwidth; x++) {
            xd = (double) (x - padwidth);
            xx = xd * xd;
            rsq = xx + yy;
            if (rsq < rangemaxsq) {
                yx = y * padwidth + x;
                if (c[yx] > maxval) {
                    xmax = x - padwidth;
                    ymax = y;
                    maxval = c[yx];
                }
            }
        }
    }
    for(y = (padheight - rangemax); y < padheight; y++) {
        yd = (double) (y - padheight);
        yy = yd * yd;
        for(x = 0; x < rangemax; x++) {
            xd = (double) x;
            xx = (double) (x * x);
            rsq = xx + yy;
            if (rsq < rangemaxsq) {
                yx = y * padwidth + x;
                if (c[yx] > maxval) {
                    xmax = x;
                    ymax = y - padheight;
                    maxval = c[yx];
                }
            }
        }
        for(x = (padwidth - rangemax); x < padwidth; x++) {
            xd = (double) (x - padwidth);
            xx = xd * xd;
            rsq = xx + yy;
            if (rsq < rangemaxsq) {
                yx = y * padwidth + x;
                if (c[yx] > maxval) {
                    xmax = x - padwidth;
                    ymax = y - padheight;
                    maxval = c[yx];
                }
            }
        }
    }

    shift[0] = xmax * (float) pwidth / (float) padwidth;
    shift[1] = ymax * (float) pwidth / (float) padwidth;

    return 0;

}
