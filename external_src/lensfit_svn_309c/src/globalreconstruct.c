#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

void
globalreconstruct(double xval, double yval, int order, int crossterm,
                  int chipvariation, int chiporder,
                  int ichip, int nchip, double *acoeffs, double *fval)
{
    /*
     * reconstructs values at specified xval,yval from the supplied
     * polynomial coefficents with decisions on order of fits, and
     * which terms to allow to be chip-dependent, matching what was
     * assumed in globalsvdfit
     */

    int i, j, k, kd = 0, nc, m;
    double ipow, jpow;

    if (order < 0) {
        order = abs(order);
        crossterm = 1;
    }

    fval[0] = 0.;

    if (crossterm == 0) {
        /*
         * global cross-terms not allowed 
         */
        if (chiporder < 0) {
            /*
             * cross-terms are allowed on chip-dependent fits 
             */
            chiporder = -chiporder;
            k = 0;
            nc = 0;
            for(i = 0; i <= order; i++) {
                ipow = (double) (i);
                for(j = 0; j <= (order - i); j++) {
                    jpow = (double) (j);
                    if (chipvariation == 1 && i <= chiporder
                        && j <= chiporder) {
                        for(m = 0; m < nchip; m++) {
                            if (m == ichip)
                                kd = nc;
                            nc++;       // array elements now indexed starting at 0
                        }
                    } else {
                        kd = nc;
                        nc++;
                    }
                    fval[0] +=
                        acoeffs[kd] * pow(xval, ipow) * pow(yval, jpow);
                    k++;
                }
            }
        } else {
            /*
             * cross-terms not allowed for chip-dependent fits 
             */
            k = 0;
            nc = 0;
            for(i = 0; i <= order; i++) {
                ipow = (double) (i);
                for(j = 0; j <= (order - i); j++) {
                    jpow = (double) (j);
                    if (chipvariation == 1 && (i + j) <= chiporder) {
                        for(m = 0; m < nchip; m++) {
                            if (m == ichip)
                                kd = nc;
                            nc++;
                        }
                    } else {
                        kd = nc;
                        nc++;   // array indexed starting at zero
                    }
                    fval[0] +=
                        acoeffs[kd] * pow(xval, ipow) * pow(yval, jpow);
                    k++;
                }
            }

        }

    } else {
        /*
         * global cross-terms allowed 
         */

        if (chiporder < 0) {
            /*
             * cross-terms are allowed on chip-dependent fits 
             */
            chiporder = -chiporder;
            k = 0;
            nc = 0;
            for(i = 0; i <= order; i++) {
                ipow = (double) (i);
                for(j = 0; j <= order; j++) {
                    jpow = (double) (j);
                    if (chipvariation == 1 && i <= chiporder
                        && j <= chiporder) {
                        for(m = 0; m < nchip; m++) {
                            if (m == ichip)
                                kd = nc;
                            nc++;       // array elements now indexed starting at 0
                        }
                    } else {
                        kd = nc;
                        nc++;
                    }
                    fval[0] +=
                        acoeffs[kd] * pow(xval, ipow) * pow(yval, jpow);
                    k++;
                }
            }
        } else {
            /*
             * cross-terms not allowed for chip-dependent fits 
             */
            k = 0;
            nc = 0;
            for(i = 0; i <= order; i++) {
                ipow = (double) (i);
                for(j = 0; j <= order; j++) {
                    jpow = (double) (j);
                    if (chipvariation == 1 && (i + j) <= chiporder) {
                        for(m = 0; m < nchip; m++) {
                            if (m == ichip)
                                kd = nc;
                            nc++;
                        }
                    } else {
                        kd = nc;
                        nc++;   // array indexed starting at zero
                    }
                    fval[0] +=
                        acoeffs[kd] * pow(xval, ipow) * pow(yval, jpow);
                    k++;
                }
            }

        }

    }

}
