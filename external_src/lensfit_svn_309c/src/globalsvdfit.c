#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define NRANSI
#include "nrutil.h"
#define TOL 1.0e-6

void
globalsvdfit(double x[], double y[], double z[], double weight[], int chip[],
             int nobj, int nchip, int order, int crossterm, int chipvariation,
             int chiporder, int ncoeff, double a[], double **u, double **v,
             double w[])
{
    void dsvbksb(double **, double[], double **, int, int, double[],
                 double[]);
    void dsvdcmp(double **, int, int, double[], double **);

    int j, i, m, n, k, p, q, col3, nc, **kchip, nelm, kpos, npos;
    double wmax, thresh, b[1000], f, ipow, jpow;

/* code to accumulate a square least-squares matrix and use the numerical
   recipes svd to solve for a polynomial fit.
   This routine assumes a global fit across many chip datasets, but with low order
   terms allowed to have independent local values for each chip.  

   x, y, z, weight are values for each object
   *chip is an array that specifies the chip that each object lies on
   nobj  is the number of objects
   nchip is the total number of chips
   order is the global fit order
   crossterm specifies whether global crossterms are allowed (if=1)
   chiporder is the chip-dependent fit order (negative implies crossterms allowed)
   chipvariation is a flag that specifies if chip variations are allowed (if=1)
   ncoeff is the total number of coefficients (i.e. dimension of arrays)
   u, v, w are arrays used by Numerical Recipes SVD routine

   If more chip-dependent terms are specified than are allowed by the number
   of global terms specified, the extra chip-dependent terms will just be ignored.
*/

    if (nobj <= 0 || nchip <= 0) {
        fprintf(stderr, " error, no data passed to globalsvdfit \n");
        exit(2);
    }

    if (order < 0) {
        fprintf(stderr,
                " negative global order specified in globalsvdfit \n");
        exit(2);
    }

    /*
     * work out number of coefficients required and allocate memory 
     */
    nelm =
        crossterm ==
        1 ? (1 + order) * (1 + order) : (1 + order) * (2 + order) / 2;
    nelm++;                     // allow for NR array convention
    kchip = (int **) calloc(nelm, sizeof(int *));
    if (kchip == NULL) {
        fprintf(stderr, " error allocating memory inside globalsvdfit \n");
        exit(2);
    }
    for(i = 0; i < nelm; i++) {
        kchip[i] = (int *) calloc(nchip, sizeof(int));
        if (kchip[i] == NULL) {
            fprintf(stderr,
                    " error allocating memory inside globalsvdfit \n");
            exit(2);
        }
    }


    /*
     * if 0 don't include high order cross-terms (normal case) 
     */
    if (crossterm == 0) {
        /*
         * work out the position in the matrix for each coefficient.
         * If a coefficient is allowed to have variation between chips,
         * there need to be nchip columns/rows for that coefficient,
         * otherwise just one row.  kchip is an array that allows
         * a lookup of which is the correct column/row for the mth
         * chip dataset 
         */
        /*
         * note on confusing feature: for consistency with the NUmerical
         * Recipes routine, array elements are numbered starting at 1. The 
         * loops p,q etc that define the powers in x,y are numbered 
         * starting at 0 
         */
        if (chiporder < 0) {
            /*
             * cross-terms are allowed in the chip-dependent coefficients 
             */
            chiporder = -chiporder;
            n = 0;
            nc = 0;
            for(p = 0; p <= order; p++) {
                for(q = 0; q <= (order - p); q++) {
                    n++;
                    /*
                     * test if this coefficient is chip-dependent 
                     */
                    if (chipvariation == 1 && p <= chiporder
                        && q <= chiporder) {
                        for(i = 0; i < nchip; i++) {
                            nc++;
                            kchip[n][i] = nc;
                        }
                    } else {
                        nc++;
                        for(i = 0; i < nchip; i++) {
                            kchip[n][i] = nc;
                        }
                    }
                }
            }
        } else {
            /*
             * chip-dependent cross-terms not allowed (such terms are global) 
             */
            n = 0;
            nc = 0;
            for(p = 0; p <= order; p++) {
                for(q = 0; q <= (order - p); q++) {
                    n++;
                    /*
                     * test if this coefficient is chip-dependent 
                     */
                    if (chipvariation == 1 && (p + q) <= chiporder) {
                        for(i = 0; i < nchip; i++) {
                            nc++;
                            kchip[n][i] = nc;
                        }
                    } else {
                        nc++;
                        for(i = 0; i < nchip; i++) {
                            kchip[n][i] = nc;
                        }
                    }
                }
            }
        }

        if (nc > ncoeff) {
            fprintf(stderr, " error in globalsvdfit \n");
            fprintf(stderr,
                    " not enough array elements declared for required order \n");
            exit(2);
        }

        /*
         * initialise svd arrays (remembering NR convention on array elemnents) 
         */
        for(j = 0; j <= nc; j++)
            b[j] = 0.;
        for(i = 0; i <= nc; i++)
            for(j = 0; j <= nc; j++)
                u[i][j] = 0.;


        /*
         * first accumulate matrices for each chip as if they were
         * being fitted independently 
         */

        /*
         * LHS 
         */
        n = 0;
        for(p = 0; p <= order; p++) {
            for(q = 0; q <= (order - p); q++) {
                n++;
                k = 0;
                for(i = 0; i <= order; i++) {
                    ipow = (double) (i + p);
                    for(j = 0; j <= (order - i); j++) {
                        jpow = (double) (q + j);
                        k++;
                        if (ipow == 1 && jpow == 0 && n == 1)
                            col3 = k;
                        for(m = 0; m < nobj; m++) {
                            f = weight[m] * pow(x[m], ipow) * pow(y[m], jpow);
                            npos = kchip[n][chip[m]];
                            kpos = kchip[k][chip[m]];
                            u[npos][kpos] += f;
                        }
                    }
                }
            }
        }

        /*
         * RHS 
         */
        k = 0;
        for(i = 0; i <= order; i++) {
            ipow = (double) (i);
            for(j = 0; j <= (order - i); j++) {
                jpow = (double) (j);
                k++;
                for(m = 0; m < nobj; m++) {
                    f = weight[m] * z[m] * pow(x[m], ipow) * pow(y[m], jpow);
                    kpos = kchip[k][chip[m]];
                    b[kpos] += f;
                }
            }

        }

    } else {
        /*
         * global cross-terms are allowed 
         */

        /*
         * work out the position in the matrix for each coefficient.
         * If a coefficient is allowed to have variation between chips,
         * there need to be nchip columns/rows for that coefficient,
         * otherwise just one row.  kchip is an array that allows
         * a lookup of which is the correct column/row for the mth
         * chip dataset 
         */
        /*
         * note on confusing feature: for consistency with the NUmerical
         * Recipes routine, array elements are numbered starting at 1. The 
         * loops p,q etc that define the powers in x,y are numbered 
         * starting at 0 
         */
        if (chiporder < 0) {
            /*
             * cross-terms are allowed in the chip-dependent coefficients 
             */
            chiporder = -chiporder;
            n = 0;
            nc = 0;
            for(p = 0; p <= order; p++) {
                for(q = 0; q <= (order - p); q++) {
                    n++;
                    /*
                     * test if this coefficient is chip-dependent 
                     */
                    if (chipvariation == 1 && p <= chiporder
                        && q <= chiporder) {
                        for(i = 0; i < nchip; i++) {
                            nc++;
                            kchip[n][i] = nc;
                        }
                    } else {
                        nc++;
                        for(i = 0; i < nchip; i++) {
                            kchip[n][i] = nc;
                        }
                    }
                }
            }
        } else {
            /*
             * chip-dependent cross-terms not allowed (such terms are global) 
             */
            n = 0;
            nc = 0;
            for(p = 0; p <= order; p++) {
                for(q = 0; q <= (order - p); q++) {
                    n++;
                    /*
                     * test if this coefficient is chip-dependent 
                     */
                    if (chipvariation == 1 && (p + q) <= chiporder) {
                        for(i = 0; i < nchip; i++) {
                            nc++;
                            kchip[n][i] = nc;
                        }
                    } else {
                        nc++;
                        for(i = 0; i < nchip; i++) {
                            kchip[n][i] = nc;
                        }
                    }
                }
            }
        }

        if (nc > ncoeff) {
            fprintf(stderr, " error in globalsvdfit \n");
            fprintf(stderr,
                    " not enough array elements declared for required order \n");
            exit(2);
        }

        /*
         * initialise svd arrays (remembering NR convention on array elemnents) 
         */
        for(j = 0; j <= nc; j++)
            b[j] = 0.;
        for(i = 0; i <= nc; i++)
            for(j = 0; j <= nc; j++)
                u[i][j] = 0.;


        /*
         * first accumulate matrices for each chip as if they were
         * being fitted independently 
         */

        /*
         * LHS 
         */
        n = 0;
        for(p = 0; p <= order; p++) {
            for(q = 0; q <= (order - p); q++) {
                n++;
                k = 0;
                for(i = 0; i <= order; i++) {
                    ipow = (double) (i + p);
                    for(j = 0; j <= (order - i); j++) {
                        jpow = (double) (q + j);
                        k++;
                        if (ipow == 1 && jpow == 0 && n == 1)
                            col3 = k;
                        for(m = 0; m < nobj; m++) {
                            f = weight[m] * pow(x[m], ipow) * pow(y[m], jpow);
                            npos = kchip[n][chip[m]];
                            kpos = kchip[k][chip[m]];
                            u[npos][kpos] += f;
                        }
                    }
                }
            }
        }

        /*
         * RHS 
         */
        k = 0;
        for(i = 0; i <= order; i++) {
            ipow = (double) (i);
            for(j = 0; j <= (order - i); j++) {
                jpow = (double) (j);
                k++;
                for(m = 0; m < nobj; m++) {
                    f = weight[m] * z[m] * pow(x[m], ipow) * pow(y[m], jpow);
                    kpos = kchip[k][chip[m]];
                    b[kpos] += f;
                }
            }

        }
    }


    /*
     * free memory 
     */
    for(i = 0; i < nelm; i++)
        free(kchip[i]);
    free(kchip);

    /*
     * use NR SVD routines to solve 
     */

    for(i = 1; i <= nc; i++)
        for(j = 1; j <= nc; j++)
            v[i][j] = 0.;
    for(j = 1; j <= nc; j++)
        w[j] = 0.;

    dsvdcmp(u, nc, nc, w, v);

    /*
     * eliminate singular values (this should take care of missing chips) 
     */
    wmax = 0.0;
    for(j = 1; j <= nc; j++)
        if (w[j] > wmax)
            wmax = w[j];
    thresh = TOL * wmax;
    for(j = 1; j <= nc; j++)
        if (w[j] < thresh)
            w[j] = 0.0;

    dsvbksb(u, w, v, nc, nc, b, a);

}

#undef TOL
#undef NRANSI
