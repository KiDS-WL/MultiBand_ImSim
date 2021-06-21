#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define NRANSI
#include "nrutil.h"
#define TOL 1.0e-6

void
svdfit2dsqc(double x[], double y[], double z[], double weight[], int nobj,
            int order, int crossterm, double a[], double **u, double **v,
            double w[])
{
    void dsvbksb(double **, double[], double **, int, int, double[],
                 double[]);
    void dsvdcmp(double **, int, int, double[], double **);

    int j, i, m, n, k, p, q;
    double wmax, thresh, b[100], f, ipow, jpow;

/* code to accumulate a square least-squares matrix and use the numerical
recipes svd to solve */

    if (crossterm == 0) {

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
                        f = 0;
                        for(m = 0; m < nobj; m++) {
                            f += weight[m] * pow(x[m], ipow) * pow(y[m],
                                                                   jpow);
                        }
                        u[n][k] = f;
                    }
                }
            }
        }


        k = 0;
        for(i = 0; i <= order; i++) {
            ipow = (double) (i);
            for(j = 0; j <= (order - i); j++) {
                jpow = (double) (j);
                k++;
                f = 0;
                for(m = 0; m < nobj; m++) {
                    f += weight[m] * z[m] * pow(x[m], ipow) * pow(y[m], jpow);
                }
                b[k] = f;
            }
        }

    } else {

        n = 0;
        for(p = 0; p <= order; p++) {
            for(q = 0; q <= order; q++) {
                n++;
                k = 0;
                for(i = 0; i <= order; i++) {
                    ipow = (double) (i + p);
                    for(j = 0; j <= order; j++) {
                        jpow = (double) (q + j);
                        k++;
                        f = 0;
                        for(m = 0; m < nobj; m++) {
                            f += weight[m] * pow(x[m], ipow) * pow(y[m],
                                                                   jpow);
                        }
                        u[n][k] = f;
                    }
                }
            }
        }


        k = 0;
        for(i = 0; i <= order; i++) {
            ipow = (double) (i);
            for(j = 0; j <= order; j++) {
                jpow = (double) (j);
                k++;
                f = 0;
                for(m = 0; m < nobj; m++) {
                    f += weight[m] * z[m] * pow(x[m], ipow) * pow(y[m], jpow);
                }
                b[k] = f;
            }
        }

    }

    dsvdcmp(u, n, n, w, v);

    wmax = 0.0;
    for(j = 1; j <= n; j++)
        if (w[j] > wmax)
            wmax = w[j];
    thresh = TOL * wmax;
    for(j = 1; j <= n; j++)
        if (w[j] < thresh)
            w[j] = 0.0;

    dsvbksb(u, w, v, n, n, b, a);

}

#undef TOL
#undef NRANSI
