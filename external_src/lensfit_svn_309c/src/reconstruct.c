#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

void
reconstruct(double xval, double yval, int order, int crossterm,
            double *acoeffs, double *fval)
{
    int i, j, k;
    double ipow, jpow;

    if (order < 0) {
        order = abs(order);
        crossterm = 1;
    }

    if (crossterm == 0) {
        fval[0] = 0;
        k = 0;
        for(i = 0; i <= order; i++) {
            ipow = (double) (i);
            for(j = 0; j <= (order - i); j++) {
                jpow = (double) (j);
                fval[0] += acoeffs[k] * pow(xval, ipow) * pow(yval, jpow);
                k++;
            }
        }
    } else {
        fval[0] = 0;
        k = 0;
        for(i = 0; i <= order; i++) {
            ipow = (double) (i);
            for(j = 0; j <= order; j++) {
                jpow = (double) (j);
                fval[0] += acoeffs[k] * pow(xval, ipow) * pow(yval, jpow);
                k++;
            }
        }
    }

}
