#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

#include <complex.h>
#include <fftw3.h>

void
convolvemodel(fftw_plan * pinv,
              fftw_complex * B, float **modelft,
              int pheight, int pwidth, fftw_complex **C)
{
    double realpart, impart;
    double fsize;
    int x, y, yx, mt, halfpwidth, halfpheight;


    fsize = (double) pheight *pwidth;
    halfpwidth = 1 + pwidth / 2;

    halfpheight = pheight / 2;

    /*
     * now calculate the relevant data-model crosscorrelation and model
     * autocorrelations 
     */

    for(mt = 0; mt < 2; mt++) {
        C[mt][0] = 0. + 0. * I;

        for(y = 0; y < pheight; y++) {
            for(x = 0; x < halfpwidth; x++) {
                yx = y * halfpwidth + x;
                /*
                 * imaginary part of model is zero 
                 */
                realpart = modelft[mt][yx] * creal(B[yx]) / fsize;
                impart = modelft[mt][yx] * cimag(B[yx]) / fsize;
                C[mt][yx] = realpart + impart * I;
            }
        }

        fftw_execute(pinv[mt]);

    }


}
