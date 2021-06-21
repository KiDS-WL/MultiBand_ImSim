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
shiftft(int pheight, int pwidth, fftw_complex * C, double *shift)
{
    double fsize, x0, y0, kx, ky, pi;
    int x, y, yx;
    int halfpwidth, halfpheight;

    /*
     * applies a shift in x,y (shift[0]. shift[1]) to the Fourier transform C 
     * of a postage stamp of size pwidth, pheight assuming C has been 
     * generated from an fftw real-to-complex transform 
     */

    pi = 3.141592654;

    fsize = (double) pheight *pwidth;
    halfpwidth = 1 + pwidth / 2;
    halfpheight = pheight / 2;

    /*
     * careful of the sign convention defined for the shift! 
     */
    x0 = -shift[0];
    y0 = -shift[1];

    if (x0 != 0. || y0 != 0.) {

        for(y = 0; y < halfpheight; y++) {
            ky = pi * (double) y / (double) (halfpheight);
            x = 0;
            yx = y * halfpwidth + x;
            C[yx] = C[yx] * cos(ky * y0) + C[yx] * sin(ky * y0) * I;
            for(x = 1; x < halfpwidth - 1; x++) {
                yx = y * halfpwidth + x;
                kx = pi * (double) x / (double) (halfpwidth - 1);
                C[yx] = C[yx] * cos(kx * x0) + C[yx] * sin(kx * x0) * I;
                C[yx] = C[yx] * cos(ky * y0) + C[yx] * sin(ky * y0) * I;
            }
            x = halfpwidth - 1;
            yx = y * halfpwidth + x;
            kx = pi;
            C[yx] = C[yx] * cos(kx * x0) + C[yx] * sin(kx * x0) * I;
            C[yx] = C[yx] * cos(ky * y0) + C[yx] * sin(ky * y0) * I;
        }
        for(y = halfpheight; y < pheight; y++) {
            ky = pi * (double) (y - pheight) / (double) (halfpheight);
            x = 0;
            yx = y * halfpwidth + x;
            C[yx] = C[yx] * cos(ky * y0) + C[yx] * sin(ky * y0) * I;
            for(x = 1; x < halfpwidth - 1; x++) {
                yx = y * halfpwidth + x;
                kx = pi * (double) x / (double) (halfpwidth - 1);
                C[yx] = C[yx] * cos(kx * x0) + C[yx] * sin(kx * x0) * I;
                C[yx] = C[yx] * cos(ky * y0) + C[yx] * sin(ky * y0) * I;
            }
            x = halfpwidth - 1;
            yx = y * halfpwidth + x;
            kx = pi;
            C[yx] = C[yx] * cos(kx * x0) + C[yx] * sin(kx * x0) * I;
            C[yx] = C[yx] * cos(ky * y0) + C[yx] * sin(ky * y0) * I;
        }
    }

    return 0;

}
