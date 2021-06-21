#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void *
ml_calloc(size_t nmemb, size_t size, float *mem, const char *format, ...)
{
    void *ptr;
    va_list ap;

    ptr = calloc(nmemb, size);
    if (!ptr) {
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        exit(EXIT_FAILURE);
    }
    if (mem)
        *mem += (float) nmemb *size;
    return ptr;
}


void
get_poly_order(int order, int chiporder, int chipvariation, int nchip,
               int *ncoeffs, int *crossterm)
{
  // force order to be positive so that no extended cross-terms are allowed
  order = abs(order);
  if (chipvariation==1) chiporder = abs(chiporder);

    if (order >= 0) {
        /*
         * number of coeffs for global fit 
         */
        *ncoeffs = (1 + order) * (2 + order) / 2;
        /*
         * number of extra coeffs for local fits of 3 coeffs each 
         */
        if (chipvariation == 1) {
            if (chiporder > 0)
                *ncoeffs +=
                    (1 + chiporder) * (2 + chiporder) * (nchip - 1) / 2;
            else
                *ncoeffs += (1 - chiporder) * (1 - chiporder) * (nchip - 1);
        }
        /*
         * no global crossterms 
         */
        *crossterm = 0;
        printf
            (" fitting polynomial surface of order %d and %d coefficients \n",
             order, *ncoeffs);
        if (chipvariation == 1)
            printf(" including chip-dependent terms of order %d \n",
                   chiporder);
    } else {
        order = abs(order);
        *ncoeffs = (1 + order) * (1 + order);
        /*
         * number of extra coeffs for local fits of 3 coeffs each 
         */
        if (chipvariation == 1) {
            if (chiporder > 0)
                *ncoeffs +=
                    (1 + chiporder) * (2 + chiporder) * (nchip - 1) / 2;
            else
                *ncoeffs += (1 - chiporder) * (1 - chiporder) * (nchip - 1);
        }
        /*
         * allow global crossterms 
         */
        *crossterm = 1;
        printf
            (" fitting polynomial surface of order %d with global cross-terms and %d coefficients \n",
             order, *ncoeffs);
        if (chipvariation == 1)
            printf(" including chip-dependent terms of order %d \n",
                   chiporder);
    }

    /*
    if (*ncoeffs >= 2000) {
        fprintf(stderr, " total number of coefficients is too high \n");
        exit(EXIT_FAILURE);
    }
    */
}
