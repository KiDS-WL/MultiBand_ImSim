#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <nrutil.h>

static int
compare(const void *ii, const void *jj)
{
    float *i, *j;
    i = (float *) ii;
    j = (float *) jj;
    if (*j > *i)
        return -1;
    if (*j < *i)
        return 1;
    return 0;
}


void
mediansub(float *apix, float *badpix, int *dim, float maxlimit,
          float *sortarray)
{
    int x, y, n, i, pixel;
    float medval;

    /*
     * make sorted index of pixel values 
     */

    n = 0;

    for(y = 0; y < dim[1]; y++) {
        for(x = 0; x < dim[0]; x++) {
            i = y * dim[0] + x;
            if (badpix[i] > 0. && apix[i] < maxlimit) {
                sortarray[n] = apix[i];
                n++;
            }
        }
    }

    if (n > 0) {
        qsort(sortarray, n, sizeof(float), compare);
        medval = sortarray[n / 2];
    } else {
        fprintf(stderr, " no valid pixels in image \n");
        exit(2);
    }

    /*
     * subtract median background value 
     */

    for(y = 0; y < dim[1]; y++) {
        for(x = 0; x < dim[0]; x++) {
            pixel = y * dim[0] + x;
            if (badpix[pixel] > 0.)
                apix[pixel] -= medval;
        }
    }

    /*
     * char filename[80];
     * fitsfile *tempptr;
     * strcpy(filename,"backsub.fits");
     * remove(filename);               
     * int status;
     * status = 0;       
     * fits_create_file(&tempptr, filename, &status);
     * int bitpix;
     * bitpix = -32;
     * int naxis;
     * long naxes[2];
     * naxes[0]=dim[0];
     * naxes[1]=dim[1];
     * naxis = 2;
     * fits_create_img(tempptr,  bitpix, naxis, naxes, &status);
     * long fpixel, nelements;
     * fpixel=1;
     * nelements = dim[0]*dim[1];
     * fits_write_img(tempptr, TFLOAT, fpixel, nelements, apix, &status);
     * fits_close_file(tempptr, &status);
     */

    return;

}
