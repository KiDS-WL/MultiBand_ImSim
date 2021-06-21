#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
extractpostagestamp(float *apix, float *badpix, int *dim, double *dpix,
                    double *dbadpix, float *temp, float *badtemp, float objx,
                    float objy, float poserror, float ximageoffset,
                    float yimageoffset, int pwidth, int pheight, float noise,
                    float maxlimit)
{
    int ix, iy, x, y, xx, yy, i, j, jj;
    int ngoodpix;
    int halfwidth, halfheight;
    int allbad;
    float negative_sigma_limit, maxval;
    double hxsize, hysize, xval, yval;

    for(y = 0; y < pheight; y++) {
        for(x = 0; x < pwidth; x++) {
            i = y * pwidth + x;
            temp[i] = 0.;
            dpix[i] = 0.;
            badtemp[i] = 0.;
        }
    }
    ngoodpix = 0;

    ix = (int) (objx + ximageoffset + 0.5) - 1;
    iy = (int) (objy + yimageoffset + 0.5) - 1;
    halfwidth = pwidth / 2;
    halfheight = pheight / 2;

    if (noise <= 0.) {
        printf(" error, noise = %g \n", noise);
        exit(1);
    }

    /*
     * reject negative points from the data 
     */
    negative_sigma_limit = -5. * noise;

    if (ix >= dim[0] || ix < 0 || iy >= dim[1] || iy < 0) {
        printf(" object in catalogue is not on CCD x=%d y=%d\n", ix, iy);
        ngoodpix = 0;
        return (ngoodpix);
    }



    /*
     * search for this object in the supplied object list 
     */

    maxval = 0.;

    for(y = 0; y < pheight; y++) {
        yy = y + iy - halfheight;
        if (yy >= 0 && yy < dim[1]) {
            for(x = 0; x < pwidth; x++) {
                xx = x + ix - halfwidth;
                if (xx >= 0 && xx < dim[0]) {
                    j = xx + yy * dim[0];
                    if (apix[j] > maxval)
                        maxval = apix[j];
                }
            }
        }
    }

    /*
     * if  any pixel exceeds the limit
     * then reject this stamp 
     */

    if (maxval > maxlimit) {
        ngoodpix = 0;
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                i = y * pwidth + x;
                dpix[i] = 0.;
                dbadpix[i] = 0.;
                temp[i] = 0.;
                badtemp[i] = 0.;
            }
        }
        return (ngoodpix);
    }


    /*
     * now read all the values into the temp array 
     */

    for(y = 0; y < pheight; y++) {
        yy = y + iy - halfheight;
        if (yy >= 0 && yy < dim[1]) {
            for(x = 0; x < pwidth; x++) {
                xx = x + ix - halfwidth;
                if (xx >= 0 && xx < dim[0]) {
                    j = xx + yy * dim[0];
                    i = y * pwidth + x;
                    temp[i] = apix[j];
                    if (temp[i] > negative_sigma_limit && badpix[j] > 0.) {
                        badtemp[i] = badpix[j];
                    } else {
                        badtemp[i] = 0.;
                    }
		    /*
                    if (region[j] > 0 && region[j] != objectcolour) {
		    badtemp[i] = 0.;        // if pixel belongs to another object set its flag to bad 
		    temp[i] = 0.;   // if pixel belongs to another object set its value to zero 
                    }
		    */
                }
            }
        }
    }


    /*
     * find if there are any badpixels in the central region - if so
     * flag the entire frame as bad 
     */

    allbad = 0;
    for(y = halfheight - poserror; y <= halfheight + poserror; y++) {
        for(x = halfwidth - poserror; x <= halfwidth + poserror; x++) {
            j = y * pwidth + x;
            if (badtemp[j] <= 0.)
                allbad++;
        }
    }

    if (allbad > 0) {
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                j = y * pwidth + x;
                dpix[j] = 0.;
                dbadpix[j] = 0.;
            }
        }
        ngoodpix = 0;
        return (ngoodpix);
    } else {
        hxsize = pwidth / 2;
        hysize = pheight / 2;
        ngoodpix = 0;
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                jj = y * pwidth + x;
                if (badtemp[jj] > 0.) {
                    xval = (double) x / hxsize - 1.;
                    yval = (double) y / hysize - 1.;
                    ngoodpix++;
                    /*
                     * if ( (xval*xval+yval*yval) <= 1.)
                     * {
                     * ngoodpix++;
                     * }
                     * else
                     * {
                     * temp[jj] = 0.;
                     * badtemp[jj] = 0.;
                     * }
                     */
                }
                /*
                 * else
                 * {
                 * temp[jj] = 0.;
                 * } 
                 */
            }
        }
    }


    /*
     * swap quadrants 
     */

    for(y = 0; y < halfheight; y++) {
        yy = y + halfheight;
        for(x = 0; x < halfwidth; x++) {
            i = y * pwidth + x;
            xx = x + halfwidth;
            j = yy * pwidth + xx;
            dpix[j] = (double) temp[i];
            dbadpix[j] = (double) badtemp[i];
        }
    }

    for(y = halfheight; y < pheight; y++) {
        yy = y - halfheight;
        for(x = 0; x < halfwidth; x++) {
            i = y * pwidth + x;
            xx = x + halfwidth;
            j = yy * pwidth + xx;
            dpix[j] = (double) temp[i];
            dbadpix[j] = (double) badtemp[i];
        }
    }

    for(y = 0; y < halfheight; y++) {
        yy = y + halfheight;
        for(x = halfwidth; x < pwidth; x++) {
            i = y * pwidth + x;
            xx = x - halfwidth;
            j = yy * pwidth + xx;
            dpix[j] = (double) temp[i];
            dbadpix[j] = (double) badtemp[i];
        }
    }

    for(y = halfheight; y < pheight; y++) {
        yy = y - halfheight;
        for(x = halfwidth; x < pwidth; x++) {
            i = y * pwidth + x;
            xx = x - halfwidth;
            j = yy * pwidth + xx;
            dpix[j] = (double) temp[i];
            dbadpix[j] = (double) badtemp[i];
        }
    }

    return (ngoodpix);

}
