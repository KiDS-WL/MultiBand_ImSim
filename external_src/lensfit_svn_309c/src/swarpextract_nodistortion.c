/*

  swarpextract

function to extract a postage stamp from an image
with correction for astrometric distortion

requires swarp wcs transformation information, 

uses swarp routines (E. Bertin)

LM 16 Jan 2009

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/* swarp include files */
#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "data.h"
#include "field.h"
#include "header.h"
#include "interpolate.h"
#include "prefs.h"
#include "projapprox.h"
#include "resample.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "weight.h"
#include "wcs/wcs.h"
#include "prefs.h"

#define pi M_PI

/* postage stamp extraction function with no distortion correction */

int
swarpextract_nodistortion(wcsstruct * wcs, float *apix, float *badpix,
                          int *dim, double *dpix, double *dbadpix, float objx,
                          float objy, float poserror, int pwidth, int pheight,
                          float noise, int *region,
                          double *scalefactor, float maxlimit)
{
    int ix, iy, x, y, xx, yy, sx, sy, nx, ny, spix, i, j, pixel;
    int ngoodpix, objectcolour;
    int halfwidth, halfheight;
    int allbad, nfound;
    int xcen, ycen;
    float negative_sigma_limit, dist, maxval;
    double rawpos[2], wcspos[2], wcsneg[2], wcscentre[2];
    double cosdist;

    for(y = 0; y < pheight; y++) {
        for(x = 0; x < pwidth; x++) {
            i = y * pwidth + x;
            dpix[i] = 0.;
            dbadpix[i] = 0.;
        }
    }
    ngoodpix = 0;

    ix = (int) (objx + 0.5) - 1;
    iy = (int) (objy + 0.5) - 1;
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

    nfound = 0;
    objectcolour = 0;
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
                    dist =
                        sqrt((x - halfwidth) * (x - halfwidth) +
                             (y - halfheight) * (y - halfheight));
                    if (dist <= poserror) {
                        if (region[j] > 0) {
                            if (region[j] != objectcolour) {
                                objectcolour = region[j];
                                nfound++;
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     * if more than one object exists inside poserror or any pixel exceeds the limit
     * then reject this stamp 
     */

    if (nfound > 1 || maxval > maxlimit) {
        ngoodpix = 0;
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                i = y * pwidth + x;
                dpix[i] = 0.;
                dbadpix[i] = 0.;
            }
        }
        // printf (" star rejected: close neighbour within distance %f \n",poserror);
        return (ngoodpix);
    }


    /*
     * get the wcs coord of the centre of the postage stamp 
     */

    rawpos[0] = (double) ix;
    rawpos[1] = (double) iy;

    raw_to_wcs(wcs, rawpos, wcscentre);

    if (*wcscentre == WCS_NOCOORD) {
        fprintf(stderr,
                " error calculating wcs position of postage stamp \n");
        fprintf(stderr, " image pixel position %d %d \n", ix, iy);
        exit(2);
    }

    /*
     * if it hasn't been set already, get the plate scale by differencing 
     * the wcs coords diagonally between two pixels near the field centre
     */
    if (scalefactor[0] <= 0.) {
        xcen = dim[0] / 2;
        ycen = dim[1] / 2;
        rawpos[0] = (double) xcen - 100;
        rawpos[1] = (double) ycen - 100;
        raw_to_wcs(wcs, rawpos, wcsneg);
        if (*wcsneg == WCS_NOCOORD) {
            fprintf(stderr,
                    " error calculating wcs position near field centre \n");
            fprintf(stderr, " image pixel position %lf %lf \n", rawpos[0],
                    rawpos[1]);
            exit(2);
        }
        rawpos[0] = (double) xcen + 100;
        rawpos[1] = (double) ycen + 100;
        raw_to_wcs(wcs, rawpos, wcspos);
        if (*wcspos == WCS_NOCOORD) {
            fprintf(stderr,
                    " error calculating wcs position near field centre \n");
            fprintf(stderr, " image pixel position %lf %lf \n", rawpos[0],
                    rawpos[1]);
            exit(2);
        }
        cosdist = sin(wcsneg[1] * pi / 180.) * sin(wcspos[1] * pi / 180.) +
            cos(wcsneg[1] * pi / 180.) * cos(wcspos[1] * pi / 180.) *
            cos((wcspos[0] - wcsneg[0]) * pi / 180.);
        scalefactor[0] = 200. * sqrt(2.) / acos(cosdist);       // units of radian/pixel
        printf(" scalefactor %lf arcsec/pixel\n",
               180. * 3600. / scalefactor[0] / pi);
    }

    /*
     * now read all the values into the dpix array 
     */

    for(y = 0; y < pheight; y++) {
        yy = y + iy - halfheight;
        if (yy >= 0 && yy < dim[1]) {
            sy = y + pheight / 2;
            if (sy >= pheight)
                sy -= pheight;
            for(x = 0; x < pwidth; x++) {
                xx = x + ix - halfwidth;
                if (xx >= 0 && xx < dim[0]) {
                    sx = x + pwidth / 2;
                    if (sx >= pwidth)
                        sx -= pwidth;
                    j = xx + yy * dim[0];
                    spix = sx + sy * pwidth;
                    dbadpix[spix] = badpix[j];
                    dpix[spix] = apix[j];
                    if (dbadpix[spix] == 0.)
                        dpix[spix] = 0.;
                    else
                        ngoodpix++;
                    /*
                     * check if this pixel or a near-neighbour belongs to another object 
                     */
                    for(j = -1; j <= 1; j++) {
                        for(i = -1; i <= 1; i++) {
                            nx = xx + i;
                            ny = yy + j;
                            if (nx < 0)
                                nx = 0;
                            if (nx >= dim[0])
                                nx = dim[0] - 1;
                            if (ny < 0)
                                ny = 0;
                            if (ny >= dim[1])
                                ny = dim[1] - 1;
                            pixel = nx + ny * dim[0];
                            if (region[pixel] > 0
                                && region[pixel] != objectcolour) {
                                dbadpix[spix] = 0.;     /* if pixel belongs to another object set its flag to bad */
                                dpix[spix] = 0.;        /* if pixel belongs to another object set its value to zero */
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     * find if there are any badpixels in the central region - if so
     * flag the entire frame as bad 
     */

    allbad = 0;
    for(y = 0; y < pheight; y++) {
        if (y < poserror || y > (pheight - poserror)) {
            for(x = 0; x < pwidth; x++) {
                if (x < poserror || x > (pwidth - poserror)) {
                    j = y * pwidth + x;
                    if (dbadpix[j] <= 0.)
                        allbad++;
                }
            }
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
        // printf (" star rejected, too many bad pixels: %d \n",allbad);
    }

    return (ngoodpix);

}
