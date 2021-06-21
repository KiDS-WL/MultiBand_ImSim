/*

  swarpextract

function to extract a postage stamp from an image
with correction for astrometric distortion

requires swarp wcs transformation information, 

uses swarp routines (E. Bertin)

LM 16 Jan 2009

arguments rationalised, background subtraction removed and background object
flagging allowed if specified  LM Sep 2010

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
#define VERBOSE 0

int
tbswarpextract_nodistortion(wcsstruct * wcs,
                            float *apix, float *badpix, float *opix, int *dim,
                            double *dpix, double *dbadpix,
                            float objx, float objy, float poserror,
                            int pwidth, int pheight,
                            float noise, int *region, double *scalefactor,
                            float maxlimit, double *rpos, double *wpos,
                            double *wneg, double *wcentre, int flagbad)
{
    int ix, iy, x, y, xx, yy, sx, sy, spix, i, j;
    int ngoodpix, objectcolour;
    int halfwidth, halfheight;
    int allbad, nfound;
    int xcen, ycen;
    float negative_sigma_limit, positive_sigma_limit;
    float dist, maxval;
    double cosdist, rsq;

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
    /*
     * reject positive points if not associate with central object 
     */
    positive_sigma_limit = 2. * noise;

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

    rpos[0] = (double) ix;
    rpos[1] = (double) iy;

    raw_to_wcs(wcs, rpos, wcentre);

    if (*wcentre == WCS_NOCOORD) {
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
        rpos[0] = (double) xcen - 100;
        rpos[1] = (double) ycen - 100;
        raw_to_wcs(wcs, rpos, wneg);
        if (*wneg == WCS_NOCOORD) {
            fprintf(stderr,
                    " error calculating wcs position near field centre \n");
            fprintf(stderr, " image pixel position %lf %lf \n", rpos[0],
                    rpos[1]);
            exit(2);
        }
        rpos[0] = (double) xcen + 100;
        rpos[1] = (double) ycen + 100;
        raw_to_wcs(wcs, rpos, wpos);
        if (*wpos == WCS_NOCOORD) {
            fprintf(stderr,
                    " error calculating wcs position near field centre \n");
            fprintf(stderr, " image pixel position %lf %lf \n", rpos[0],
                    rpos[1]);
            exit(2);
        }
        cosdist = sin(wneg[1] * pi / 180.) * sin(wpos[1] * pi / 180.) +
            cos(wneg[1] * pi / 180.) * cos(wpos[1] * pi / 180.) *
            cos((wpos[0] - wneg[0]) * pi / 180.);
        scalefactor[0] = 200. * sqrt(2.) / acos(cosdist);       // units of radian/pixel
        if (VERBOSE == 1)
            printf(" scalefactor %lf arcsec/pixel\n",
                   180. * 3600. / scalefactor[0] / pi);
    }

    for (y = 0; y < pheight; y++) {
      for (x = 0; x < pwidth; x++) {
	spix = x + y * pwidth;
	dbadpix[spix] = 0.;
	dpix[spix] = 0.;
      }}

    // read main image values into postage stamp, swapping quadrants
    // and flag pixels as bad if they contain another object
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
		    // measure distance from nominal star position and either flag
		    // using the corrected badpixel mask inside poserror or the nominal
		    // pixel mask outside
		    rsq = (y-halfheight)*(y-halfheight)+(x-halfwidth)*(x-halfwidth);
		    if (rsq > poserror*poserror)
		      dbadpix[spix] = badpix[j];
		    else
		      dbadpix[spix] = opix[j];
                    dpix[spix] = apix[j];
		    // flag out background objects if specified
		    if (flagbad==1 && region[j] > 0 && region[j] != objectcolour)
		      dbadpix[spix] = 0.;
		    // if flagged, set pixel value to zero
                    if (dbadpix[spix] == 0.)
                        dpix[spix] = 0.;
                    else
                        ngoodpix++;
                }
            }
        }
    }

    // return if all pixels have now been flagged as bad
    if (ngoodpix <= 0) {
        ngoodpix = 0;
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                i = y * pwidth + x;
                dpix[i] = 0.;
                dbadpix[i] = 0.;
            }
        }
        return (ngoodpix);
    }

    /*
     * find if there are any badpixels in the central region - if so
     * flag the entire frame as bad 
     */

    allbad = 0;
    for(y = 0; y < pheight; y++) 
      {
	yy = y>=halfheight ? y-pheight : y;
	for(x = 0; x < pwidth; x++) 
	  {
	    xx = x>=halfwidth ? x-pwidth : x;
	    rsq = xx*xx + yy*yy;
	    if (rsq < poserror*poserror)
	      {
		j = y * pwidth + x;
		if (dbadpix[j] <= 0.)
		  allbad++;
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
        return (ngoodpix);
        // printf (" star rejected, too many bad pixels: %d \n",allbad);
    }

    /*
     * run through good pixels and accumulate values for background fitting 
     */

    /*

    nfit = 0;
    hsize = (double) pwidth / 2.;

    for(y = 0; y < pheight; y++) {
        sy = y + pheight / 2;
        if (sy >= pheight)
            sy -= pheight;
        yy = sy + iy - halfheight;
        for(x = 0; x < pwidth; x++) {
            j = y * pwidth + x;
            // quadrant swap
            sx = x + pwidth / 2;
            if (sx >= pwidth)
                sx -= pwidth;
            xx = sx + ix - halfwidth;
            if (yy >= 0 && yy < dim[1] && xx >= 0 && xx < dim[0]) {
                jpix = xx + yy * dim[0];
                if (dbadpix[j] > 0. && region[jpix] == 0) {
                    yd = (double) y;
                    if (y >= pheight / 2)
                        yd = (double) (y - pheight);
                    xd = (double) x;
                    if (x >= pwidth / 2)
                        xd = (double) (x - pwidth);
                    xfit[nfit] = xd / hsize;
                    yfit[nfit] = yd / hsize;
                    wfit[nfit] = 1.;
                    zfit[nfit] = dpix[j];
                    nfit++;
                }
            }
        }
    }

    if (nfit < 10) {
        ngoodpix = 0;
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                j = y * pwidth + x;
                dpix[j] = 0.;
                dbadpix[j] = 0.;
            }
        }
        return ngoodpix;
    }

    // fit linear surface 

    order = 1;
    ncoeffs = 3;
    crossterm = 0;

    for(i = 0; i <= ncoeffs; i++) {
        w[i] = 0.;
        for(j = 0; j <= ncoeffs; j++) {
            u[i][j] = 0.;
            v[i][j] = 0.;
        }
    }

    svdfit2dsqc(xfit, yfit, zfit, wfit, nfit, order, crossterm, avals, u, v,
                w);

    
    //  subtract background from good pixels 

    nfit = 0;
    for(y = 0; y < pheight; y++) {
        for(x = 0; x < pwidth; x++) {
            j = y * pwidth + x;
            if (dbadpix[j] > 0.) {
                yd = (double) y;
                if (y >= pheight / 2)
                    yd = (double) (y - pheight);
                xd = (double) x;
                if (x >= pwidth / 2)
                    xd = (double) (x - pwidth);
                xfit[nfit] = xd / hsize;
                yfit[nfit] = yd / hsize;
                zfit[nfit] =
                    avals[1] + avals[2] * yfit[nfit] + avals[3] * xfit[nfit];
                dpix[j] -= zfit[nfit];
            }
        }
    }

    */

    // check for negative pixels and flag any found
    for(y = 0; y < pheight; y++) {
        for(x = 0; x < pwidth; x++) {
            j = y * pwidth + x;
            if (dbadpix[j] > 0. && dpix[j] < negative_sigma_limit) {
	      dpix[j] = 0.;
	      dbadpix[j] = 0.;
	      ngoodpix--;
            }
        }
    }

    return (ngoodpix);

}
