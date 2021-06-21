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
#define VERBOSE 0

/* static kernel function */
static void make_kernel(double, double *, interpenum);

void




svdfit2dsqc(double x[], double y[], double z[], double weight[], int nobj,
            int order, int crossterm, double a[], double **u, double **v,
            double w[]);


/* tangent plane transformation function */
void
xytoradec(double *rpos, double *wcentre, double scalefactor, double *wpos)
{
    double ax, ay, denom, sindecentre, cosdecentre, relra, ra, dec,
        cosdist, racentre, decentre;

/* function to convert X and Y tangent-plane coordinates into
   RA and DEC, given the coords
   of the tangent-plane projection point, and the plate scale factor.
   Input coords are assumed in degrees, but plate scale in radians/pixel.
   Based on old fortran function, LM. */

/* convert position to radians */
    racentre = wcentre[0] * pi / 180.;
    decentre = wcentre[1] * pi / 180.;

    sindecentre = sin(decentre);
    cosdecentre = cos(decentre);

/*   take account of the plate scale assuming it has been input in pixels/radian */

    if (scalefactor > 0.) {
        ax = rpos[0] / scalefactor;
        ay = rpos[1] / scalefactor;
    } else {
        fflush(stdout);
        fprintf(stderr,
                " error in value of astrometric scale factor in function xytoradec \n");
        exit(2);
    }

/*    calculate relative RA */

    denom = (cosdecentre - ay * sindecentre);

    if (denom != 0.) {
        relra = atan(ax / denom);
    } else {
        if (ax == 0.)
            relra = 0.;         // dec=90, ra undefined
        else
            relra = ax > 0. ? pi / 2. : -pi / 2.;
    }

    //    calculate dec allowing for difficult cases
    if (ax == 0. && ay == 0.)
        dec = decentre;
    else {
        if (sindecentre == 0.)
            dec = atan(ay * cos(relra));
        else {
            if (sindecentre == 1.)
                dec =
                    ax != 0. ? atan(sin(relra) / ax) : atan(-cos(relra) / ay);
            else {
                if (sindecentre == -1.)
                    dec =
                        ax !=
                        0. ? atan(-sin(relra) / ax) : atan(-cos(relra) / ay);
                else {
                    if (denom != 0.)
                        dec =
                            atan((ay * cosdecentre +
                                  sindecentre) * cos(relra) / denom);
                    else if (cos(relra) != 0.)
                        dec = decentre > 0. ? pi / 2. : -pi / 2.;
                    else if (ax != 0.)
                        dec = atan((1. + ay * ay) * sindecentre / ax);
                    else
                        dec = decentre > 0. ? pi / 2. : -pi / 2.;
                }
            }
        }
    }


    /*
     * now there is a 180-degree ambiguity in the above solution,
     * in that the calculated position can be diametrically opposite
     * its true position on the celestial sphere.  The easiest way
     * to test for this is to measure the angle between the calculated
     * position and the tangent-plane point, and if that is greater
     * than 90 degrees then assume that the 180-degree ambiguity has occurred.
     * This is reasonable since we wouldn't expect to project points onto
     * a tangent plane from the reverse side of the celestial sphere. 
     */

    cosdist = sindecentre * sin(dec) + cosdecentre * cos(dec) * cos(relra);
    if (cosdist < 0.) {
        dec = -dec;
        if (relra != 0.)
            relra += pi;
        if (relra > pi)
            relra -= 2. * pi;
    }
    //    convert relative RA into actual RA
    ra = racentre - relra;
    if (ra > 2. * pi)
        ra -= 2. * pi;
    if (ra < 0.)
        ra += 2. * pi;

    // convert back to degrees

    wpos[0] = ra * 180. / pi;
    wpos[1] = dec * 180. / pi;

    return;

}


/* postage stamp extraction function. 
   this version doesn't flag background objects as this is now done later in lensfit.
   Pixels of any object are ignored for fitting the background however.
   LM Dec 2009
*/


int
tbswarpextract_noflag(wcsstruct * wcs,
                      float *apix, float *badpix, int *dim, double *dpix,
                      double *dbadpix, float objx, float objy,
                      float poserror, int pwidth, int pheight,
                      ikernelstruct * ikernel, float noise,
                      double *scalefactor, float maxlimit, double *rpos,
                      double *wpos, double *wneg, double *wcentre,
                      double **kernel)
//double *xfit, double *yfit,
//double *zfit, double *wfit, double **u, double **v)
{
    int ix, iy, x, y, xx, yy, sx, sy, spix, i, j, n, pixel;
    int ngoodpix;
    int halfwidth, halfheight, ival[2];
    int allbad;
    int naxis = 2;
    int kwidth[2], xpix0, ypix0, xcen, ycen;
    float negative_sigma_limit, positive_sigma_limit;
    float maxval;
    double cosdist, val, swarpshift[2];

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
        fflush(stdout);
        fprintf(stderr, " error, noise = %g \n", noise);
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
        //printf(" object in catalogue is not on CCD x=%d y=%d\n",ix,iy);
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
     * if more than one object exists inside poserror or any pixel exceeds the limit
     * then reject this stamp 
     */

    if (maxval > maxlimit) {
        ngoodpix = 0;
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                i = y * pwidth + x;
                dpix[i] = 0.;
                dbadpix[i] = 0.;
            }
        }
        // printf (" star rejected: saturation limit %f exceeded \n",maxlimit);
        return (ngoodpix);
    }


    /*
     * get the wcs coord of the centre of the postage stamp 
     */

    rpos[0] = (double) ix;
    rpos[1] = (double) iy;

    raw_to_wcs(wcs, rpos, wcentre);

    if (*wcentre == WCS_NOCOORD) {
        fflush(stdout);
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
            fflush(stdout);
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
            fflush(stdout);
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

    if (scalefactor[0] < (180. * 3600. / pi)
        || scalefactor[0] > (180. * 3600. * 10. / pi)) {
        fflush(stdout);
        fprintf(stderr,
                " error in scalefactor in swarpextract, out of expected range \n");
        fprintf(stderr, " scalefactor %lf arcsec/pixel\n",
                180. * 3600. / scalefactor[0] / pi);
        fprintf(stderr, " likely mismatch with PSF creation \n");
        exit(2);
    }

    /*
     * now read all the values into the dpix array 
     */

    for(y = 0; y < pheight; y++) {
        for(x = 0; x < pwidth; x++) {
            /*
             * pixel location in postagestamp including quadrant swap 
             */
            sx = x - pwidth / 2;
            sy = y - pheight / 2;
            if (sx < 0)
                sx += pwidth;
            if (sy < 0)
                sy += pheight;
            spix = sx + sy * pwidth;
            /*
             * transform coordinates to allow for astrometric distortion 
             */
            rpos[0] = 1. + (double) (x - pwidth / 2);   // pixel location wrt centre
            rpos[1] = 1. + (double) (y - pheight / 2);  // pixel location wrt centre
            /*
             * this is the key step.  The current x,y locations are
             * considered to be points in a tangent plane projection
             * centred on the current postage stamp and with the same
             * pixel scale as the image, locally.  The pixel values
             * will be interpolated onto this grid using the swarp
             * interpolation.  The grid is set so that no
             * interpolation should happen in the centre
             */
            xytoradec(rpos, wcentre, scalefactor[0], wpos);     // convert tangentplane xy to wcs
            wcs_to_raw(wcs, wpos, rpos);        // convert wcs back to pixel on image
            if (*rpos == WCS_NOCOORD) {
                fprintf(stderr,
                        " error converting coords in swarpextract \n");
                fprintf(stderr, " pixel wcs  = %lf %lf \n", wpos[0], wpos[1]);
                fprintf(stderr, " centre wcs = %lf %lf \n", wcentre[0],
                        wcentre[1]);
                fprintf(stderr, " pixel %d %d \n", x, y);
                exit(2);
            }
            /*
             * work out interpolation kernel 
             */
            for(n = 0; n < naxis; n++) {
                val = rpos[n];
                ival[n] =
                    (ikernel->interptype[n] ==
                     INTERP_NEARESTNEIGHBOUR) ? (int) (val -
                                                       0.50001) : (int) val;
                swarpshift[n] = val - ival[n];
                make_kernel(swarpshift[n], kernel[n], ikernel->interptype[n]);
                kwidth[n] = ikernel->width[0];
            }
            xpix0 = kwidth[0] / 2;
            ypix0 = kwidth[1] / 2;
            /*
             * create postagestamp pixel value by swarp interpolation 
             */
            dpix[spix] = 0.;
            dbadpix[spix] = 0.;
            /*
             * check the pixel is on the actual image 
             */
            if (ival[0] >= 0 && ival[0] < dim[0] && ival[1] >= 0
                && ival[1] < dim[1]) {
                /*
                 * set the badpix value 
                 */
                pixel = ival[0] + ival[1] * dim[0];
                dbadpix[spix] = badpix[pixel];
                for(j = 0; j < kwidth[1]; j++) {
                    for(i = 0; i < kwidth[0]; i++) {
                        xx = ival[0] + i - xpix0;       // ival[] holds the value of the nearest image pixel
                        /*
                         * dont allow the kernel to spill off the edge of the actual image 
                         */
                        if (xx < 0)
                            xx = 0;
                        if (xx >= dim[0])
                            xx = dim[0] - 1;
                        yy = ival[1] + j - ypix0;
                        if (yy < 0)
                            yy = 0;
                        if (yy >= dim[1])
                            yy = dim[1] - 1;
                        pixel = xx + yy * dim[0];
                        /*
                         * set the postagestamp pixel to bad if ANY neighbour 
                         * being used in the interpolation is bad 
                         */
                        if (badpix[pixel] > 0.) {
                            dpix[spix] +=
                                apix[pixel] * kernel[0][i] * kernel[1][j];
                        } else {
                            dbadpix[spix] = 0.;
                        }
                    }
                }
            }
            if (dbadpix[spix] == 0.) {
                dpix[spix] = 0.;
            } else {
                // count numbers of good pixels
                ngoodpix++;
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
        //printf(" no good pixels \n"); fflush(stdout);
        return (ngoodpix);
    }

    /*
     * find if there are any badpixels in the central region - if so
     * flag the entire frame as bad 
     */

    allbad = 0;
    for(y = 0; y < pheight; y++) {
        if (y < poserror || y > (pheight - poserror)) {
            for(x = 0; x < pwidth; x++) {
                j = y * pwidth + x;
                if (x < poserror || x > (pwidth - poserror)) {
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
        //printf (" star rejected, too many bad pixels: %d \n",allbad); fflush(stdout);
        return (ngoodpix);
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
        //printf (" nfit = %d \n",nfit);
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


    // subtract background from good pixels 

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



/****** make_kernel **********************************************************
PROTO	void make_kernel(double pos, double *kernel, interpenum interptype)
PURPOSE	Conpute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/01/2008
 ***/
void
make_kernel(double pos, double *kernel, interpenum interptype)
{
    double x, val, sinx1, sinx2, sinx3, cosx1;

    if (interptype == INTERP_NEARESTNEIGHBOUR)
        *kernel = 1;
    else if (interptype == INTERP_BILINEAR) {
        *(kernel++) = 1.0 - pos;
        *kernel = pos;
    } else if (interptype == INTERP_LANCZOS2) {
        if (pos < 1e-5 && pos > -1e5) {
            *(kernel++) = 0.0;
            *(kernel++) = 1.0;
            *(kernel++) = 0.0;
            *kernel = 0.0;
        } else {
            x = -PI / 2.0 * (pos + 1.0);
#ifdef HAVE_SINCOS
            sincos(x, &sinx1, &cosx1);
#else
            sinx1 = sin(x);
            cosx1 = cos(x);
#endif
            val = (*(kernel++) = sinx1 / (x * x));
            x += PI / 2.0;
            val += (*(kernel++) = -cosx1 / (x * x));
            x += PI / 2.0;
            val += (*(kernel++) = -sinx1 / (x * x));
            x += PI / 2.0;
            val += (*kernel = cosx1 / (x * x));
            val = 1.0 / val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *kernel *= val;
        }
    } else if (interptype == INTERP_LANCZOS3) {
        if (pos < 1e-5 && pos > -1e5) {
            *(kernel++) = 0.0;
            *(kernel++) = 0.0;
            *(kernel++) = 1.0;
            *(kernel++) = 0.0;
            *(kernel++) = 0.0;
            *kernel = 0.0;
        } else {
            x = -PI / 3.0 * (pos + 2.0);
#ifdef HAVE_SINCOS
            sincos(x, &sinx1, &cosx1);
#else
            sinx1 = sin(x);
            cosx1 = cos(x);
#endif
            val = (*(kernel++) = sinx1 / (x * x));
            x += PI / 3.0;
            val += (*(kernel++) =
                    (sinx2 = -0.5 * sinx1 - 0.866025403785 * cosx1)
                    / (x * x));
            x += PI / 3.0;
            val += (*(kernel++) =
                    (sinx3 = -0.5 * sinx1 + 0.866025403785 * cosx1)
                    / (x * x));
            x += PI / 3.0;
            val += (*(kernel++) = sinx1 / (x * x));
            x += PI / 3.0;
            val += (*(kernel++) = sinx2 / (x * x));
            x += PI / 3.0;
            val += (*kernel = sinx3 / (x * x));
            val = 1.0 / val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *kernel *= val;
        }
    } else if (interptype == INTERP_LANCZOS4) {
        if (pos < 1e-5 && pos > -1e5) {
            *(kernel++) = 0.0;
            *(kernel++) = 0.0;
            *(kernel++) = 0.0;
            *(kernel++) = 1.0;
            *(kernel++) = 0.0;
            *(kernel++) = 0.0;
            *(kernel++) = 0.0;
            *kernel = 0.0;
        } else {
            x = -PI / 4.0 * (pos + 3.0);
#ifdef HAVE_SINCOS
            sincos(x, &sinx1, &cosx1);
#else
            sinx1 = sin(x);
            cosx1 = cos(x);
#endif
            val = (*(kernel++) = sinx1 / (x * x));
            x += PI / 4.0;
            val += (*(kernel++) = -(sinx2 = 0.707106781186 * (sinx1 + cosx1))
                    / (x * x));
            x += PI / 4.0;
            val += (*(kernel++) = cosx1 / (x * x));
            x += PI / 4.0;
            val += (*(kernel++) =
                    -(sinx3 = 0.707106781186 * (cosx1 - sinx1)) / (x * x));
            x += PI / 4.0;
            val += (*(kernel++) = -sinx1 / (x * x));
            x += PI / 4.0;
            val += (*(kernel++) = sinx2 / (x * x));
            x += PI / 4.0;
            val += (*(kernel++) = -cosx1 / (x * x));
            x += PI / 4.0;
            val += (*kernel = sinx3 / (x * x));
            val = 1.0 / val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *(kernel--) *= val;
            *kernel *= val;
        }
    } else
        error(EXIT_FAILURE,
              "*Internal Error*: Unknown interpolation type in ",
              "make_kernel()");

    return;
}





/****** load_field ************************************************************
PROTO	fieldstruct *load_field(catstruct *cat, int frameno, int fieldno)
PURPOSE	Initialize a field structure (in read mode)
INPUT	Cat structure,
	FITS extension number in file (0=primary)
	Field number in coaddition (for various config settings)
	Field flags,
OUTPUT	The new field pointer if OK, NULL otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/06/2007
 ***/
void
load_dfield(fieldstruct * field, catstruct * cat, char *headfilename,
            int frameno, int fieldno)
{
    tabstruct *tab, *intab;
    //   fieldstruct     *field;
    char *pstr;
    int i;

/* First allocate memory for the new field (and nullify pointers) 
  QCALLOC(field,fieldstruct, 1); */
    field->flags = FIELD_READ;

/* Set conversion flags */
    field->cflags = (prefs.interp_flag[fieldno] ? CONVERT_INTERP : 0)
        | (prefs.subback_flag[fieldno] ? CONVERT_BACKSUB : 0);
    field->frameno = frameno;
    field->fieldno = fieldno;
    field->cat = new_cat(1);
    strcpy(field->cat->filename, cat->filename);
    strcpy(field->filename, cat->filename);

    intab = cat->tab;
    for(i = frameno; i--;)
        intab = intab->nexttab;
    copy_tab_fromptr(intab, field->cat, 0);

    tab = field->tab = field->cat->tab;
    tab->cat = field->cat;

/* A short, "relative" version of the filename */
    if (!(field->rfilename = strrchr(field->filename, '/')))
        field->rfilename = field->filename;
    else
        field->rfilename++;

    sprintf(gstr, "Looking for %s", field->rfilename);

/* Create a file name with a "header" extension */
    strcpy(field->hfilename, headfilename);
    if (!(pstr = strrchr(field->hfilename, '.')))
        pstr = field->hfilename + strlen(field->hfilename);
    sprintf(pstr, "%s", prefs.head_suffix);

    // printf(" header file %s \n",field->hfilename);

    NFPRINTF(OUTPUT, gstr);


/* Insert additional header informations from the "header" file */
    field->headflag = !read_aschead(field->hfilename, frameno, tab);

    if (tab->naxis < 1)
        error(EXIT_FAILURE, "*Error*: Zero-dimensional table in ",
              field->filename);

/* Force data to be at least 2D */
    if (tab->naxis < 2) {
        tab->naxis = 2;
        QREALLOC(tab->naxisn, int, 2);
        tab->naxisn[1] = 1;
    }

/* Some defaults */
    field->fascale = 1.0;
    field->fscale = prefs.fscale_default[fieldno];
    field->gain = prefs.gain_default[fieldno];
    field->saturation = prefs.sat_default[fieldno];

/* Force input celestial system to "PIXEL" if requested by user */
    if (prefs.celsys_type == CELSYS_PIXEL)
        for(i = 0; i < tab->naxis; i++) {
            sprintf(gstr, "CTYPE%-3d", i + 1);
            fitswrite(tab->headbuf, gstr, "PIXEL", H_STRING, T_STRING);
        }

/* Read WCS information in FITS header */
    field->wcs = read_wcs(tab);

/* Read additional field-related information in FITS header */
    readfitsinfo_field(field, tab);

/* Set field width and field height (the latter can be "virtual") */
    field->width = tab->naxisn[0];
    field->height = 1;
    for(i = 1; i < tab->naxis; i++)
        field->height *= tab->naxisn[i];
    field->npix = field->width * field->height;

/*-- Background */
    field->backw = prefs.back_size[fieldno] < field->width ?
        prefs.back_size[fieldno]
        : field->width;
    field->backh = prefs.back_size[fieldno] < field->height ?
        prefs.back_size[fieldno]
        : field->height;
    field->nbackp = field->backw * field->backh;
    if (field->backw > 0) {
        if ((field->nbackx = (field->width - 1) / field->backw + 1) < 1)
            field->nbackx = 1;
        if ((field->nbacky = (field->height - 1) / field->backh + 1) < 1)
            field->nbacky = 1;
    } else {
        field->nbackx = 1;
        field->nbacky = 1;
    }
    field->nback = field->nbackx * field->nbacky;
    field->nbackfx = field->nbackx > 1 ? prefs.back_fsize[fieldno] : 1;
    field->nbackfy = field->nbacky > 1 ? prefs.back_fsize[fieldno] : 1;
/* Set the back_type flag if absolute background is selected */
    field->back_type = prefs.back_type[fieldno];
    field->backdefault = prefs.back_default[fieldno];

/* Check flux scale */
    if (field->fscale <= 0.0) {
        // warning(field->filename, " has flux scale = 0: I will take 1 instead");
        field->fscale = 1.0;
    }
    //  printf(" flux scale factor %lf \n",field->fscale);

    return;
}
