/*

  swarpextract

function to extract a postage stamp from an image
with correction for astrometric distortion

requires swarp wcs transformation information, 

uses swarp routines (E. Bertin)

LM 16 Jan 2009

this version does not fit and subtract local background

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

/* postage stamp extraction function with no swarp distortion correction
   LM Dec 2009 - Sep 2010 */

void xytoradec (double *, double *, double, double *);

int tbswarpextract_distortion_measure (wcsstruct *wcs, 
		    float *apix, float *badpix, int *dim, double *dpix, double *dbadpix,
		    float objx, float objy, float poserror, int pwidth, int pheight,
		    float noise, float maxlimit, double *edist, double *wcentre)
{
 int ix, iy, x, y, xx, yy, sx, sy, spix, i, j, k, kk;
 int ngoodpix;
 int halfwidth, halfheight;
 int allbad;
 int xcen, ycen;
 float negative_sigma_limit, positive_sigma_limit;
 float maxval;
 double xcent, ycent;
 double cosdist[4];
 double wtest[8][2], theta;
 double rpos[8][2], cpos[2];
 double a,b,c,d,denom;
 double sxx,sxy,syy,sxdx,sxdy,sydx,sydy;
 double rpost[8][2], wpos[2], rcirc, scalefactor;

 for (y=0; y<pheight; y++)
   {
     for (x=0; x<pwidth; x++)
       {
	 i = y*pwidth + x;
	 dpix[i]=0.;
	 dbadpix[i]=0.;
       }
   }
 ngoodpix  = 0;

 ix = (int)(objx+0.5)-1;
 iy = (int)(objy+0.5)-1;
 halfwidth = pwidth/2;
 halfheight = pheight/2;

 if (noise <= 0.)
   {
     fflush(stdout);
     fprintf (stderr," error, noise = %g \n", noise);
     exit (2);
   }

 /* reject negative points from the data */
 negative_sigma_limit = -5.*noise;
 /* reject positive points if not associate with central object */
 positive_sigma_limit = 2.*noise;

 if (ix>=dim[0] || 
     ix<0 || 
     iy>=dim[1] || 
     iy<0)
  {
   printf(" object in catalogue is not on CCD x=%d y=%d\n",ix,iy);
   ngoodpix = 0;
   return (ngoodpix);
  }

 /* reject any postage stamp with an excessively high value anywhere on it (usually >half saturation level) */
 maxval = 0.;
 for (y=0; y<pheight; y++)
   {
     yy = y + iy - halfheight;
     if (yy>=0 && yy<dim[1])
       {
	 for (x=0; x<pwidth; x++)
	   {
	     xx = x + ix - halfwidth;
	     if (xx>=0 && xx<dim[0])
	       {
		 j = xx + yy*dim[0];
		 if (apix[j]>maxval) maxval = apix[j];
	       }
	   }
       }
   }

 /* if any pixel exceeds the limit
    then reject this stamp */

 if (maxval>maxlimit)
   {
     ngoodpix=0;
     for (y=0; y<pheight; y++)
       {
	 for (x=0; x<pwidth; x++)
	   {
	     i = y*pwidth + x;
	     dpix[i]=0.;
	     dbadpix[i]=0.;
	   }
       }
     //printf (" postage stamp rejected: maximum value %f exceeds saturation limit %f \n",maxval,maxlimit);
     return (ngoodpix);
   }

 /* get the pixel scale and the e1 and e2 shear at this location by testing a set
    of points around the location */
 xcen = ix;
 ycen = iy;

 rcirc = 20.;

 // shift the centre of the test point if its too close to the edge of the image
 if (xcen+rcirc >= dim[0]-10) xcen = dim[0]-rcirc-11;
 if (xcen-rcirc < 10) xcen = 11+rcirc;
 if (ycen+rcirc >= dim[1]-10) ycen = dim[1]-rcirc-11;
 if (ycen-rcirc < 10) ycen = 11+rcirc;

 // get the wcs of this test position
 cpos[0]=(double)xcen;
 cpos[1]=(double)ycen;
 raw_to_wcs(wcs, cpos, wcentre);
 if (*wcentre == WCS_NOCOORD)
     {
       fprintf(stderr," error calculating wcs position of postage stamp \n");
       fprintf(stderr," image pixel position %d %d \n",ix,iy);
       exit(2);
     }

 // lay down points in a circle and measure separation between them
 for (k=0; k<8; k++)
   {
     theta = k*2.*pi/8.;
     rpos[k][0]=(double)xcen + rcirc*sin(theta);
     rpos[k][1]=(double)ycen + rcirc*cos(theta);
     raw_to_wcs(wcs, rpos[k], wtest[k]);
       if (*wtest[k] == WCS_NOCOORD)
	 {
	   fprintf(stderr," error calculating wcs position near field centre \n");
	   fprintf(stderr," image pixel position %lf %lf \n",rpos[k][0],rpos[k][1]);
	   exit(2);
	 }
   }
 // find angles between opposing pairs of points
 scalefactor=0.;
 for (k=0; k<4; k++)
   {
     kk = k+4;
     cosdist[k] = sin(wtest[k][1]*pi/180.)*sin(wtest[kk][1]*pi/180.) + 
       cos(wtest[k][1]*pi/180.)*cos(wtest[kk][1]*pi/180.)*cos((wtest[k][0]-wtest[kk][0])*pi/180.);
     cosdist[k] = acos(cosdist[k]); // convert to separation in radians
   }
 // get approximate scalefactor from geometric mean of two separations
 scalefactor = sqrt(4.*rcirc*rcirc/(cosdist[0]*cosdist[2])); // units of pixels/radian
 //printf(" approx scalefactor %f arcsec/pixel \n",180.*3600./pi/scalefactor);
 // convert celestial coordinates to tangent plane and accumulate summations for full transformation
 sxx=sxy=syy=sxdx=sxdy=sydx=sydy=0.;
 xcent=ycent=0.;
 for (k=0; k<8; k++)
   {
     // adjust positions relative to centre of test circle
     rpos[k][0] -= xcen;
     rpos[k][1] -= ycen;
     // create celestial positions from ideal tangent plane
     xytoradec(rpos[k],wcentre,scalefactor,wpos); 
     // convert wcs back to pixel on image
     wcs_to_raw(wcs, wpos, rpost[k]); 
     if (*rpost[k] == WCS_NOCOORD)
       {
	 fprintf(stderr," error converting coords in swarpextract \n");
	 fprintf(stderr," pixel wcs  = %lf %lf \n",wpos[0],wpos[1]);
	 fprintf(stderr," centre wcs = %lf %lf \n",wcentre[0],wcentre[1]);
	 fprintf(stderr," pixel %lf %lf \n",rpost[k][0],rpost[k][1]);
	 exit(2);
       }
     xcent += rpost[k][0];
     ycent += rpost[k][1];
   }
 xcent /= 8.;
 ycent /= 8.;
 for (k=0; k<8; k++)
   {
     // convert back to relative positions
     rpost[k][0] -= xcent;
     rpost[k][1] -= ycent;
     // linear transform summations
     sxx += rpos[k][0]*rpos[k][0]/rcirc/rcirc;
     sxy += rpos[k][0]*rpos[k][1]/rcirc/rcirc;
     syy += rpos[k][1]*rpos[k][1]/rcirc/rcirc;
     sxdx += rpos[k][0]*rpost[k][0]/rcirc/rcirc;
     sxdy += rpos[k][0]*rpost[k][1]/rcirc/rcirc;
     sydx += rpos[k][1]*rpost[k][0]/rcirc/rcirc;
     sydy += rpos[k][1]*rpost[k][1]/rcirc/rcirc;    
     //printf("sums %f %f %f %f %f %f %f \n",sxx,sxy,syy,sxdx,sxdy,sydx,sydy);   
   }
 // hence get transformation coefficients
 denom = sxx*syy-sxy*sxy;
 if (denom != 0.)
   {
     a = (syy*sxdx - sxy*sydx)/denom;
     b = (sxx*sydx - sxy*sxdx)/denom;
     c = (syy*sxdy - sxy*sydy)/denom;
     d = (sxx*sydy - sxy*sxdy)/denom;
     //printf("sums %f %f %f %f %f \n",a,b,c,d,denom);
   }
 else
   {
     fflush(stdout);
     fprintf(stderr," denom = 0 \n");
     fprintf(stderr," %f %f %f \n",sxx,syy,sxy);
     exit(2);
   }
 edist[0]=a*scalefactor;
 edist[1]=b*scalefactor;
 edist[2]=c*scalefactor;
 edist[3]=d*scalefactor;
 /* testing only from here on
 // derive accurate scalefactor
 scf = sqrt((a+d)*(a+d)+(c-b)*(c-b))/2.;  
 // remove scalefactor from transformation
 a /= scf;
 b /= scf;
 c /= scf;
 d /= scf;
 scalefactor *= scf;
 //scalefactor = 1./scalefactor;  // convert to radian/pixel
 //scalefactor *= 180.*3600./pi;   // convert to arcsec/pixel
 //printf(" scale %f %f \n",scf,180.*3600./scalefactor/pi);
 // rotation
 costheta = (a+d)/2.;
 sintheta = (c-b)/2.;
 //printf(" cos, sin %f %f \n",costheta,sintheta);
 // rotation-corrected shear
 ap = a*costheta - b*sintheta;
 bp = a*sintheta + b*costheta;
 cp = c*costheta - d*sintheta;
 dp = c*sintheta + d*costheta;
 //printf(" shear %f %f %f %f \n",1.-ap,dp-1.,-bp,-cp);
 // final shears
 ap = (ap-dp)/2.;
 bp = (bp+cp)/2.;
 a = cosdist[0]/cosdist[2];
 a = (a-1.)/(a+1.);
 b = cosdist[3]/cosdist[1];
 b = (b-1.)/(b+1.);
 printf(" comparison %f %f %f %f \n",ap,bp,a,b);
 // check transformation
 xcent=ycent=0.;
 for (k=0; k<8; k++)
   {
     theta = k*2.*pi/8.;
     rpos[k][0]= rcirc*sin(theta);
     rpos[k][1]= rcirc*cos(theta);
     xytoradec(rpos[k],wcentre,scalefactor,wpos); 
     // convert wcs back to pixel on image
     wcs_to_raw(wcs, wpos, rpost[k]); 
     if (*rpost[k] == WCS_NOCOORD)
       {
	 fprintf(stderr," error converting coords in swarpextract \n");
	 fprintf(stderr," pixel wcs  = %lf %lf \n",wpos[0],wpos[1]);
	 fprintf(stderr," centre wcs = %lf %lf \n",wcentre[0],wcentre[1]);
	 fprintf(stderr," pixel %d %d \n",rpost[k][0],rpost[k][1]);
	 exit(2);
       }
     xcent += rpost[k][0];
     ycent += rpost[k][1];
   }
 xcent /= 8.;
 ycent /= 8.;
 for (k=0; k<8; k++)
   {
     // convert back to relative positions
     rpost[k][0] -= xcent;
     rpost[k][1] -= ycent;
     a = rpos[k][0]*(1.+e1dist[0]) +  rpos[k][1]*e2dist[0];
     b = rpos[k][0]*e2dist[0]      +  rpos[k][1]*(1.-e1dist[0]);
     wpos[0] = a*costheta - b*sintheta;
     wpos[1] = a*sintheta + b*costheta;
     // printf(" position %d %f %f %f %f %f %f \n",k,wpos[0],wpos[1],rpost[k][0],rpost[k][1],wpos[0]-rpost[k][0],wpos[1]-rpost[k][1]);
   }
 */

 // printf(" distortion %lf %lf %f %f %f \n",wcentre[0],wcentre[1],scalefactor,e1dist[0],e2dist[0]);

 for (y=0; y<pheight; y++)
   {
     for (x=0; x<pwidth; x++)
       {
	 spix = x + y*pwidth;
	 dbadpix[spix] = 0.;
	 dpix[spix] = 0.;	 
       }
   }

 // read main image values into postage stamp, swapping quadrants
 for (y=0; y<pheight; y++)
   {
     yy = y + iy - halfheight;
     if (yy>=0 && yy<dim[1])
       {
	 sy = y+pheight/2;
	 if (sy>=pheight) sy-=pheight;
	 for (x=0; x<pwidth; x++)
	   {
	     xx = x + ix - halfwidth;
	     if (xx>=0 && xx<dim[0])
	       {
		 sx = x+pwidth/2;
		 if (sx>=pwidth) sx -= pwidth;
		 j = xx + yy*dim[0];
		 spix = sx+sy*pwidth;
		 dbadpix[spix] = badpix[j];
		 dpix[spix] = apix[j];
		 if (dbadpix[spix]==0.) 
		     dpix[spix]=0.;
		 else
		     ngoodpix++;
	       }
	   }
       }
   }

 // return if all pixels have now been flagged as bad
 if (ngoodpix <= 0) 
   {
     ngoodpix=0;
     for (y=0; y<pheight; y++)
       {
	 for (x=0; x<pwidth; x++)
	   {
	     i = y*pwidth + x;
	     dpix[i]=0.;
	     dbadpix[i]=0.;
	   }
       }
     //printf (" star %f %f rejected: no good pixels remaining \n",objx,objy);
     return(ngoodpix);
   }

 /* find if there are any badpixels in the central region - if so
     flag the entire frame as bad */

  allbad = 0;
  for (y=0; y<pheight; y++)
    {
      if ( y<poserror || y>(pheight-poserror) )
	{
	  for (x=0; x<pwidth; x++)
	    {
	      if ( x<poserror || x>(pwidth-poserror) )
		{
		  j = y*pwidth+x;
		  if (dbadpix[j]<=0.) allbad++;
		}
	    }
	}
    }
	    
  if (allbad > 0 )
    {
      for (y=0; y<pheight; y++)
	{
	  for (x=0; x<pwidth; x++)
	    {
	      j = y*pwidth+x;
	      dpix[j]=0.;
	      dbadpix[j]=0.;
	    }
	}
      ngoodpix = 0;
      //printf (" star %f %f rejected, too many bad pixels: %d \n",objx,objy,allbad);
      return (ngoodpix);
    }

  /* run through good pixels and accumulate values for background fitting */

  /*

  nfit=0;
  hsize=(double)pwidth/2.;

  for (y=0; y<pheight; y++)
    {
      sy = y + pheight/2;
      if (sy>=pheight) sy -= pheight;
      yy = sy + iy - halfheight;
      for (x=0; x<pwidth; x++)
	{
	  j = y*pwidth+x;
	  // quadrant swap
	  sx = x + pwidth/2;
	  if (sx>=pwidth) sx -= pwidth;
	  xx = sx + ix - halfwidth;
	  if (yy>=0 && yy<dim[1] && xx>=0 && xx<dim[0])
	    {
	      jpix = xx + yy*dim[0];
	      if (dbadpix[j]>0. && region[jpix]==0)
		{
		  yd = (double)y;
		  if (y>=pheight/2) yd = (double)(y-pheight);
		  xd = (double)x;
		  if (x>=pwidth/2) xd = (double)(x-pwidth);
		  xfit[nfit]=xd/hsize;
		  yfit[nfit]=yd/hsize;
		  wfit[nfit]=1.;
		  zfit[nfit]=dpix[j];
		  nfit++;
		}
	    }
	}
    }

  if (nfit<10) 
    {
      ngoodpix=0;
      for (y=0; y<pheight; y++)
	{
	  for (x=0; x<pwidth; x++)
	    {
	      j = y*pwidth+x;
	      dpix[j]=0.;
	      dbadpix[j]=0.;
	    }
	}
      return ngoodpix;
    }

  // fit linear surface 
  order=1;
  ncoeffs=3;
  crossterm=0;

  for (i=0; i<=ncoeffs; i++)
       {
	 w[i] = 0.;
	 for (j=0; j<=ncoeffs; j++)
	      {
		u[i][j] = 0.;
		v[i][j] = 0.;
	      }
       }

  svdfit2dsqc(xfit,yfit,zfit,wfit,nfit,order,crossterm,avals,u,v,w);

  // subtract background from good pixels 

  nfit=0;
  for (y=0; y<pheight; y++)
    {
      for (x=0; x<pwidth; x++)
	{
	  j = y*pwidth+x;
	  if (dbadpix[j]>0.)
	    {
	      yd = (double)y;
	      if (y>=pheight/2) yd = (double)(y-pheight);
	      xd = (double)x;
	      if (x>=pwidth/2) xd = (double)(x-pwidth);
	      xfit[nfit] = xd/hsize;
	      yfit[nfit] = yd/hsize;
	      zfit[nfit] = avals[1]+avals[2]*yfit[nfit]+avals[3]*xfit[nfit];
	      dpix[j] -= zfit[nfit];
	    }
	}
    }

  */

  // check for negative pixels and flag any found
  for (y=0; y<pheight; y++)
    {
      for (x=0; x<pwidth; x++)
	{
	  j = y*pwidth+x;
	  if (dbadpix[j]>0. && dpix[j]<negative_sigma_limit)
	    {
	      dpix[j]=0.;
	      dbadpix[j]=0.;
	      ngoodpix--;
	    }
	}
    }

 return (ngoodpix);

}

