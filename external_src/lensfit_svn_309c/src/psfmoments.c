#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

// version corrected by Massimo Viola June 2015

void
psfmoments(double *psf, int pheight, int pwidth, double *psfe,
           double *centroid, double *moments)
{
    double Q1, Q2, Q11, Q12, Q22, Q_norm, denom;
    double dist, distsq, xxcen, yycen, xx, yy;
    double pval, sigma, sigmasq, rsq, w;
    int i, x, y, iter, niter;

    // initialise variables and set default return values
    Q1 = 0.;
    Q2 = 0.;
    Q11 = 0.;
    Q12 = 0.;
    Q22 = 0.;
    Q_norm = 0.;
    psfe[0] = 0.;
    psfe[1] = 0.;
    moments[0] = moments[1] = moments[2] = 0.;
    centroid[0] = centroid[1] = 0.;

    // sigma of Gaussian weight function for astrometry and moments
    sigma = 2.5; 
    sigmasq = sigma*sigma;
    niter = 3; // number of moments iterations for centroid
    // distance limit for moments summations
    dist = (double)pwidth/2; 
    distsq = dist*dist;

    // to match with the original lensfit tophat moments, set
    // niter=1, dist=8, sigma very large

    // iteratively find the weighted model centroid
    
    for (iter=0; iter<niter; iter++)
      {
        Q1 = Q2 = Q_norm = 0.; // initialise summations for this iteration
        for (y=0; y<pheight; y++)
          {
            // This what you have to change if you want to run on images with quadrants NOT swapped!
            yy = y>pheight/2 ? y-pheight : y; // allow for swapped quadrants
            //yy = y //If quadrants are not swapped !
            yycen = yy - centroid[1]; // centre the weight function on the current best guess centroid
            for (x=0; x<pwidth; x++)
              {
                i = y*pwidth + x;
                // This what you have to change if you want to run on images with quadrants NOT swapped!
                xx = x>pwidth/2 ? x-pwidth : x; // allow for swapped quadrants
                //xx = x //If quadrants are not swapped
                xxcen = xx - centroid[0]; // centre the weight function on the current best guess centroid                    
                dist = xxcen*xxcen + yycen*yycen;
                if (dist < distsq) // force outer boundary to be circular centred on zero
                  {
                    rsq = xxcen*xxcen + yycen*yycen; // radius squared defined in oversampled units relative to best-guess centroid
                    w = exp(-rsq/2./sigmasq); // weight function
                    pval = psf[i];
                    Q1 += xxcen*pval*w;
                    Q2 += yycen*pval*w;
                    Q_norm += pval*w;
                    //if (iter==3) {printf("%d\t%f\t%f\t%f\t\%f\t%f\n",i,xxcen,yycen,pval,Q1,Q2); }
                  }
              }
          }
        // printf("Q1 = %f\n", Q1);
        // printf("Q2 = %f\n", Q2);
        
        if (Q_norm>0.)
          {
            //The new centroid is the old one plus the shift computed in the last iteration!
            centroid[0] += Q1/Q_norm; 
            centroid[1] += Q2/Q_norm;
            // printf("xc = %f\n", centroid[0]);
            // printf("yc = %f\n", centroid[1]);
            // printf("Q_norm = %f\n", Q_norm);
            
          }
      }
    
    if (Q_norm > 0.) {
      /*
       * measure second moments 
       */
      Q_norm = 0.; // initialise summation
      for (y=0; y<pheight; y++)
        {
          // This what you have to change if you want to run on images with quadrants NOT swapped!
          yy = y>pheight/2 ? y-pheight : y; // allow for swapped quadrants
          //yy = y //If quadrants are not swapped
          yycen = yy - centroid[1]; // centre the weight function on the centroid
          for (x=0; x<pwidth; x++)
            {
              i = y*pwidth + x;
              // This what you have to change if you want to run on images with quadrants NOT swapped!
              xx = x>pwidth/2 ? x-pwidth : x; // allow for swapped quadrants
              //xx = x //If quadrants are not swapped
              
              xxcen = xx - centroid[0]; // centre the weight function on the centroid
              dist = xxcen*xxcen + yycen*yycen;
              
              if (dist < distsq) // force outer boundary to be circular centred on zero
                {
                  rsq = xxcen*xxcen + yycen*yycen; // radius squared relative to best-guess centroid
                  w = exp(-rsq/2./sigmasq); // weight function
                  pval = psf[i];
                  Q11 += pval * xxcen*xxcen * w;
                  Q12 += pval * xxcen*yycen * w;
                  Q22 += pval * yycen*yycen * w;
                  Q_norm += pval * w;
                  //printf("%d\t%f\t%f\t%f\t%f\t%f\n",i,xxcen,yycen,pval,Q11,Q22); 
                  
                }
            }
        }
      
      denom = Q11 * Q22 - Q12 * Q12;
      if (denom > 0.) {
        denom = Q11 + Q22 + 2. * sqrt(denom);
      } else {
        denom = Q11 + Q22;
      }
      if (denom != 0.) {
        psfe[0] = (Q11 - Q22) / denom;
        psfe[1] = (2. * Q12) / denom;
      }
      
      if (psfe[0] > 1.) psfe[0] = 1.;
      if (psfe[1] > 1.) psfe[1] = 1.;
      if (psfe[0] < -1.) psfe[0] = -1.;
      if (psfe[1] < -1.) psfe[1] = -1.;
      
      if (Q_norm > 0.) {
        moments[0] = Q11 / Q_norm;
        moments[1] = Q22 / Q_norm;
        moments[2] = Q12 / Q_norm;
      }
      
    }
    
    // printf("e1 = %f\n", psfe[0]);
    // printf("e2 = %f\n", psfe[1]);
    // printf("q11 = %f\n", moments[0]);
    // printf("q22 = %f\n", moments[1]);
    // printf("q12 = %f\n", moments[2]);
    // printf("Q_norm = %f\n", Q_norm);
    
}
