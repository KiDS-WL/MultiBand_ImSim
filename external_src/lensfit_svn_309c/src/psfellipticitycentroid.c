#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

void
psfellipticitycentroid(double *psf, int pheight, int pwidth, double *psfe,
                       double *centroid)
{
    double Q1, Q2, Q11, Q12, Q22, Q_norm, denom, xd, yd;
    double radius, dist, distsq;
    int x, y, x1, y1, pixel;

    Q1 = 0.;
    Q2 = 0.;
    Q11 = 0.;
    Q12 = 0.;
    Q22 = 0.;
    Q_norm = 0.;

    radius = 8.;
    distsq = radius * radius;

    /*
     * measure first moments 
     */

    for(y = 0; y < pheight / 2; y++) {
        for(x = 0; x < pwidth / 2; x++) {
            pixel = x + y * pwidth;
            dist = (double) (x * x + y * y);
            if (dist <= distsq) {
                Q1 += psf[pixel] * x;
                Q2 += psf[pixel] * y;
                Q_norm += psf[pixel];
            }
        }
        for(x = pwidth / 2; x < pwidth; x++) {
            pixel = x + y * pwidth;
            x1 = x - pwidth;
            dist = (double) (x1 * x1 + y * y);
            if (dist <= distsq) {
                Q1 += psf[pixel] * x1;
                Q2 += psf[pixel] * y;
                Q_norm += psf[pixel];
            }
        }
    }

    for(y = pheight / 2; y < pheight; y++) {
        y1 = y - pheight;
        for(x = 0; x < pwidth / 2; x++) {
            pixel = x + y * pwidth;
            dist = (double) (x * x + y1 * y1);
            if (dist <= distsq) {
                Q1 += psf[pixel] * x;
                Q2 += psf[pixel] * y1;
                Q_norm += psf[pixel];
            }
        }
        for(x = pwidth / 2; x < pwidth; x++) {
            pixel = x + y * pwidth;
            x1 = x - pwidth;
            dist = (double) (x1 * x1 + y1 * y1);
            if (dist <= distsq) {
                Q1 += psf[pixel] * x1;
                Q2 += psf[pixel] * y1;
                Q_norm += psf[pixel];
            }
        }
    }

    psfe[0] = 0.;
    psfe[1] = 0.;

    if (Q_norm > 0.) {
        Q1 = Q1 / Q_norm;
        Q2 = Q2 / Q_norm;

        centroid[0] = Q1;
        centroid[1] = Q2;

        /*
         * measure second moments 
         */

        Q_norm = 0.;

        for(y = 0; y < pheight / 2; y++) {
            yd = y - Q2;
            for(x = 0; x < pwidth / 2; x++) {
                pixel = x + y * pwidth;
                xd = x - Q1;
                dist = xd * xd + yd * yd;
                if (dist <= distsq) {
                    Q11 += psf[pixel] * xd * xd;
                    Q12 += psf[pixel] * xd * yd;
                    Q22 += psf[pixel] * yd * yd;
                    Q_norm += psf[pixel];
                }
            }
            for(x = pwidth / 2; x < pwidth; x++) {
                pixel = x + y * pwidth;
                x1 = x - pwidth;
                xd = x1 - Q1;
                dist = xd * xd + yd * yd;
                if (dist <= distsq) {
                    Q11 += psf[pixel] * xd * xd;
                    Q12 += psf[pixel] * xd * yd;
                    Q22 += psf[pixel] * yd * yd;
                    Q_norm += psf[pixel];
                }
            }
        }

        for(y = pheight / 2; y < pheight; y++) {
            y1 = y - pheight;
            yd = y1 - Q2;
            for(x = 0; x < pwidth / 2; x++) {
                pixel = x + y * pwidth;
                xd = x - Q1;
                dist = xd * xd + yd * yd;
                if (dist <= distsq) {
                    Q11 += psf[pixel] * xd * xd;
                    Q12 += psf[pixel] * xd * yd;
                    Q22 += psf[pixel] * yd * yd;
                    Q_norm += psf[pixel];
                }
            }
            for(x = pwidth / 2; x < pwidth; x++) {
                pixel = x + y * pwidth;
                x1 = x - pwidth;
                xd = x1 - Q1;
                dist = xd * xd + yd * yd;
                if (dist <= distsq) {
                    Q11 += psf[pixel] * xd * xd;
                    Q12 += psf[pixel] * xd * yd;
                    Q22 += psf[pixel] * yd * yd;
                    Q_norm += psf[pixel];
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

        if (psfe[0] > 1.)
            psfe[0] = 1.;
        if (psfe[1] > 1.)
            psfe[1] = 1.;
        if (psfe[0] < -1.)
            psfe[0] = -1.;
        if (psfe[1] < -1.)
            psfe[1] = -1.;

    }

}
