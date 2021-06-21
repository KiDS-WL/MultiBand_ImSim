#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lensutil.h"  // contains a custom calloc routine, can be replaced if it's a pain

/* code to identify objects as groups of connected pixels within an image, 
   based on promam code varylpthreshold 
   LM 2008 

function return value: int, on success, number of separate objects detected. On failure, -1.

function arguments:

*image  1D input float array containing the 2D image to be segmented, where each
        element in the array is identified with respect to its integer x,y position
        as
          array_element = y*width + x
        noting C convention of all indices, starting at 0

*region 1D output int array containing the corresponding segmentation map. The pixel
        values in this map are 0 for pixels with no associated object detection,
        otherwise integer values, unique for each separate object that has been identified,
        in the range 1 - N, where N is the total number of objects detected, as given
        by the function return value. 

noise   input float value, the rms noise on the image (assumed the same for all pixels)

width   input int width of 2D image

height  input int height of 2D image

lower_area_limit  input int, number of pixels below which objects will be merged if other
                  connected pixels are found, or eliminated from output if not possible to merge. 
                  Typical value 8.

lower_merging_limit input int, maximum number of pixels which an object is allowed to have
                    if it is allowed to merge with another object, if fewer than this then objects
                    remain separated (deblended).  Typical value 15

arealimit         input int, maximum allowed size of any segmented object. Typical value 100000

fintensity_limit  input float, minimum intensity value used in the segmentation process, to prevent
                  segmentation going too deeply into the noise. Typical value, 2 x noise.

fmax_intensity    input float, maximum expected pixel intensity value. Typical value 0.25 times CCD saturation level in image units

imax_intensity    input int, maximum size of arrays needed in making linked lists - set to (int)(saturation-level/noise/2)

*/

int
varylpthresholdf2(float *image, int *region, float noise, int width,
                  int height, int lower_area_limit, int lower_merging_limit,
                  int arealimit, float fintensity_limit, float fmax_intensity,
                  int imax_intensity)
{

/*  this function is the main function for identifying the candidate
    opacities, using the method of varying thresholds to identify
    connected regions within the specified range of sizes */

    int i, x, y, xentry, centry, pvalue, intensity, nx, ny, touchingcolour,
        pixelcolour, ncolour, number, imagesize, xmin, xmax, ymin, ymax,
        ntouch, swap, i1, touchcolour[8], temp, intensity_limit,
        max_intensity;
    int *Istart, *Ilast, *Ilist;
    int *Cstart, *Clist;
    int *num, *oldnum, *mergenum;
    int *colour, *found, *last;
    int *thresh;
    float scale;

    if (noise > 0.) {
        scale = 2. / noise;
        max_intensity = (int) (fmax_intensity * scale);
        intensity_limit = (int) (fintensity_limit * scale);
    } else {
        printf(" noise value passed to varylpthresholdf < 0: value = %f \n",
               noise);
        printf(" test:  image[0] = %f \n", image[0]);
        printf(" test:  region[0] = %d \n", region[0]);
        printf(" test:  width = %d \n", width);
        printf(" test:  height = %d \n", height);
        return -1;
    }

    imagesize = width * height;

    Istart = ml_calloc(1 + imax_intensity, sizeof(int), NULL, "Istart");
    Ilast = ml_calloc(1 + imax_intensity, sizeof(int), NULL, "Ilast");
    Ilist = ml_calloc(imagesize, sizeof(int), NULL, "Ilist");
    colour = ml_calloc(imagesize, sizeof(int), NULL, "colour");
    num = ml_calloc(imagesize, sizeof(int), NULL, "num");
    oldnum = ml_calloc(imagesize, sizeof(int), NULL, "oldnum");
    mergenum = ml_calloc(imagesize, sizeof(int), NULL, "oldnum");
    found = ml_calloc(imagesize, sizeof(int), NULL, "found");
    Cstart = ml_calloc(imagesize, sizeof(int), NULL, "Cstart");
    Clist = ml_calloc(imagesize, sizeof(int), NULL, "Clist");
    last = ml_calloc(imagesize, sizeof(int), NULL, "last");
    thresh = ml_calloc(imagesize, sizeof(int), NULL, "thresh");


    ncolour = 0;

/*  set start values for linked list */
    for(i = 0; i < imagesize; i++) {
        Ilist[i] = -1;          /* links from one pixel to another */
        colour[i] = 0;
        found[i] = 0;
        num[i] = 0;
        oldnum[i] = 0;
        mergenum[i] = 0;
        Cstart[i] = 0;
        last[i] = 0;
        Clist[i] = 0;
        thresh[i] = 0;
        region[i] = 0;
    }


    for(i = 0; i < (1 + max_intensity); i++) {
        Istart[i] = -1;         /* first pixel in the linked list at this intensity */
        Ilast[i] = -1;          /* last pixel in the linked list at this intensity */
    }

/* work through every pixel, make linked list in intensity */

    for(y = 0; y < height; y++) {
        for(x = 0; x < width; x++) {
            xentry = y * width + x;
            pvalue = 0;
            if (image[xentry] > 0.)
                pvalue = (int) (image[xentry] * scale);
            if (pvalue > max_intensity)
                pvalue = max_intensity;
            if (Ilast[pvalue] < 0) {    /* if this intensity not found before, start new list */
                Istart[pvalue] = xentry;
                Ilast[pvalue] = xentry;
            } else {
                centry = Ilast[pvalue];
                Ilist[centry] = xentry;
                Ilast[pvalue] = xentry;
            }
        }
    }


/*  work through this list in order, starting with the brightest pixel.
    if a pixel is isolated at that intensity, assign it a colour.
    If it is has a coloured neighbour, assign it the colour of its
    neighbour.
    If it has more than one coloured neighbour, set all pixels of those
    colours to be the same colour.
    If the sum of coloured pixels for each colour will exceed the limit,
    stop the colouring for those colours.  */

/*  set start values for linked list */
    for(i = 0; i < imagesize; i++) {
        Clist[i] = -1;          /* links from one pixel to another */
    }

    for(i = 0; i < imagesize; i++) {
        Cstart[i] = -1;         /* first pixel in the linked list for this colour */
        last[i] = -1;           /* last pixel in the linked list for this colour */
    }



    for(intensity = max_intensity; intensity >= intensity_limit; intensity--) {
        if (Istart[intensity] > 0) {

            xentry = Istart[intensity];


            while (xentry > 0) {
                y = xentry / width;
                x = xentry - y * width;

/*  find whether there are any coloured adjoining pixels.  Make a list of
    unique colours */
                ntouch = 0;

                ymin = y - 1;
                if (ymin < 0)
                    ymin = 0;
                ymax = y + 1;
                if (ymax >= height)
                    ymax = height - 1;
                xmin = x - 1;
                if (xmin < 0)
                    xmin = 0;
                xmax = x + 1;
                if (xmax >= width)
                    xmax = width - 1;

                for(ny = ymin; ny <= ymax; ny++) {
                    for(nx = xmin; nx <= xmax; nx++) {
                        centry = ny * width + nx;
                        touchingcolour = colour[centry];
                        if (touchingcolour > 0 && found[touchingcolour] == 0) {
                            touchcolour[ntouch] = touchingcolour;
                            found[touchingcolour] = 1;
                            ntouch++;
                        }
                    }
                }

                /*
                 * reset "found" flags 
                 */
                for(i = 0; i < ntouch; i++)
                    found[touchcolour[i]] = 0;

/*  if the number of pixels touching is zero, then start a new colour for this pixel */
                if (ntouch == 0) {
                    ncolour++;
                    colour[xentry] = ncolour;
                    num[ncolour] = 1;
                    Cstart[ncolour] = xentry;   /* set first pixel for this colour in linked list */
                    Clist[xentry] = -1; /* start linked list for this colour */
                    last[ncolour] = xentry;     /* set last pixel in linked list */
                }

/*  if the number of previously-coloured pixels touching is just one, then
    the current pixel is simply coloured the same colour and the statistics updated */

/*  if the number of previously-coloured pixels touching is greater than one, then
    regions will have to be tested for merging properties and a decision
    taken about whether to merge them or not.  */

/*  The first part of this process is to order the touching regions by
    size, smallest first.  The current pixel will then be merged with these
    regions in that order, with a test at each stage on the size of the two
    regions being merged.  If BOTH of the two regions being considered
    at any one time are above the lower_size_limit then the whole merging
    process will be stopped and the regions will appear in the output as
    separate regions.  If one region is smaller than the lower_area_limit the
    merging will proceed */

                if (ntouch > 1) {
                    for(swap = 0; swap < (ntouch - 1); swap++) {
                        for(i = 0; i < (ntouch - 1); i++) {
                            i1 = i + 1;
                            if (num[touchcolour[i]] > num[touchcolour[i1]]) {
                                temp = touchcolour[i1];
                                touchcolour[i1] = touchcolour[i];
                                touchcolour[i] = temp;
                            }
                        }
                    }
                }

/* if the number of pixels touching is one or greater, then colour this pixel 
   the same as the first touching region (i.e. the smallest adjoining region) */
                if (ntouch >= 1) {
                    pixelcolour = touchcolour[0];       /* find the colour of the adjoining region */
                    colour[xentry] = pixelcolour;       /* set this pixel to that colour */
                    num[pixelcolour]++; /* count number of pixels in region */
                    Clist[last[pixelcolour]] = xentry;  /* add pixel into linked list */
                    last[pixelcolour] = xentry; /* reset last member in list */
                }


/*  if the number of different-coloured pixels touching is greater than one, then consider each in
    (size-ordered) turn, and set the counter mergenum at each merge.  Note that the merging still
    takes place, but what happens is that the array mergenum holds the size of the smallest region
    taking part.  That is later tested and if it is above the lower_area_limit then the
    regions are frozen as they were before the merging took place */

                if (ntouch > 1) {
                    for(i = 1; i < ntouch; i++) {
                        touchingcolour = touchcolour[i];

                        /*
                         * re-colour all pixels of the current colour (pixelcolour) to 
                         * the new colour (touchingcolour) 
                         */
                        centry = Cstart[pixelcolour];
                        colour[centry] = touchingcolour;
                        while (Clist[centry] >= 0) {
                            centry = Clist[centry];
                            colour[centry] = touchingcolour;
                        }
                        Clist[last[touchingcolour]] = Cstart[pixelcolour];      /* join lists */
                        last[touchingcolour] = last[pixelcolour];       /* reset last list member */
                        Cstart[pixelcolour] = -1;       /* disable linked list for original colour */

                        /*
                         * set mergenum to be the size of the smallest region taking part in the merging
                         * at this particular moment, but allow this value to increase between
                         * successive mergers.  Thus the value of mergenum should finally be the size
                         * of the second-largest region to have merged.  If this value is larger than the
                         * critical size, then the merger will be prevented once all pixels at
                         * this intensity have been considered 
                         */
                        if (num[pixelcolour] > mergenum[touchingcolour]) {
                            mergenum[touchingcolour] = num[pixelcolour];
                        }

                        /*
                         * add numbers to create statistics for the merged region 
                         */
                        num[touchingcolour] += num[pixelcolour];
                        num[pixelcolour] = 0;

                        /*
                         * update current pixel colour to the new colour 
                         */
                        pixelcolour = touchingcolour;
                    }
                }


/*  set the next pixel in the list at this intensity value */
                xentry = Ilist[xentry];
            }

/*  before moving to the next lowest intensity level, go through the image
    and update the list of valid coloured pixels (i.e. those lying
    within the specified size limits) */

/*  work through colours which have changed at this threshold level.  If each region satisfies 
    the selection criteria on size, merging, and surface brightness (i.e. surface brightness
    should not decrease), then update the map of regions */
            for(pixelcolour = 1; pixelcolour <= ncolour; pixelcolour++) {
                number = num[pixelcolour];
                if (number > oldnum[pixelcolour] &&     /* region size has changed */
                    number > lower_area_limit &&        /* region is big enough */
                    number < arealimit &&       /* region is not too big */
                    mergenum[pixelcolour] < lower_merging_limit) {      /* anti-merging criterion */
                    thresh[pixelcolour] = intensity;    /* store current threshold */
                    oldnum[pixelcolour] = number;       /* store the size of the region */


                    /*
                     * update pixel map of regions 
                     */
                    centry = Cstart[pixelcolour];
                    region[centry] = pixelcolour;
                    while (Clist[centry] >= 0) {
                        centry = Clist[centry];
                        region[centry] = pixelcolour;
                    }

                }
            }

/*  next lowest intensity level */
        }
    }


/*  go through map of final regions, identify which colours remain, and reallocate
    new colours so that there are no gaps in the colour table */
    for(centry = 0; centry < imagesize; centry++)
        colour[centry] = 0;
    ncolour = 0;

    for(centry = 0; centry < imagesize; centry++) {
        pixelcolour = region[centry];
        if (pixelcolour > 0 && colour[pixelcolour] == 0) {      /* test if this colour has already been found */
            ncolour++;          /* increment list of valid colours */
            colour[pixelcolour] = ncolour;      /* allocate a new colour to the old colour */
        }

        /*
         * set this pixel to have the new colour, and store parameters 
         */
        if (colour[pixelcolour] > 0) {
            region[centry] = colour[pixelcolour];
            num[colour[pixelcolour]] = oldnum[pixelcolour];
        }
    }

    ncolour++;
    free(Istart);
    free(Ilast);
    free(Ilist);
    free(colour);
    free(num);
    free(oldnum);
    free(mergenum);
    free(found);
    free(Cstart);
    free(Clist);
    free(last);
    free(thresh);

    /*
     * return the number of regions found 
     */
    return ncolour;
}
