int
varylpthresholdf3(int nt, int imagenum)
{

/* code to identify objects as groups of connected pixels within an image, 
   based on promam code varylpthreshold 
   LM 2008 */

/*  this function is the main function for identifying the candidate
    opacities, using the method of varying thresholds to identify
    connected regions within the specified range of sizes */

    int i, x, y, xentry, centry, pvalue, intensity, nx, ny, touchingcolour,
        pixelcolour, ncolour, number, xmin, xmax, ymin, ymax,
        ntouch, swap, i1, touchcolour[8], temp, sintensity_limit;

    ncolour = 0;

    //printf(" thread %d image %d arraysize %d \n",nt,imagenum,nmaxregion[nt]);
    //fflush(stdout);

/*  set start values for linked list */
    for(i = 0; i < imagesize[nt]; i++) {
        iwork[nt][0][i] = 0;    // region
        iwork[nt][1][i] = 0;    // colour
        iwork[nt][2][i] = -1;   // links from one pixel to another 
        iwork[nt][3][i] = 0;    // Cstart
        iwork[nt][4][i] = 0;    // Clist
        iwork[nt][5][i] = 0;    // last
    }
    for(i = 0; i < nmaxregion[nt]; i++) {
        swork[nt][0][i] = 0;    // num
        swork[nt][1][i] = 0;    // oldnum
        swork[nt][2][i] = 0;    // mergenum
        swork[nt][3][i] = 0;    // found
    }

    //printf(" thread %d image %d arraysize %d \n",nt,imagenum,nmaxregion[nt]);
    //fflush(stdout);

    if (noise[imagenum] > 0.) {
        sintensity_limit = (int) (intensity_limit[nt] * scale[nt]);
    } else {
        printf(" noise value passed to varylpthresholdf < 0: value = %f \n",
               noise[imagenum]);
        printf(" test:  image[0] = %f \n", apix[nt][0]);
        printf(" test:  region[0] = %d \n", iwork[nt][0][0]);
        printf(" test:  width = %d \n", dim[imagenum][0]);
        printf(" test:  height = %d \n", dim[imagenum][1]);
        return -1;
    }

    for(i = 0; i < (1 + imax_intensity[nt]); i++) {
        Istart[nt][i] = -1;     /* first pixel in the linked list at this intensity */
        Ilast[nt][i] = -1;      /* last pixel in the linked list at this intensity */
    }

/* work through every pixel, make linked list in intensity */

    for(y = 0; y < dim[imagenum][1]; y++) {
        for(x = 0; x < dim[imagenum][0]; x++) {
            xentry = y * dim[imagenum][0] + x;
            pvalue = 0;
            if (apix[nt][xentry] > 0.)
                pvalue = (int) (apix[nt][xentry] * scale[nt]);
            if (pvalue > imax_intensity[nt])
                pvalue = imax_intensity[nt];
            if (Ilast[nt][pvalue] < 0) {        /* if this intensity not found before, start new list */
                Istart[nt][pvalue] = xentry;
                Ilast[nt][pvalue] = xentry;
            } else {
                centry = Ilast[nt][pvalue];
                iwork[nt][2][centry] = xentry;
                Ilast[nt][pvalue] = xentry;
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
    for(i = 0; i < imagesize[nt]; i++) {
        iwork[nt][4][i] = -1;   /* links from one pixel to another */
    }
    for(i = 0; i < imagesize[nt]; i++) {
        iwork[nt][3][i] = -1;   /* first pixel in the linked list for this colour */
        iwork[nt][5][i] = -1;   /* last pixel in the linked list for this colour */
    }


    for(intensity = imax_intensity[nt]; intensity >= sintensity_limit;
        intensity--) {
        if (Istart[nt][intensity] > 0) {

            xentry = Istart[nt][intensity];


            while (xentry > 0) {
                y = xentry / dim[imagenum][0];
                x = xentry - y * dim[imagenum][0];

/*  find whether there are any coloured adjoining pixels.  Make a list of
    unique colours */
                ntouch = 0;

                ymin = y - 1;
                if (ymin < 0)
                    ymin = 0;
                ymax = y + 1;
                if (ymax >= dim[imagenum][1])
                    ymax = dim[imagenum][1] - 1;
                xmin = x - 1;
                if (xmin < 0)
                    xmin = 0;
                xmax = x + 1;
                if (xmax >= dim[imagenum][0])
                    xmax = dim[imagenum][0] - 1;

                for(ny = ymin; ny <= ymax; ny++) {
                    for(nx = xmin; nx <= xmax; nx++) {
                        centry = ny * dim[imagenum][0] + nx;
                        touchingcolour = iwork[nt][1][centry];
                        if (touchingcolour > 0
                            && swork[nt][3][touchingcolour] == 0) {
                            touchcolour[ntouch] = touchingcolour;
                            swork[nt][3][touchingcolour] = 1;
                            ntouch++;
                        }
                    }
                }

                /*
                 * reset "found" flags 
                 */
                for(i = 0; i < ntouch; i++)
                    swork[nt][3][touchcolour[i]] = 0;

/*  if the number of pixels touching is zero, then start a new colour for this pixel */
                if (ntouch == 0) {
                    ncolour++;
                    iwork[nt][1][xentry] = ncolour;
                    swork[nt][0][ncolour] = 1;
                    iwork[nt][3][ncolour] = xentry;     /* set first pixel for this colour in linked list */
                    iwork[nt][4][xentry] = -1;  /* start linked list for this colour */
                    iwork[nt][5][ncolour] = xentry;     /* set last pixel in linked list */
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
                            if (swork[nt][0][touchcolour[i]] >
                                swork[nt][0][touchcolour[i1]]) {
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
                    iwork[nt][1][xentry] = pixelcolour; /* set this pixel to that colour */
                    if (swork[nt][0][pixelcolour] < 65535)
                        swork[nt][0][pixelcolour]++;    /* count number of pixels in region */
                    iwork[nt][4][iwork[nt][5][pixelcolour]] = xentry;   /* add pixel into linked list */
                    iwork[nt][5][pixelcolour] = xentry; /* reset last member in list */
                }


/*  if the number of different-coloured pixels touching is greater than one, then consider each in
    (size-ordered) turn, and set the counter mergenum at each merge.  Note that the merging still
    takes place, but what happens is that the array swork[nt][2] holds the size of the smallest region
    taking part.  That is later tested and if it is above the lower_area_limit then the
    regions are frozen as they were before the merging took place */

                if (ntouch > 1) {
                    for(i = 1; i < ntouch; i++) {
                        touchingcolour = touchcolour[i];

                        /*
                         * re-colour all pixels of the current colour (pixelcolour) to 
                         * the new colour (touchingcolour) 
                         */
                        centry = iwork[nt][3][pixelcolour];
                        iwork[nt][1][centry] = touchingcolour;
                        while (iwork[nt][4][centry] >= 0) {
                            centry = iwork[nt][4][centry];
                            iwork[nt][1][centry] = touchingcolour;
                        }
                        iwork[nt][4][iwork[nt][5][touchingcolour]] = iwork[nt][3][pixelcolour]; /* join lists */
                        iwork[nt][5][touchingcolour] = iwork[nt][5][pixelcolour];       /* reset last list member */
                        iwork[nt][3][pixelcolour] = -1; /* disable linked list for original colour */

                        /*
                         * set mergenum to be the size of the smallest region taking part in the merging
                         * at this particular moment, but allow this value to increase between
                         * successive mergers.  Thus the value of mergenum should finally be the size
                         * of the second-largest region to have merged.  If this value is larger than the
                         * critical size, then the merger will be prevented once all pixels at
                         * this intensity have been considered 
                         */
                        if (swork[nt][0][pixelcolour] >
                            swork[nt][2][touchingcolour]) {
                            swork[nt][2][touchingcolour] =
                                swork[nt][0][pixelcolour];
                        }

                        /*
                         * add numbers to create statistics for the merged region 
                         */
                        if ((int) swork[nt][0][touchingcolour] +
                            (int) swork[nt][0][pixelcolour] < 65535)
                            swork[nt][0][touchingcolour] +=
                                swork[nt][0][pixelcolour];
                        swork[nt][0][pixelcolour] = 0;

                        /*
                         * update current pixel colour to the new colour 
                         */
                        pixelcolour = touchingcolour;
                    }
                }


/*  set the next pixel in the list at this intensity value */
                xentry = iwork[nt][2][xentry];
            }

/*  before moving to the next lowest intensity level, go through the image
    and update the list of valid coloured pixels (i.e. those lying
    within the specified size limits) */

/*  work through colours which have changed at this threshold level.  If each region satisfies 
    the selection criteria on size, merging, and surface brightness (i.e. surface brightness
    should not decrease), then update the map of regions */
            for(pixelcolour = 1; pixelcolour <= ncolour; pixelcolour++) {
                number = swork[nt][0][pixelcolour];
                if (number > swork[nt][1][pixelcolour] &&       /* region size has changed */
                    number > lower_area_limit &&        /* region is big enough */
                    number < arealimit &&       /* region is not too big */
                    swork[nt][2][pixelcolour] < merging_area_limit) {   /* anti-merging criterion */
                    swork[nt][1][pixelcolour] = number; /* store the size of the region */


                    /*
                     * update pixel map of regions 
                     */
                    centry = iwork[nt][3][pixelcolour];
                    iwork[nt][0][centry] = pixelcolour;
                    while (iwork[nt][4][centry] >= 0) {
                        centry = iwork[nt][4][centry];
                        iwork[nt][0][centry] = pixelcolour;
                    }

                }
            }

/*  next lowest intensity level */
        }
    }


/*  go through map of final regions, identify which colours remain, and reallocate
    new colours so that there are no gaps in the colour table */
    for(centry = 0; centry < imagesize[nt]; centry++)
        iwork[nt][1][centry] = 0;
    ncolour = 0;

    for(centry = 0; centry < imagesize[nt]; centry++) {
        pixelcolour = iwork[nt][0][centry];
        if (pixelcolour > 0 && iwork[nt][1][pixelcolour] == 0) {        /* test if this colour has already been found */
            ncolour++;          /* increment list of valid colours */
            iwork[nt][1][pixelcolour] = ncolour;        /* allocate a new colour to the old colour */
        }

        /*
         * set this pixel to have the new colour, and store parameters 
         */
        if (iwork[nt][1][pixelcolour] > 0) {
            iwork[nt][0][centry] = iwork[nt][1][pixelcolour];
            swork[nt][0][iwork[nt][0][pixelcolour]] =
                swork[nt][1][pixelcolour];
        }
    }

    ncolour++;

/*  return the number of regions found */

    return ncolour;


}

