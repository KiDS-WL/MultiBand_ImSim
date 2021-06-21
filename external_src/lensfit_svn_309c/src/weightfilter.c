/* function based on weightimagefilter to remove potential spurious
   bad pixel masking.  Use in conjunction with star positions to decide
   if a bad pixel flag is likely to be correct or not */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

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


int weightfilter(float *apix, float *opix, float *changed, float *array, int *dim, int fheight)
{
  int more, nc, niter;
  int y1;
  int x, y, dy, pixel, qpixel, n, i;
  float val;

  nc=0;

  /* run a vertically-oriented median erosion filter of height 1+ 2*fheight and
     width 1 over the pixel data, and any bad pixels that get eroded are reflagged */

  for (y=0; y<dim[1]; y++)
    {
      for (x=0; x<dim[0]; x++)
	{
	  pixel = x + y*dim[0];
	  opix[pixel]=apix[pixel];
	  if (apix[pixel]<=0.)
	    {
	      n=0;
	      for (dy=-fheight; dy<=fheight; dy++)
		{
		  y1 = y+dy;
		  if (y1>=0 && y1<dim[1])
		    {
		      qpixel = x + y1*dim[0];
		      array[n] = apix[qpixel];
		      n++;
		    }
		}
	      qsort (array, n, sizeof(float), compare); 
	      opix[pixel] = array[(n/2)];
	      if (opix[pixel]>0.)
		{
		  changed[nc]=pixel;
		  nc++;
		}
	    }
	}
    }

  /* go through and reattach any deleted pixels that contact bad ones */

  more=1;
  niter=0;
  while (more==1)
    {
      more=0;
      for (i=0; i<nc; i++)
	{
	  pixel=changed[i];
	  val=opix[pixel];
	  y=pixel/dim[0];
	  x=pixel-y*dim[0];
	  if(x>0)
	    if(opix[pixel-1]<=0.) opix[pixel]=apix[pixel];
	  if(x<dim[0]-1)
	    if(opix[pixel+1]<=0.) opix[pixel]=apix[pixel];
	  if(y>0)
	    if(opix[pixel-dim[0]]<=0.) opix[pixel]=apix[pixel];
	  if(y<dim[1]-1)
	    if(opix[pixel+dim[0]]<=0.) opix[pixel]=apix[pixel];
	  if (opix[pixel]!=val) 
	    {
	      more=1;
	      niter++;
	    }
	}
    }

  nc -= niter;
  return(nc);

}
