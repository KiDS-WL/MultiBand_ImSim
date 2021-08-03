/*
standalone code to read an input FITS image of a PSF and output as
a lensfit polynomial coefficients FITS file, assuming that the 
PSF should be invariant with position

Uses the lensfit function for outputting PSF coefficients (included here).

Compile with:

  cc -O3 -o psf_meimage2coeffs psf_meimage2coeffs.c -I/opt/local/include -L/opt/local/lib -lcfitsio

using -I, -L for wherever your cfitsio include and lib files are located

invoke with

  psf_meimage2coeffs <input me FITS file name> <output coeffs FITS file name>


LM 4/2/2014
updated for multiextension PSF file with coefficients per CCD  LM 26/7/2021

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>
#include <errno.h>

// define default pixel scale
#define DEFAULT_PIXEL_SCALE 0.214

long anaxes[2]; // input FITS image axis sizes
double *apix; 
double **spix;
float gs_scale[1];

int
writeglobalcoeffs(char *psfname, int correct_distortion, double **f, int ncoeffs, int order,
                  int chipvariation, int chiporder, double *scalefactor,
                  int pwidth, int pheight, int nxchip, int nychip,
                  int *numstar, int xchipsampling, int ychipsampling,
                  int big_gap, double hxsize, double hysize, float poserror)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[3], fpixel[3] = { 1, 1, 1 };
    int i, j, ii, bitpix;
    double *psf3d;
    int psfsize, psfsize3d, row;

    int tfields = 2;            /* number of table columns in table */
    char extname[] = "PSF star numbers";        /* extension name */
    char *ttype[] = { "CHIP", "NSTAR" };        /* names of columns */
    char *tform[] = { "1J", "1J" };     /* data types of columns */
    char *tunit[] = { "", "" }; /* units */

    /*
     * write out psf data as a fits file 
     */

    if (pwidth != pheight) {
        fprintf(stderr, " PSF subimages must be square! \n");
        return 1;
    }

    if (pwidth <= 0 || pheight <= 0) {
        fprintf(stderr, " error in subimage dimensions \n");
        return 1;
    }

    bitpix = -64;
    anaxis = 3;
    anaxes[0] = pwidth;
    anaxes[1] = pheight;
    anaxes[2] = ncoeffs;

    psfsize = anaxes[0] * anaxes[1];
    psfsize3d = ncoeffs * psfsize;
    psf3d = (double *) calloc(psfsize3d, sizeof(double));
    if (psf3d == NULL) {
        fprintf(stderr,
                " error allocating memory for 3D psf coeffs in writecoeffs_fits \n");
        fprintf(stderr, " memory requested = %d doubles \n", psfsize3d);
        exit(2);
    }

    fits_create_file(&afptr, psfname, &status);
    fits_create_img(afptr, bitpix, anaxis, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    for(i = 0; i < ncoeffs; i++) {
        for(j = 0; j < psfsize; j++) {
            ii = i * psfsize + j;
            psf3d[ii] = f[j][i];
        }
    }

    /*
     * write all PSF data into image array 
     */

    if (fits_write_pix(afptr, TDOUBLE, fpixel, psfsize3d, psf3d, &status)) {
        printf(" error writing psf coefficient data \n");
        exit(2);
    }

    /*
     * write keywords into header 
     */

    if (fits_update_key(afptr, TINT, "CORRECTD", &correct_distortion,
                        "distortion correction", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TDOUBLE, "PIXSCALE", scalefactor,
                        "pixel scalefactor", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "ORDER", &order,
                        "order of global fit", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "CHIPVAR", &chipvariation,
                        "chip-dependency allowed?", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_update_key(afptr, TINT, "CHORDER", &chiporder,
                        "chip fit order", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_update_key(afptr, TINT, "PWIDTH", &pwidth,
                        "postage stamp width", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "PHEIGHT", &pheight,
                        "postage stamp height", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NXCHIP", &nxchip,
                        "number of chips along x axis", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "NYCHIP", &nychip,
                        "number of chips along y axis", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "XSAMPL", &xchipsampling,
                        "x chip spacing", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "YSAMPL", &ychipsampling,
                        "y chip spacing", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "BIG_GAP", &big_gap,
                        "size of extra gap", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TDOUBLE, "HXSIZE", &hxsize,
                        "x size scaling", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TDOUBLE, "HYSIZE", &hysize,
                        "y size scaling", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TFLOAT, "TOL", &poserror,
                        "tolerance radius", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    fits_create_tbl(afptr, BINARY_TBL, nxchip * nychip, tfields, ttype, tform,
                    tunit, extname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < nxchip * nychip; i++) {
        row = i + 1;
        //printf(" chip %d numstars %d \n", row, numstar[i]);
        fits_write_col(afptr, TINT, 1, row, 1, 1, &row, &status);
        fits_write_col(afptr, TINT, 2, row, 1, 1, &numstar[i], &status);
    }
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }


    /*
     * close the fits file 
     */

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    free(psf3d);

    return (0);

}


int readimage(char *imagename)
{
  // read input image into array
  fitsfile *afptr; /* FITS file pointers */
  int j, next;
  int x, y, xx, yy;
  int pwidth, pheight;
  int status = 0;  /* CFITSIO status value must be initialized to zero */
  int anaxis;
  int hdunum, hdutype;
  long fpixel[2]={1,1};
  long nsize;
  char comment[FLEN_COMMENT];

/* open input image */
    fits_open_file(&afptr, imagename, READONLY, &status); 
    if (status) {
       fits_report_error(stderr, status); /* print error message */
       exit(EXIT_FAILURE);
    }

/* get number of extensions */
    fits_get_num_hdus(afptr, &hdunum, &status);
    if (status) {
       fits_report_error(stderr, status); /* print error message */
       exit(EXIT_FAILURE);
    }

    // remove primary hdu from hdu count
    if (hdunum>1)
      {
	hdunum -= 1;
      }

	// fits_movrel_hdu(afptr, 1, &hdutype, &status);
	// if (status) {
	//   fits_report_error(stderr, status); /* print error message */
	//   exit(EXIT_FAILURE);
	// }
	// if (hdutype != IMAGE_HDU)
	//   {
	//     fprintf(stderr,"extension does not image type \n");
	//     exit(EXIT_FAILURE);
	//   }    
  // try to get pixel scale keyword from primary header, set default 0.2 arcsec if not found
    char *key = "GS_SCALE";
    ffgky(afptr, TFLOAT, key, gs_scale, comment, &status);

    if (status) {
      // fprintf(stderr," failed to read keyword %s from fits file\n",key);
      printf("setting default pixel scale %f \n",(float)DEFAULT_PIXEL_SCALE);
      *gs_scale = (float)DEFAULT_PIXEL_SCALE;  // set default arcsec per pixel
      status=0;
    }

    // go through the HDUs and read the image data
    for (next=0; next<hdunum; next++)
      {
    
	// if more than one hdu, assume multi-extension image format
	// and move forward to the next extension
	if (hdunum>1)
	  {
	    printf("moving to hdu %d \n",next);
	    fits_movrel_hdu(afptr, 1, &hdutype, &status);
	    if (status) {
	      fits_report_error(stderr, status); /* print error message */
	      exit(EXIT_FAILURE);
	    }
	    if (hdutype != IMAGE_HDU)
	      {
		fprintf(stderr,"extension does not image type \n");
		exit(EXIT_FAILURE);
	      }
	  }

	/* read dimensions */
	fits_get_img_dim(afptr, &anaxis, &status);  
	fits_get_img_size(afptr, 2, anaxes, &status);
	if (status) {
	  fits_report_error(stderr, status); /* print error message */
	  exit(EXIT_FAILURE);
	}

	if (anaxis != 2) {
	  fprintf(stderr," readimage: number of axes = %d \n",anaxis);
	  printf("Error: images with other than 2 dimensions are not supported\n");
	  exit (2);
	}

	nsize = anaxes[0]*anaxes[1];
	pwidth = anaxes[0];
	pheight = anaxes[1];

	// allocate memory on first image read and assume all images have same size
	if (next==0)
	  {
	    apix = (double*)calloc(nsize, sizeof(double));
	    spix = (double**)calloc(nsize, sizeof(double*));
	    for (j=0; j<nsize; j++)
	      spix[j] = (double*)calloc(hdunum, sizeof(double));
	  }

	/* read input data into image array (convert to double) */
	if (fits_read_pix(afptr, TDOUBLE, fpixel, nsize, NULL, apix,
			  NULL, &status) )
	  {
	    fits_report_error(stderr, status); /* print error message */
	    fprintf(stderr," error reading pixel data \n");
	    exit (EXIT_FAILURE);
	  }

	// swap quadrants and store in spix array
	for (y=0; y<pheight; y++)
	  {
	    yy = y + pheight/2;
	    if (yy >= pheight) yy -= pheight;
	    for (x=0; x<pwidth; x++)
	      {
		xx = x + pwidth/2;
		if (xx >= pwidth) xx -= pwidth;
		spix[yy*pwidth + xx][next] = apix[y*pwidth+x];
	      }
	  }
	
      }
	
/* close image file */

    fits_close_file (afptr,&status);

    if (status) {
      fits_report_error(stderr, status); /* print error message */
      exit(EXIT_FAILURE);
    }

    // return number of extensions found
    return (hdunum);
    
}


/* opens a global PSF coefficients file and appends a table of shift coefficients */
int
updateglobalcoeffs(char *psfname, int ncoeffs, int order,
                   int chipvariation, int chiporder, double *xc, double *yc)
{
    fitsfile *afptr;
    int status = 0;
    int i;
    int row;
    int hdunum, hdutype;

    int tfields = 2;            /* number of table columns in table */
    char extname[] = "star shift coefficients"; /* extension name */
    char *ttype[] = { "XC", "YC" };     /* names of columns */
    char *tform[] = { "1D", "1D" };     /* data types of columns */
    char *tunit[] = { "", "" }; /* units */

    fits_open_file(&afptr, psfname, READWRITE, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * write keywords into header 
     */

    if (fits_update_key(afptr, TINT, "SORDER", &order,
                        "order of shifts global fit", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    if (fits_update_key(afptr, TINT, "SCHIPVAR", &chipvariation,
                        "shifts chip-dependency allowed?", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (fits_update_key(afptr, TINT, "SCHORDER", &chiporder,
                        "shifts chip fit order", &status)) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    hdunum = 2;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        status = 0;
        fprintf(stderr, " error moving to table of star information \n");
        fits_report_error(stderr, status);      /* print error message */
        exit(2);
    }


    /*
     * create new table 
     */
    fits_create_tbl(afptr, BINARY_TBL, ncoeffs, tfields, ttype, tform, tunit,
                    extname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }
    for(i = 0; i < ncoeffs; i++) {
        row = i + 1;
        fits_write_col(afptr, TDOUBLE, 1, row, 1, 1, &xc[i], &status);
        fits_write_col(afptr, TDOUBLE, 2, row, 1, 1, &yc[i], &status);
    }
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * close the fits file 
     */

    fits_close_file(afptr, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    return (0);

}




int main(int argc, char *argv[])
{
  int pwidth, pheight;
  int hdunum;
  double pscale;

  if (argc == 3)
    {
      hdunum = readimage(argv[1]);
    }
  else
    {
      printf(" %s <input image name> <output coefficients name>\n", argv[0]);
      exit(0);
    }

  pwidth = anaxes[0];
  pheight = anaxes[1];

  // set keyword values that lensfit expects to see, such that these values are for a 
  // polynomial of order 0 (=constant across the chip)

  int correct_distortion = 0;
  int ncoeffs = hdunum;
  int order = 0;
  int chipvariation = 1;
  int chiporder = 0;
  int nxchip, nychip;
  nxchip = hdunum==1 ? 1 : hdunum/4;  // assume four rows of chips
  nychip = hdunum==1 ? 1 : hdunum/nxchip;
  int numstar[1] = {10000};
  int xchipsampling = 2110;
  int ychipsampling = 4160;
  int big_gap = 0;
  double hxsize = xchipsampling*nxchip/2.;
  double hysize = ychipsampling*nychip/2.;
  float poserror = 6.;

  // convert scale factor to pixels per radian
  pscale = 180.*3600./M_PI/gs_scale[0];

  // write out to a lensfit PSF coefficients file using the standard lensfit function
  writeglobalcoeffs(argv[2], correct_distortion, spix, ncoeffs, order, chipvariation, chiporder, &pscale,
		    pwidth, pheight, nxchip, nychip, numstar, xchipsampling, ychipsampling,
		    big_gap, hxsize, hysize, poserror);


  // now write the table of astrometric shifts to the PSF coefficients file - in fact these are all zero here
  int sorder = 0;
  int sncoeffs = 1;
  int schipvariation = 0;
  int schiporder = 0;
  double xcoeffs[1], ycoeffs[1];

  xcoeffs[0] = ycoeffs[0] = 0.;

  updateglobalcoeffs(argv[2], sncoeffs, sorder, schipvariation,
                       schiporder, xcoeffs, ycoeffs);

  return 0;

}
