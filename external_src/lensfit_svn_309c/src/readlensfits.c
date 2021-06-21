/*-------------------------------------------------------------

readlensfits

example program to read the FITS table output from lensfits and 
do some simple analysis of the results.  

Lance Miller January 2008

modified for CARS  May 2008
modified for reading adaptive-length columns  LM  June 2008
modified for reading enhanced disk.bulge output LM July 2008
modified for reading PSF ellipticity  LM July 2008
this version fits sensitivity as a function of flux and size
   and for each galaxy calculates a flux-size weight 
   (independent of galaxy orientation)  LM Sep 2008
this twocomp version retains the two-component weight in the output
also allows for overlap between catalogues

new version calculates likelihood-based weights, weight components are equal but
still output as two-component LM March 2009

corrected weight bug, added mag, rh to output, LM 6 May 2009

changed to new prior, also reads catalogue magnitudes from file, now no need
for either input prior file or input catalogue.  LM 24 Sep 2010

Added option to filter on number of exposures for each galaxy, LM 29 Sep 2010.

Added output of version date string. LM 29 Sep 2010.

Modified prior.  Dec 2010. LM.

Catalogue format change 18 Dec 2010. LM

Write out psfevals and the exposure order.  CH 5th August 2011  

Change to weights and writing out ML estimator LM 25th Sept 2011

Changed to allow KiDS filenames and output of individual ellipticities to file, 2013.

Miscellaneous changes to output format, now in SVN control.

SVN Revision $Rev: 277 $ 
last changed by $Author: miller $
on $LastChangedDate: 2015-11-24 11:16:27 +0000 (Tue, 24 Nov 2015) $

--------------------------------------------------------------*/

#define DATE_STRING "$LastChangedDate: 2015-11-24 11:16:27 +0000 (Tue, 24 Nov 2015) $"
#define REV_STRING "$Rev: 277 $"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_erf.h>


/* define whether to write out FITS cube of posterior surfaces for each object
   (not recommended for large files > 1000 objects) */
#define WRITE_POSTERIOR_SURFACES 1

#define pi M_PI

void printerror(int);
int lensfits_tablesize(char *, int *, int *, double *);
int lensfits_read(char *, int, int, long *, int *, char ***, int**, float *, float *,
                  double *, double *, double **, double **, double **, int **, 
                  double *, int *, double *,
                  double *, double *, double *, double *, int *, double **,
                  int *, int **, float **, float*, int*);
int writeposterior_fits(char *, float *, int, int, int);
double efunc(double);
double bulgefunc(double);
double rfunc (double, float);
void svdfit2dsqc(double *, double *, double *, double *, int, int, int,
                 double *, double **, double **, double *);
void reconstruct(double, double, int, int, double *, double *);

// string for lensfit FITS-file version number
char lversion[10];
int versionflag;  // 0 for old format, 1 for new format tables
// OK to write out PSF ellipticities?
int write_psfs = 0; // 0=yes
int iversion; // lensfit version number
// default ellipticity prior parameter values (only used if not in fits table)
double efunc_a = 0.0256;    // intrinsic ellipticity term
double efunc_emax = 0.804;  // cut-off ellipticity 
double efunc_sigma = 0.2539; // cut-off sigma term
double bulgepopfrac = 0.1;   // fraction of pure-bulge galaxies
double bulgefunc_b = 2.3677769; // bulge prior parameter
double bulgefunc_c = 6.6907854; // bulge prior parameter
// default R-band prior parameter values (only used if not in fits table)
double afit[2]={-1.319698, -0.277907};  
double scalefactor = 0.186;
float refmag={23.0};

int
main(int argc, char *argv[])
{
    int tablesize[2];
    char tablename[500];
    int xsize, psize, ij;
    // int jj;
    int ifile, nfiles;
    int nume = 0;
    int *nobj, ntotal, i, j, ngood, nout, np, npmax, maxlen;
    int nwcut, nscut, necut, nbadin, npriorcut;
    int ix, iy, ns, hns, output_sampling;
    int **fitclass, **sampling, **numvals, ***evalues;
    int **duplicatefile, **duplicateobj;
    int **numexposures; //, minexposurenumber;
    int ***chip;
    int **id;
    int pdim, pval;
    //int iee1, iee2;
    int ***fflag;
    float **mag;
    double pstep, **eprior, prior, likel, post, edgef;
    double wngood, Pthreshold, Pthresholdlimit, nwecut;
    double **wcsx, **wcsy, ***psfe, ***psfe1, ***psfe2, ***probvals, **size, **bfrac, **flux;
    double **profileSN, **strehl, **snratio;
    double sumfv;
    double sume1, sume2, sump, sumpmoment, maxweight;
    double sump1, sump2, sumpp, peake1, peake2;
    double rwratio, iwratio;
    complex double complexpsf, wratio;
    double w1 ; //, emodlimit, sizelimit;
    double e_interval, norm, moment, wmoment, maxmoment, pmoment;
    double emod, emodp, emax, emax2, rmin, weightcut;
    double maxlikel, maxpost; // , maxloglikel;
    double x,y;
    // double rval, rprior, rnorm, eprob, pemod, Rarcsec;
    // int k, order, ncoeffs, ncoeffs2, nfit, nfitmax, nfitprevious;
    //double **u, **v, *w, *acoeff;
    //double *xfit, *yfit, *zfit, *wfit;
    // double xfitmin, xfitmax, yfitmin, yfitmax;
    // double edist, searchradius, fval, previous_maxl_e1, previous_maxl_e2;
    // double maxl_e1, maxl_e2, maxfitval_e1, maxfitval_e2;
    double tinterval;
    float *f, *sumf;
    float *e1 = NULL, *e2 = NULL;
    float ee1, ee2, ***posterior, *number;
    char *options, *fitsname, *sname, *pstr;
    char ****imagename;
    long nmaximages;
    FILE *ret, *ffile, *sfile, *evals=NULL;

  // output version information based on SVN info
    char version[100];
    bzero(version,100);
    strcpy(version, REV_STRING);
    int len; 
    len = strlen(version);
    strcpy(&version[len-2],"\0");
    printf("\n readlensfit version SVN:%s",&version[len-5]);
    char datedelims[]="()";
    char *dstring = NULL;
    char date_string[200];
    bzero(date_string,200);
    strcpy(date_string, DATE_STRING);
    dstring  = strtok(date_string, datedelims);
    if (dstring != NULL) dstring = strtok(NULL, datedelims);
    if (dstring != NULL) printf(" %s ",dstring);
    printf("\n");


    /*
     * set maximum allowed ellipticity 
     */
    emax = emax2 = 1.;
    /*
     * set lower limit to valid weight 
     */
    weightcut = 0.;
    /*
     * set minimum size 
     * this should be larger than the minimum size sampled in lensfit (Rmin)
     */
    rmin = 0.5;

    // joint size/ellipticity cut (not enabled)
    /*
    sizelimit = 2.3;
    emodlimit = 1.;  // emodlimit=1 => no cut applied
    */

    // Posterior threshold cut
    Pthresholdlimit = 0.01;
    printf (" Posterior threshold limit %f \n",Pthresholdlimit);

    if (argc < 3) {
        printf
            (" %s <fits file name(s)...> <output ascii filename>\n",
             argv[0]);
        exit(2);
    }

    /*
    // define maximum likelihood fit order and allocate memory
    order = 2;
    ncoeffs = (1 + order) * (2 + order) / 2;
    ncoeffs2 = 1000;
    acoeff = (double*)calloc(1+ncoeffs, sizeof(double));
    u = (double**)calloc(1+ncoeffs, sizeof(double*));
    v = (double**)calloc(1+ncoeffs, sizeof(double*));
    w = (double*)calloc(1+ncoeffs, sizeof(double));
    for (k=0; k<=ncoeffs; k++)
      {
	u[k] = (double*)calloc(1+ncoeffs,sizeof(double));
	v[k] = (double*)calloc(1+ncoeffs,sizeof(double));
      }
    nfitmax = 1000;
    xfit = (double*)calloc(nfitmax, sizeof(double));
    yfit = (double*)calloc(nfitmax, sizeof(double));
    zfit = (double*)calloc(nfitmax, sizeof(double));
    wfit = (double*)calloc(nfitmax, sizeof(double));
    */
    
    /*
   minexposurenumber = atoi(argv[1]);
    if (minexposurenumber <= 0)
      {
	fflush(stdout);
	fprintf(stderr," error specifying minimum number of exposures on command-line \n");
        fprintf(stderr,
            " %s <min number of exposures> <fits file name(s)...> <output ascii filename>\n",
             argv[0]);
	fprintf(stderr," minimum number must be >= 1 \n");
        exit(2);
      }
    printf(" specified minimum number of exposures per galaxy = %d \n",minexposurenumber);
    */

    nfiles = argc - 2;
    if (nfiles <= 0)
      {
	fflush(stdout);
	fprintf(stderr," error specifying input fits files on command-line \n");
        printf
            (" %s <min number of exposures> <fits file name(s)...> <output ascii filename>\n",
             argv[0]);
        exit(2);
      }


    nobj = (int *) calloc(nfiles, sizeof(int));
    ntotal = 0;

    /*
     * open the master file(s) 
     */

    fitsname = (char *) calloc(500, sizeof(char));

    if (access(argv[argc - 1], F_OK) == 0) {
        printf(" last named file on command-line exists already \n");
        exit(2);
    }

    if ((ret = fopen(argv[argc - 1], "w")) == NULL) {
        printf(" output file %s could not be opened\n", argv[argc - 1]);
        exit(3);
    }
    printf(" writing to file %s\n", argv[argc - 1]);
    
    // set the output sampling to choose for object list and fits cube
    output_sampling = 8;

    if (WRITE_POSTERIOR_SURFACES == 1)
      {
	sname = (char*)calloc(500, sizeof(char));
	strcpy(sname, argv[argc-1]);
	if (!(pstr = strrchr(sname, '.')))
	  pstr = sname + strlen(sname);
	sprintf(pstr, ".sampling%d.lis",output_sampling);
	if ((sfile = fopen(sname,"w"))==NULL){
	  fflush(stdout);
	  fprintf(stderr," sampled output file %s could not be opened\n",sname);
	  exit(2);
	}
	printf(" writing sampling=%d file to %s \n",output_sampling,sname);
      }
    
    if (write_psfs==0)
      {
	sname = (char*)calloc(500, sizeof(char));
	strcpy(sname, argv[argc-1]);
	if (!(pstr = strrchr(sname, '.')))
	  pstr = sname + strlen(sname);
	sprintf(pstr, ".evals.lis");
	if ((evals = fopen(sname,"w"))==NULL){
	  fflush(stdout);
	  fprintf(stderr," ellipticities output file %s could not be opened\n",sname);
	  exit(2);
	}
	printf(" writing ellipticities to file %s \n",sname);
      }

	fprintf(ret,"# lensfit catalogue %s \n",argv[1]);
	fprintf(ret,"#  1  wcsx\n");
	fprintf(ret,"#  2  wcsy\n");
	fprintf(ret,"#  3  <e1>\n");
	fprintf(ret,"#  4  <e2>\n");
	fprintf(ret,"#  5  weight \n");
	fprintf(ret,"#  6  fitclass \n");
	fprintf(ret,"#  7  scalelength /pixels \n");
	fprintf(ret,"#  8  bulge-fraction \n");
	fprintf(ret,"#  9  model-flux \n");
	fprintf(ret,"# 10  pixel SNratio \n");
	fprintf(ret,"# 11  model SNratio \n");
	fprintf(ret,"# 12  <PSF>-e1 \n");
	fprintf(ret,"# 13  <PSF>-e2 \n");
	fprintf(ret,"# 14  <PSF>-Strehl-ratio \n");
	fprintf(ret,"# 15  star-galaxy f-probability\n");
	fprintf(ret,"# 16  fit probability (not to be used for object selection)\n");
	fprintf(ret,"# 17  measurement variance \n");
	fprintf(ret,"# 18  mean-likelihood |e| \n");
	fprintf(ret,"# 19  mean likelihood <e1> = col 3 \n");
	fprintf(ret,"# 20  mean likelihood <e2> = col 4 \n");
	fprintf(ret,"# 21  catmag \n");
	fprintf(ret,"# 22  n-exposures-used \n");
	fprintf(ret,"# 23  cat-ID \n");

    ngood = 0;
    wngood = 0.;
    nwecut = 0.;
    npriorcut = 0;
    nbadin = nwcut = nscut = necut = 0;
    nout = 0;
    np = 0;

    eprior = (double **) calloc(nfiles, sizeof(double *));
    wcsx = (double **) calloc(nfiles, sizeof(double *));
    wcsy = (double **) calloc(nfiles, sizeof(double *));
    chip = (int ***) calloc(nfiles, sizeof(int **));
    fflag= (int ***) calloc(nfiles, sizeof(int **));
    psfe1 = (double ***) calloc(nfiles, sizeof(double **));
    psfe2 = (double ***) calloc(nfiles, sizeof(double **));
    psfe = (double ***) calloc(nfiles, sizeof(double **));
    strehl = (double **) calloc(nfiles, sizeof(double *));
    fitclass = (int **) calloc(nfiles, sizeof(int *));
    size = (double **) calloc(nfiles, sizeof(double *));
    bfrac = (double **) calloc(nfiles, sizeof(double *));
    flux = (double **) calloc(nfiles, sizeof(double *));
    snratio = (double **) calloc(nfiles, sizeof(double *));
    profileSN = (double **) calloc(nfiles, sizeof(double *));
    sampling = (int **) calloc(nfiles, sizeof(int *));
    numvals = (int **) calloc(nfiles, sizeof(int *));
    probvals = (double ***) calloc(nfiles, sizeof(double **));
    evalues = (int ***) calloc(nfiles, sizeof(int **));
    posterior = (float ***) calloc(nfiles, sizeof(float **));
    duplicatefile = (int **) calloc(nfiles, sizeof(int *));
    duplicateobj = (int **) calloc(nfiles, sizeof(int *));
    imagename = (char ****) calloc(nfiles, sizeof(char ***));
    mag = (float **) calloc(nfiles, sizeof(float *));
    id = (int **) calloc(nfiles, sizeof(int *));
    numexposures = (int**)calloc(nfiles, sizeof(int*));

    for(ifile = 0; ifile < nfiles; ifile++) {
        printf(" opening %s\n", argv[ifile + 1]);
        options = (char *) "r";
        /*
         * read fits file names in
         */
        ffile = fopen(argv[ifile + 1], options);
        if (ffile == NULL)
            perror(" Error input opening file");

        strcpy(tablename, argv[ifile + 1]);

        if (lensfits_tablesize(tablename, tablesize, &maxlen, &e_interval)) {
            fprintf(stderr,
                    " error reading size of lensfit output table %s\n",
                    tablename);
            exit(2);
        }

        if (WRITE_POSTERIOR_SURFACES == 1) {
            strcpy(fitsname, tablename);
            strcat(fitsname, "_posterior.fits");
        }


        if (ifile > 0) {
            if (nume != tablesize[0]) {
                fprintf(stderr,
                        " input fits tables have different ellipticity grids \n");
                exit(2);
            }
        } else {
            nume = tablesize[0];

            e1 = (float *) calloc(nume, sizeof(float));
            e2 = (float *) calloc(nume, sizeof(float));
            number = (float *) calloc(nume, sizeof(float));
            if (e1 == NULL || e2 == NULL) {
                fprintf(stderr, " memory allocation error for e1,e2 \n");
                exit(2);
            }
        }

        nobj[ifile] = tablesize[1];
	if (WRITE_POSTERIOR_SURFACES == 1)
	  {  // output postage stamp size
	    xsize = 1 + (int) (0.1 + 2. / (output_sampling*e_interval));
	    psize = xsize * xsize;
	    npmax = 2000;
	    f = (float *) calloc(psize * (npmax+1), sizeof(float));
	    sumf = (float *) calloc(psize, sizeof(float));
	  }

        wcsx[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        wcsy[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
	chip[ifile] = (int **) calloc(nobj[ifile], sizeof(int *));
        psfe[ifile] = (double **) calloc(nobj[ifile], sizeof(double *));
	fflag[ifile] = (int **) calloc(nobj[ifile], sizeof(int *));
        psfe1[ifile] = (double **) calloc(nobj[ifile], sizeof(double *));
        psfe2[ifile] = (double **) calloc(nobj[ifile], sizeof(double *));

	strehl[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        fitclass[ifile] = (int *) calloc(nobj[ifile], sizeof(int));
        size[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        bfrac[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        flux[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        snratio[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        profileSN[ifile] = (double *) calloc(nobj[ifile], sizeof(double));
        sampling[ifile] = (int *) calloc(nobj[ifile], sizeof(int));
        numvals[ifile] = (int *) calloc(nobj[ifile], sizeof(int));
        probvals[ifile] = (double **) calloc(nobj[ifile], sizeof(double *));
        evalues[ifile] = (int **) calloc(nobj[ifile], sizeof(int *));
        posterior[ifile] = (float **) calloc(nobj[ifile], sizeof(float *));
        duplicatefile[ifile] = (int *) calloc(nobj[ifile], sizeof(int));
        duplicateobj[ifile] = (int *) calloc(nobj[ifile], sizeof(int));
        imagename[ifile] = (char ***) calloc(nobj[ifile], sizeof(char **));
	mag[ifile] = (float *) calloc(nobj[ifile], sizeof(float));
	id[ifile] = (int *) calloc(nobj[ifile], sizeof(int));
	numexposures[ifile] = (int*)calloc(nobj[ifile], sizeof(int));

        if (wcsx[ifile] == NULL || wcsy[ifile] == NULL || id[ifile] == NULL
            || fitclass[ifile] == NULL || sampling[ifile] == NULL
            || strehl[ifile] == NULL || numvals[ifile] == NULL
            || probvals[ifile] == NULL || evalues[ifile] == NULL
            || posterior[ifile] == NULL || size[ifile] == NULL
            || bfrac[ifile] == NULL || flux[ifile] == NULL
	    || mag[ifile] == NULL || fflag[ifile] == NULL
            || profileSN[ifile] == NULL || psfe1[ifile] == NULL) {
            fprintf(stderr, " memory allocation error for second table \n");
            exit(2);
        }
        for(i = 0; i < nobj[ifile]; i++) 
	  {
	    psfe[ifile][i] = (double *) calloc(2, sizeof(double));
            probvals[ifile][i] = (double *) calloc(2, sizeof(double));
            evalues[ifile][i] = (int *) calloc(maxlen, sizeof(int));
            posterior[ifile][i] = (float *) calloc(maxlen, sizeof(float));
            if (psfe[ifile][i] == NULL || probvals[ifile][i] == NULL
                || posterior[ifile][i] == NULL) {
                fprintf(stderr, " memory allocation error, second table \n");
                exit(2);
            }
        }

        if (lensfits_read
            (tablename, nume, nobj[ifile], &nmaximages, numexposures[ifile], 
	     imagename[ifile], chip[ifile], e1,
             e2, wcsx[ifile], wcsy[ifile], psfe[ifile], psfe1[ifile], psfe2[ifile], fflag[ifile], strehl[ifile],
             fitclass[ifile], size[ifile], bfrac[ifile], flux[ifile],
             snratio[ifile], profileSN[ifile], sampling[ifile], probvals[ifile],
             numvals[ifile], evalues[ifile], posterior[ifile],
	     mag[ifile],id[ifile])) {
            fprintf(stderr, " error reading lensfit output table %s\n",
                    tablename);
            exit(2);
        }

        ntotal += nobj[ifile];

        fclose(ffile);

	// create prior grid at sampling 0.001
	// (assume all files have the same prior)
	pstep = 0.001;
	pdim = 1 + (int)(1./pstep);
	eprior[ifile] = (double*)calloc(pdim, sizeof(double));
	norm = pmoment = 0.;
	for (j=0; j<pdim; j++)
	  {
	    emod = j*pstep;
	    eprior[ifile][j] = 
	      (1.-bulgepopfrac)*exp(efunc(emod)) + 
		  bulgepopfrac*exp(bulgefunc(emod));
	    // printf(" %d %f %f \n",j,emod,eprior[ifile][j]); fflush(stdout);
	    // calculate moments
	    norm += eprior[ifile][j]*emod;
	    pmoment += eprior[ifile][j]*pow(emod,3);
	  }
	pmoment /= norm;
	printf(" prior sigma = %f \n",sqrt(pmoment/2.));

	/*
	if ((int)nmaximages < minexposurenumber)
	  {
	    fflush(stdout);
	    fprintf(stderr," specified minimum number of exposures %d > actual number %d \n",
		    minexposurenumber,(int)nmaximages);
	    exit(2);
	  }
	*/
    }

/*-----------------------------------------------------------------*/


    /*
     * now go through lensfit and catalogue files, check they match, and
     * calculate the final weights etc 
     */

    for(ifile = 0; ifile < nfiles; ifile++) {
      for(i = 0; i < nobj[ifile]; i++) {
	sume1 = sume2 = sump = 0.;
	sump1 = sump2 = sumpp = sumpmoment = 0.;
	peake1 = peake2 = 0.;
	moment = 0.;
	// maxl_e1 = maxl_e2 = maxfitval_e1 = maxfitval_e2 = 0.;

	tinterval = e_interval*sampling[ifile][i];

	ns = sampling[ifile][i];
	hns = ns/2;

	// calculate mean from likelihood
	if (numvals[ifile][i] > 0) {
	  maxlikel = maxpost = -1.e10;
	  for(j = 0; j < numvals[ifile][i]; j++) 
	    { // work through tabulated values and calculate mean ellipticities etc
	      ee1 = e1[evalues[ifile][i][j]];
	      ee2 = e2[evalues[ifile][i][j]];

	      // evaluate prior to calculate likelihood from posterior
	      emod = sqrt(ee1 * ee1 + ee2 * ee2);
	      pval = (int)(0.5 + emod/pstep);
	      prior = eprior[ifile][pval];

	      // test for table version and dal with each appropriately
	      if(versionflag==1)
		{ // later version, assume likelihood surface
		  likel = exp(posterior[ifile][i][j]);
		  post = likel*prior;
		}
	      else
		{ // assume posterior surface
		  post = exp(posterior[ifile][i][j]);
		  if (prior > 0.)
		    {
		      likel = post/prior;		  
		    }
		  else
		    {
		      likel = 0.;
		    }
		}

	      // find maximum posterior
	      if (post > maxpost) 
		{
		  maxpost = post;
		  //maxl_e1 = ee1;
		  //maxl_e2 = ee2;
		}

	    }

	  // make threshold on posterior to be applied when calculating likelihood estimator
	  Pthreshold = Pthresholdlimit*maxpost;

	  for(j = 0; j < numvals[ifile][i]; j++) 
	    { // work through tabulated values and calculate mean ellipticities etc
	      ee1 = e1[evalues[ifile][i][j]];
	      ee2 = e2[evalues[ifile][i][j]];

	      // evaluate prior to calculate likelihood from posterior
	      emod = sqrt(ee1 * ee1 + ee2 * ee2);
	      pval = (int)(0.5 + emod/pstep);
	      prior = eprior[ifile][pval];

	      // test for table version and dal with each appropriately
	      if(versionflag==1)
		{ // later version, assume likelihood surface
		  likel = exp(posterior[ifile][i][j]);
		  post = likel*prior;
		}
	      else
		{ // assume posterior surface
		  post = exp(posterior[ifile][i][j]);
		  if (prior > 0.)
		    {
		      likel = post/prior;		  
		    }
		  else
		    {
		      likel = 0.;
		    }
		}

	      // calculate mean and moment
	      if (likel > 0.) 
		{
		  // finely sample the array down to e_interval and calculate means
		  for (iy=-hns; iy<=hns; iy++)
		    {
		      y = (int)(1000+(ee2/e_interval)+0.5) + iy - 1000;
		      y *= e_interval;
		      for (ix=-hns; ix<=hns; ix++)
			{
			  x = (int)(1000+(ee1/e_interval)+0.5) + ix - 1000;
			  x *= e_interval;
			  // recalculate prior for this fine grid
			  emod = sqrt(x*x + y*y);
			  pval = (int)(0.5 + emod/pstep);
			  prior = eprior[ifile][pval];
			  // edge/corner-weighting factor
			  edgef=1.;
			  if ((iy==-hns || iy==hns) && hns>0) edgef *= 0.5;
			  if ((ix==-hns || ix==hns) && hns>0) edgef *= 0.5;
			  // apply posterior cut
			  if (post >= Pthreshold)
			    {
			      // accumulate mean likelihood
			      sume1 += edgef * ee1 * likel;
			      sume2 += edgef * ee2 * likel;
			      sump += edgef * likel;
			      // accumulate mean posterior
			      sump1 += edgef * ee1 * likel * prior;
			      sump2 += edgef * ee2 * likel * prior;
			      sumpp += edgef * likel * prior;
			    }
			  // accumulate moment including excluded pixels
			  moment += edgef*emod*emod*likel;
			  sumpmoment += edgef*likel;
			}
		    }
		}
	    }

	  /*
	  // try walking uphill in likelihood from the mean likelihood position
	  // if available (otherwise max posterior position)
	  // to find the best maximum likelihood position
	  if (sump>0.)
	    {
	      maxl_e1 = sume1/sump;
	      maxl_e2 = sume2/sump;
	    }
	  tinterval = e_interval*sampling[ifile][i];
	  previous_maxl_e1 = previous_maxl_e2 = -1.;
	  while (maxl_e1 != previous_maxl_e1 || maxl_e2 != previous_maxl_e2)
	    {
	      previous_maxl_e1 = maxl_e1;
	      previous_maxl_e2 = maxl_e2;
	      for(j = 0; j < numvals[ifile][i]; j++) 
		{ // work through tabulated values and calculate mean ellipticities etc
		  ee1 = e1[evalues[ifile][i][j]];
		  ee2 = e2[evalues[ifile][i][j]];
		  // evaluate prior to calculate likelihood from posterior
		  emod = sqrt(ee1 * ee1 + ee2 * ee2);
		  pval = (int)(0.5 + emod/pstep);
		  prior = eprior[ifile][pval];
		  // test for table version and dal with each appropriately
		  if(versionflag==1)
		    { // later version, assume likelihood surface
		      likel = exp(posterior[ifile][i][j]);
		      post = likel*prior;
		    }
		  else
		    { // assume posterior surface
		      post = exp(posterior[ifile][i][j]);
		      if (prior > 0.)
			{
			  likel = post/prior;		  
			}
		      else
			{
			  likel = 0.;
			}
		    }
		  // apply posterior threshold
		  if (likel > 0. && post >= Pthreshold && emod<efunc_emax) 
		    {
		      // test if this point has a higher likelihood
		      if ( fabs(ee1-maxl_e1)<=tinterval+e_interval/2. 
			   && fabs(ee2-maxl_e2)<=tinterval+e_interval/2. )
			{
			  // identify the maximum likelihood
			  if (likel>maxlikel)
			    {
			      maxlikel = likel;
			      maxl_e1 = ee1;
			      maxl_e2 = ee2;
			    }
			}
		    }
		}
	    }


	  // go through values again and find ellipticity pixels that connect
	  // contiguously to the maximum pixel
	  // start with some small search radius and increase this until sufficient
	  // samples have been obtained for a fit of given order

	  nfitprevious = -1;
	  nfit = 0;
	  searchradius = tinterval;
	  xfitmax = yfitmax = -1.;
	  xfitmin = yfitmin = 1.;
	  // make sure there are enough points which cover a minimum range in both x and y
	  // but check that we have not reached the end of the measured set
	  // and make sure that new points are beign added at each step
	  while ( (nfit<ncoeffs || (xfitmax-xfitmin)<(1+order)*tinterval ||
		   (yfitmax-yfitmin)<(1+order)*tinterval ) &&
		  nfit<numvals[ifile][i] && nfit<nfitmax && nfit>nfitprevious)
	    {
	      searchradius += tinterval;
	      //printf(" nfit %d %d search %f \n",nfit,numvals[ifile][i],searchradius);
	      nfitprevious = nfit;
	      nfit = 0;
	      for(j = 0; j < numvals[ifile][i]; j++) 
		{ // work through tabulated values and calculate mean ellipticities etc
		  ee1 = e1[evalues[ifile][i][j]];
		  ee2 = e2[evalues[ifile][i][j]];
		  edist = sqrt(pow((ee1-maxl_e1),2) + pow((ee2-maxl_e2),2));
		  if (edist <= searchradius)
		    {
		      // evaluate prior to calculate likelihood from posterior
		      emod = sqrt(ee1 * ee1 + ee2 * ee2);
		      pval = (int)(0.5 + emod/pstep);
		      prior = eprior[ifile][pval];
		      // test for table version and deal with each appropriately
		      if(versionflag==1)
			{ // later version, assume likelihood surface
			  likel = exp(posterior[ifile][i][j]);
			  post = likel*prior;
			}
		      else
			{ // assume posterior surface
			  post = exp(posterior[ifile][i][j]);
			  if (prior > 0.)
			    {
			      likel = post/prior;		  
			    }
			  else
			    {
			      likel = 0.;
			    }
			}
		      // apply posterior threshold
		      if (likel > 0. && post >= Pthreshold) 
			{
			  if (nfit<nfitmax)
			    {
			      xfit[nfit] = ee1;
			      yfit[nfit] = ee2;
			      zfit[nfit] = log(likel);
			      wfit[nfit] = 1.;
			      if (ee1<xfitmin) xfitmin = ee1;
			      if (ee1>xfitmax) xfitmax = ee1;
			      if (ee2<yfitmin) yfitmin = ee2;
			      if (ee2>yfitmax) yfitmax = ee2;
			    }
			  nfit++;
			}
		    }
		}
	    }
	  // set default maximum to maximum discrete location
	  maxfitval_e1 = maxl_e1;
	  maxfitval_e2 = maxl_e2;
	  maxloglikel = log(maxlikel);
	  // if enough samples and enough x,y range, fit a bicubic surface
	  if ( nfit>=ncoeffs && (xfitmax-xfitmin)>=(1+order)*tinterval &&
	       (yfitmax-yfitmin)>=(1+order)*tinterval )
	    {
	      // SVD fit
	      svdfit2dsqc(xfit,yfit,zfit,wfit,nfit,order,0,acoeff,u,v,w);
	      // find maximum near discrete maximum
	      for (iee1=-50; iee1<=50; iee1++)
		{
		  ee1 = maxl_e1 + iee1*tinterval/50.;
		  for (iee2=-50; iee2<=50; iee2++)
		    {
		      ee2 = maxl_e2 + iee2*tinterval/50.;
		      emod = sqrt(ee1*ee1+ee2*ee2);
		      if (emod < 1.)
			{
			  // reconstruct remembering the shift in array convention for acoeff
			  reconstruct(ee1,ee2,order,0,&acoeff[1],&fval);
			  if (fval>maxloglikel)
			    {
			      maxfitval_e1 = ee1;
			      maxfitval_e2 = ee2;
			      maxloglikel = fval;
			    }
			}
		    }
		}

	    }
	  // otherwise try to fit a quadratic if enough points and ncoeffs2<=ncoeffs
	  else if (nfit>=ncoeffs2)
	    {
	      // SVD fit
	      svdfit2dsqc(xfit,yfit,zfit,wfit,nfit,2,0,acoeff,u,v,w);
	      // find maximum near discrete maximum
	      for (iee1=-50; iee1<=50; iee1++)
		{
		  ee1 = maxl_e1 + iee1*tinterval/50.;
		  for (iee2=-50; iee2<=50; iee2++)
		    {
		      ee2 = maxl_e2 + iee2*tinterval/50.;
		      // reconstruct remembering the shift in array convention for acoeff
		      reconstruct(ee1,ee2,2,0,&acoeff[1],&fval);
		      if (fval>maxloglikel)
			{
			  maxfitval_e1 = ee1;
			  maxfitval_e2 = ee2;
			  maxloglikel = fval;
			}
		    }
		}

	    }
	  */

	}
	else 
	  { // case of no measured likelihood values
	    //nfit = 0;
	    //peake1 = 0.;
	    //peake2 = 0.;
	    sume1 = 0.;
	    sume2 = 0.;
	    sump = 0.;
	    moment = 0.;
	    //maxl_e1=maxl_e2=maxfitval_e1=maxfitval_e2=0.;
	  }

	maxmoment = efunc_emax*efunc_emax/2.;
	wmoment = maxmoment;
	if (sump > 0. && sumpmoment > 0.) {
	  sume1 /= sump;
	  sume2 /= sump;
	  moment /= sumpmoment;
	  moment -= (sume1*sume1 + sume2*sume2);
	  if (moment < tinterval*tinterval) moment = tinterval*tinterval;
	  wmoment = moment;
	  if (moment < maxmoment)
	    {
	      moment = moment*maxmoment/(maxmoment-moment);
	      moment += pmoment;
	    }
	  else
	    {
	      //printf ("object %d moment %f > maxmoment %f\n",i+1,moment,maxmoment);
	      moment = 0.;
	    }
	} else {
	  sume1 = 0.;
	  sume2 = 0.;
	  moment = 0.;
	}

	if (sumpp > 0.) {
	  sump1 = sump1 / sumpp;
	  sump2 = sump2 / sumpp;
	} else {
	  sump1 = 0.;
	  sump2 = 0.;
	}

	emod = sqrt(sume1*sume1 + sume2*sume2);
	emodp = sqrt(sump1*sump1 + sump2*sump2);


	moment /= 2.; // convert to 1D moment

	// set the inverse-variance weight
	if (moment!=0. && fitclass[ifile][i]==0)
	  w1 = 1./moment;
	else
	  w1 = 0.;

	// cut high weight objects if they have low S/N
	// (probably caused by blending with background objects or artefacts)
	if (iversion > 228)
	  {
	    if (wmoment < 1./pow(profileSN[ifile][i],2) ) 
	      {
		w1 = 0.;
		if (fitclass[ifile][i]==0)
		  {
		    nwcut++;
		    fitclass[ifile][i] = -7;
		  }
	      }
	  }

	// evaluate prior probability of a galaxy being more elliptical than the observed value
	/*
	j = pdim-1;
	pemod = j*pstep;
	eprob = 0.;
	if (pemod >= emod)
	  {
	    eprob =  eprior[ifile][j]*pemod;
	    while (j>0 && pemod > emod) 
	      {
		j--;
		pemod = j*pstep;
		if (pemod >= emod)
		  eprob +=  eprior[ifile][j]*pemod;
	      }
	    eprob /= norm;
	  }
	*/

	// similarly for prior probability of scalelength
	/*
	rprior = 0.;
	rnorm = 0.;
	for (rval=5.; rval>=rmin; rval-=0.01)
	  {
	    Rarcsec = rval*scalefactor;
	    rnorm += exp(rfunc(Rarcsec, mag[ifile][i]));
	    if (rval >= size[ifile][i])
	      rprior += exp(rfunc(Rarcsec, mag[ifile][i]));
	  }
	rprior /= rnorm;
	*/

	// remove unlikely objects
	/*
	if (rprior*eprob < 1.e-2 || eprob < 1.e-2 || rprior < 1.e-2) 
	  {
	    if (w1>0.)
	      {
		printf("mag-size %f %f %f\n",mag[ifile][i],size[ifile][i],w1);
		fitclass[ifile][i] = -7;
		npriorcut++;
		w1 = 0.;
	      }
	  }
	*/

	// apply the size cut
	if (size[ifile][i] < rmin)
	  {
	    if (w1 > 0.) nscut++;
	    w1 = 0.;
	  }

	/*
	// apply the weight cut
	if (w1 < weightcut) 
	  {
	    if (w1 > 0.) nwcut++;
	    w1 = 0.;
	  }
	*/

	// add up numbers and weighted numbers of good galaxies
	if (w1>0.) 
	  {
	    ngood++;
	    wngood += w1;
	  }

	if (WRITE_POSTERIOR_SURFACES==1 && sampling[ifile][i] == output_sampling && w1>0.)
	  {
	    fprintf(sfile,"%10.6lf %10.6lf %5.2f %8d\n",
		    wcsx[ifile][i],wcsy[ifile][i],mag[ifile][i],id[ifile][i]);
	  }

	rwratio = iwratio = 0.;
	wratio = 0.+I*0.;
	complexpsf = psfe[ifile][i][0] + I*psfe[ifile][i][1];
	if (cabs(complexpsf)>0.) 
	  {
	    wratio = (sume1 - I*sume2)*complexpsf;
	    rwratio = creal(wratio);
	  }

	/*
	if (w1>0. && size[ifile][i]>sizelimit && emod>emodlimit)
	  {
	    nwecut += w1;
	    w1=0.;
	    necut++;
	  }
	*/

	/*
	if (w1>10. && mag[ifile][i] > 23.5)
	  {
	    printf("%10.6lf %10.6lf %5.2f %8d\n",
		    wcsx[ifile][i],wcsy[ifile][i],mag[ifile][i],id[ifile][i]);
	  }
	*/

	//emod = sqrt(maxfitval_e1*maxfitval_e1+maxfitval_e2*maxfitval_e2);
	emod = sqrt(sume1*sume1 + sume2*sume2);

	fprintf(ret,
		"%10.6lf %10.6lf %7.4f %7.4f %7.4f %2d %7.4f %6.4f %8.1f %8.1f %8.6f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %5.2f %3d %8d", // %7.4f %7.4f %7.4f %7.4f %8d",
		// wcsx[ifile][i], wcsy[ifile][i], maxfitval_e1, maxfitval_e2, moment,
		wcsx[ifile][i], wcsy[ifile][i], sume1, sume2, w1,
		fitclass[ifile][i], size[ifile][i], bfrac[ifile][i],
		flux[ifile][i], snratio[ifile][i], profileSN[ifile][i],
		psfe[ifile][i][0], psfe[ifile][i][1],
		strehl[ifile][i], probvals[ifile][i][0], probvals[ifile][i][1],
		wmoment, emod,
		sume1, sume2,
                mag[ifile][i],	
		numexposures[ifile][i],
		id[ifile][i]
		// ,maxl_e1, maxl_e2, 
		// tinterval, nfit
		);
	
	/* write out the individual PSF ellipticities */
	if (write_psfs==0 && evals != NULL)
	  {
	    fprintf(evals," %d ",i+1);
	    for (j=0; j<nmaximages; j++)
	      {
		if (fflag[ifile][i][j]==1)
		  {
		    fprintf(evals," %8.5f %8.5f ",psfe1[ifile][i][j],psfe2[ifile][i][j]);
		  }
		else
		  {
		    fprintf(evals," 99.00000 99.00000 ");
		  }
	      }
	    fprintf(evals,"\n");
	  }

	/*
	for(jj = 0; jj < nmaximages; jj++)
	  fprintf(ret, " %s", imagename[ifile][i][jj]);
	*/
	fprintf(ret, "\n");
	  

	if (WRITE_POSTERIOR_SURFACES == 1) {
	  if (sampling[ifile][i] >= output_sampling && 
	      w1 > 0. &&
	      np < npmax )
	    {
	    sumfv = 0.;
	    for(j = 0; j < numvals[ifile][i]; j++) {
	      ee1 = e1[evalues[ifile][i][j]];
	      ee2 = e2[evalues[ifile][i][j]];
	      emod = sqrt(ee1 * ee1 + ee2 * ee2);
	      pval = (int)(0.5 + emod/pstep);
	      prior = eprior[ifile][pval];
	      if(versionflag==0)
		{ // posterior surface, calculate likelihood
		  if (prior > 0.) {
		    if (posterior[ifile][i][j] > -1.e10) {
		      sumfv +=
			exp(posterior[ifile][i][j]) / prior;
		    }
		  }
		}
	      else
		{ // likelihood surface
		  sumfv += exp(posterior[ifile][i][j]);
		}
	    }
	    for(j = 0; j < numvals[ifile][i]; j++) {
	      ee1 = e1[evalues[ifile][i][j]];
	      ee2 = e2[evalues[ifile][i][j]];
	      emod = sqrt(ee1 * ee1 + ee2 * ee2);
	      pval = (int)(0.5 + emod/pstep);
	      prior = eprior[ifile][pval];
	      iy = (int)(xsize + 0.5 + ee2/(sampling[ifile][i]*e_interval))-xsize;
	      ix = (int)(xsize + 0.5 + ee1/(sampling[ifile][i]*e_interval))-xsize;
	      ij = np*psize + (iy+xsize/2)*xsize + ix + xsize/2;
	      printf("%d %d %f %f %d %d %g \n",i,j,ee1,ee2,ix,iy,posterior[ifile][i][j]);
		if (posterior[ifile][i][j] > -1.e10) {
		  if(versionflag==0)
		    { // posterior surface, calculate likelihood
		      if (prior > 0.)
			{
			  f[ij] = exp(posterior[ifile][i][j]) /prior /sumfv;
			  sumf[(iy+xsize/2)*xsize + ix + xsize/2] += f[ij] /sumfv;
			}
		      else
			{
			  f[ij] = 0.;
			}
		    }
		  else
		    { // likelihood surface
		      f[ij] = exp(posterior[ifile][i][j]);
		      sumf[(iy+xsize/2)*xsize + ix + xsize/2] += f[ij] /sumfv;
		    }
		}
		else 
		  {
		    f[ij] = 0.;
		  }
	    }
	    np++;
	    printf(" %d %d %lf %lf %f \n",np,id[ifile][i],wcsx[ifile][i],wcsy[ifile][i],w1);
	    }
	}
	// printf (" %d %d \n",nout,i);
	nout++;
	    
	}
      }
    
    if (WRITE_POSTERIOR_SURFACES == 1) {
        for(i = 0; i < xsize * xsize; i++) {
            f[np * psize + i] = sumf[i];
        }
        np++;
        writeposterior_fits(fitsname, f, np, xsize, xsize);
	fclose(sfile);
    }

    fclose(ret);
    // printf(" %d objects set to zero weight on input \n", nbadin);

    // find maximum allowed weight for galaxies
    moment = e_interval*e_interval;
    moment = pmoment + moment*maxmoment/(maxmoment-moment);
    maxweight = 2./moment;
    printf(" maximum possible galaxy weight = %f \n",maxweight);

    /*
    if (emodlimit<1.)
      {
	printf(" %d galaxies with size > %f cut by maximum ellipticity %f \n", necut,sizelimit,emodlimit);
	printf(" %f weighted galaxies cut by above\n",nwecut/maxweight);
      }
    */
    if (rmin>0.)
      printf(" %d galaxies cut by minimum size %f \n", nscut,rmin);
    if (nwcut>0)
      printf(" %d high-weight galaxies cut by S/N limit \n", nwcut);
    if (npriorcut > 0)
      printf(" %d galaxies cut by joint size/ellipticity prior cut\n",npriorcut);
    printf(" %d objects written out \n", nout);
    printf(" %d good galaxies \n", ngood);
    printf(" %f weighted galaxies \n", wngood/maxweight);

    if (ngood <= 0)
        exit(2);

    return 0;
    
}


double bulgefunc (double eval)
{
  /* function for log(ellipticity prior) for bulge-dominated galaxies, assumed
     independent of magnitude, based on fits of Simard et al. 2002
  */
  double prior, pnorm, erfarg1, erfarg2;
  gsl_sf_result xc, yc;

  if (eval<1.)
    {

  // basic log(prior) function
  prior = -bulgefunc_c*eval*eval - bulgefunc_b*eval;

  // normalisation from mathematica
  erfarg1 = bulgefunc_b/2./sqrt(bulgefunc_c);
  erfarg2 = (bulgefunc_b+2.*bulgefunc_c)/2./sqrt(bulgefunc_c);
  if (gsl_sf_erf_e(erfarg1, &yc)==0 && gsl_sf_erf_e(erfarg2, &xc)==0)
    {
      pnorm = (yc.val-xc.val);
      if (pnorm < 0.)
	{
	  pnorm *= sqrt(pi)*bulgefunc_b*
	    exp(bulgefunc_b*bulgefunc_b/4./bulgefunc_c)/
	    4./bulgefunc_c/sqrt(bulgefunc_c);
	  pnorm += 1./(1.+1./(exp(bulgefunc_b+bulgefunc_c)-1.))/2./bulgefunc_c;
	  pnorm *= 2.*pi;
	  if (pnorm > 0.)
	    {
	      pnorm = log(pnorm);
	    }
	  else
	    {
	      fflush(stdout);
	      fprintf(stderr," negative norm in bulgefunc \n");
	      exit(2);
	    }
	}
      else
	{
	  fflush(stdout);
	  fprintf(stderr," rounding error in erf in bulgefunc\n");
	  exit(2);
	}
    }
  else
    {
      fflush(stdout);
      fprintf(stderr," error calculating erf in bulgefunc \n");
      exit(2);
    }

  // return log(normalised prior)
  prior -= pnorm;

    }
  else
    {
      prior = -1.e10;
    }

  return (prior);

}

double efunc (double eval)
{
  /* 
   *  function for log(ellipticity prior), for disk galaxies, assumed independent 
   *  of magnitude, based on fit to SDSS disk axis ratios following Unterborn
   *  and Ryden 2008.  This fit is actually based on the distribution of
   *  axis ratios in the g-band DR7 data (see programs qdist.for and mockprior.for)
   */
  /* NB this is only approximately normalised */
  double prior, pnorm, arg1, arg2;
  gsl_sf_result xc, yc;

  if (eval < efunc_emax)
    {
      // high ellipticity cut-off (note revised functional form with (1+e) term)
      prior = (1.-exp((eval-efunc_emax)/efunc_sigma))/(1.+eval);
      // intrinsic ellipticity term
      prior *= 1./sqrt((eval*eval + efunc_a*efunc_a));
    }
  else
    {
      prior = 1.e-10;
    }

  // approx normalisation assuming efunc_a is small
  //pnorm = 1.-efunc_sigma*(exp((1.-efunc_emax)/efunc_sigma)-exp(-efunc_emax/efunc_sigma));
  // normalisation of revised function with extra (1+e) term
  arg1 = -(1.+efunc_emax)/efunc_sigma;
  arg2 = -1./efunc_sigma;
  if (gsl_sf_expint_E1_e(arg1, &yc)==0 && gsl_sf_expint_E1_e(arg2, &xc)==0)
    {
      pnorm = log(1.+efunc_emax) + exp(arg1)*(yc.val-xc.val);
    }
  else
    {
      fflush(stdout);
      fprintf(stderr," error normalising efunc prior \n");
      exit(2);
    }
  pnorm *= (sqrt(efunc_a*efunc_a+1.)-efunc_a);
  pnorm *= 2.*pi;
	       
  return log(prior/pnorm) ;
}



int
lensfits_tablesize(char *tablename, int *tablesize, int *colmaxlen,
                   double *keyvals)
{
    fitsfile *fptr;             /* pointer to the FITS file, defined in fitsio.h */
    int status, hdunum, hdutype, ncol, maxlen;
    long nume, nobj;
    char comment[64], tform19[64], *vlen = NULL;
    double r_afit0, r_afit1, pscale;
    float rfmag;

    status = 0;

    if (fits_open_file(&fptr, tablename, READONLY, &status))
        printerror(status);

    /*
     * move to the HDU 
     */
    hdunum = 2;
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status))
        printerror(status);

    if (hdutype != BINARY_TBL) {
        printf(" not a binary table \n");
        exit(2);
    }

    fits_get_num_cols(fptr, &ncol, &status);
    if (ncol != 2) {
        fprintf(stderr, " error in table format, no. columns = %d \n", ncol);
        exit(2);
    }

    fits_get_num_rows(fptr, &nume, &status);

    if (fits_read_key
        (fptr, TDOUBLE, "DELTA-E", &keyvals[0], comment, &status)) {
        fits_report_error(stderr, status);
        printf(" keyword DELTA-E\n");
        return (status);
    }

    if (fits_read_key
        (fptr, TSTRING, "VERSION", lversion, comment, &status)) 
      {
	printf(" lensfit version earlier than 7.1 \n");
	bzero(lversion,10);
	status = 0;
	versionflag=0;
	iversion = 0;
      }
    else
      {
	versionflag=1;
	vlen = strtok(lversion,":");
	vlen = strtok(NULL,":");
	iversion = atoi(vlen);
	printf(" reading FITS file from lensfit version %d \n",iversion);
	if (fits_read_key
	    (fptr, TDOUBLE, "EMAX", &efunc_emax, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword EMAX\n");
	  return (status);
	}
	if (fits_read_key
	    (fptr, TDOUBLE, "EFUNCA", &efunc_a, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword EFUNCA\n");
	  return (status);
	}
	if (fits_read_key
	    (fptr, TDOUBLE, "EFUNCS", &efunc_sigma, comment, &status)) {
        fits_report_error(stderr, status);
        printf(" keyword EFUNCS\n");
        return (status);
	}
	if (fits_read_key
	    (fptr, TDOUBLE, "BPOPFRAC", &bulgepopfrac, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword BPOPFRAC\n");
	  //return (status);
	  status=0;
	}
	if (fits_read_key
	    (fptr, TDOUBLE, "BFUNCB", &bulgefunc_b, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword BFUNCB\n");
	  return (status);
	}
	if (fits_read_key
	    (fptr, TDOUBLE, "BFUNCC", &bulgefunc_c, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword BFUNCC\n");
	  return (status);
	}
	if (fits_read_key
	    (fptr, TDOUBLE, "RFIT0", &r_afit0, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword RFIT0\n");
	  status = 0;
	}
	else
	  {
	    afit[0] = r_afit0;
	  }
	if (fits_read_key
	    (fptr, TDOUBLE, "RFIT1", &r_afit1, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword RFIT1\n");
	  status = 0;
	}
	else
	  {
	    afit[1] = r_afit1;
	  }
	if (fits_read_key
	    (fptr, TFLOAT, "REFMAG", &rfmag, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword REFMAG\n");
	  status = 0;
	}
	else
	  {
	    refmag = rfmag;
	  }
	if (fits_read_key
	    (fptr, TDOUBLE, "PSCALE", &pscale, comment, &status)) {
	  fits_report_error(stderr, status);
	  printf(" keyword PSCALE\n");
	  status = 0;
	}
	else
	  {
	    scalefactor = pscale;
	  }
      }

    /*
     * move to the next HDU 
     */
    hdunum = 3;
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status))
        printerror(status);

    if (hdutype != BINARY_TBL) {
        printf(" not a binary table \n");
        exit(2);
    }


    if (fits_read_key(fptr, TSTRING, "TFORM19", tform19, comment, &status)) {
        fits_report_error(stderr, status);
        printf(" keyword TFORM19\n");
        return (status);
    }

    maxlen = 0;
    vlen = strtok(tform19, "(");
    if (vlen != NULL) {
        vlen = strtok(NULL, ")");
        if (vlen != NULL) {
            maxlen = atoi(vlen);
            printf(" maximum column size = %d \n", maxlen);
        }
    }
    if (maxlen == 0) {
        printf(" error reading variable-length column size \n");
        printf(" keyword %s \n", tform19);
	exit(2);
    }
    colmaxlen[0] = maxlen;

    fits_get_num_cols(fptr, &ncol, &status);
    if (ncol != 20) {
        fprintf(stderr, " error in table format, no. columns = %d \n", ncol);
        exit(2);
    }

    fits_get_num_rows(fptr, &nobj, &status);

    printf(" %ld objects in table \n", nobj);


    if (fits_close_file(fptr, &status))
        printerror(status);

    tablesize[0] = (int) nume;
    tablesize[1] = (int) nobj;

    return status;
}


int
lensfits_read(char *tablename, int nume, int nobj, long *nmaximages, int *numexposures,
              char ***imagename, int **chip, float *e1, float *e2, double *wcsx,
              double *wcsy, double **psfe, double **psfe1, double **psfe2, int **fflag, 
              double *strehl, int *fitclass,
              double *size, double *bfrac, double *flux, double *snratio,
              double *profileSN, int *sampling, double **probvals, int *numvals,
              int **evalues, float **posterior, float *mag, int *id)
{

    fitsfile *fptr;             /* pointer to the FITS file, defined in fitsio.h */
    int status, hdunum, hdutype, anynull, intnull, typecode;
    long frow, felem, row, repeat, width, offset;
    long nume2, nobj2;
    float floatnull;
    int ncol, i, j, jj;
    int found, nfound;
    char **expname;
    double *temppsfe1, *temppsfe2;
    char delimiter[2];
    char *delimiterposn;
    int comparisonlen;

    // delimiter string that splits off the exposure name from the chip number in an image filename
    strcpy(delimiter,"_");

    status = 0;

    if (fits_open_file(&fptr, tablename, READONLY, &status))
        printerror(status);

    /*
     * move to the HDU 
     */
    hdunum = 2;
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status))
        printerror(status);

    if (hdutype != BINARY_TBL) {
        printf(" not a binary table \n");
        exit(2);
    }

    fits_get_num_cols(fptr, &ncol, &status);
    if (ncol != 2) {
        fprintf(stderr, " error in table format, no. columns = %d \n", ncol);
        exit(2);
    }

    fits_get_num_rows(fptr, &nume2, &status);

    if ((int) nume2 != nume) {
        printf(" unexpected table size! \n");
        printf(" expected nume = %d, got %d \n", nume, (int) nume2);
        exit(2);
    }

    frow = 1;
    felem = 1;
    floatnull = 0.;

    /*
     * read the columns 
     */
    fits_read_col(fptr, TFLOAT, 1, frow, felem, nume, &floatnull, e1,
                  &anynull, &status);
    fits_read_col(fptr, TFLOAT, 2, frow, felem, nume, &floatnull, e2,
                  &anynull, &status);

    // printf(" %f %f %f %f \n",e1[0],e1[nume-1],e2[0],e2[nume-1]);

    /*
     * move to the next HDU 
     */
    hdunum = 3;
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status))
        printerror(status);

    if (hdutype != BINARY_TBL) {
        printf(" not a binary table \n");
        exit(2);
    }

    fits_get_num_cols(fptr, &ncol, &status);
    if (ncol != 20) {
        fprintf(stderr, " error in table format, no. columns = %d \n", ncol);
        exit(2);
    }

    fits_get_num_rows(fptr, &nobj2, &status);

    if ((int) nobj2 != nobj) {
        printf(" unexpected table size! \n");
        printf(" expected nobj = %d, got %d \n", nobj, (int) nobj2);
        exit(2);
    }

    frow = 1;
    felem = 1;
    intnull = -4;
    floatnull = 0.;

    /*
     * read the columns 
     */
    if (fits_read_col(fptr, TINT, 1, frow, felem, nobj2, &anynull, id, 
		      &anynull, &status)) printerror(status);
    if (fits_read_col(fptr, TDOUBLE, 2, frow, felem, nobj2, &floatnull, wcsx,
		      &anynull, &status)) printerror(status);
    if (fits_read_col(fptr, TDOUBLE, 3, frow, felem, nobj2, &floatnull, wcsy,
		      &anynull, &status)) printerror(status);
    if (fits_read_col(fptr, TFLOAT, 4, frow, felem, nobj2, &floatnull, mag,
		      &anynull, &status)) printerror(status);

    typecode=0;
    fits_get_coltype(fptr, 5, &typecode, &repeat, &width, &status);
    if (typecode != 16) {
        fflush(stdout);
        fprintf(stderr,
                "error in table format,  column 5 is not a TSTRING: typecode %d \n",
                typecode);
        exit(2);
    }
    nmaximages[0] = repeat / width;
    printf(" %ld images found per object \n", nmaximages[0]);

    temppsfe1 = (double*)calloc(nmaximages[0], sizeof(double));
    temppsfe2 = (double*)calloc(nmaximages[0], sizeof(double));
    expname = (char**)calloc(nmaximages[0], sizeof(char*));
    nfound = 0;

    for (j=0; j<nmaximages[0]; j++)
      {
        expname[j] = (char *) calloc(width+1, sizeof(char));
	//printf (" %s h ",expname[j]);
      }

    for(i = 0; i < nobj; i++) {
        row = i + 1;
        imagename[i] = (char **) calloc(nmaximages[0], sizeof(char *));
        for(j = 0; j < nmaximages[0]; j++)
	  {
            imagename[i][j] = (char *) calloc(width+1, sizeof(char));
	    //printf (" %s h ",imagename[i][j]);
	  }
        if (fits_read_col(fptr, TSTRING, 5, row, 1, nmaximages[0], &anynull,
                          imagename[i], &anynull, &status))
            printerror(status);
	numexposures[i] = 0;
	
        for(j = 0; j < nmaximages[0]; j++)
	  {
	    //	    printf("after fitsread %s \n",imagename[i][j]);
	    if (strlen(imagename[i][j])>0) 
	      {
		if (strcmp(imagename[i][j]," ")!=0)
		  numexposures[i]++;
	      }
	  }

        psfe1[i] = (double*)calloc(nmaximages[0], sizeof(double));
        psfe2[i] = (double*)calloc(nmaximages[0], sizeof(double));
        // flag to identify if an exposure has been found or not for this object                   
        fflag[i] = (int*)calloc(nmaximages[0], sizeof(int));

	//printf(" numexp %d \n",numexposures[i]);

	if (numexposures[i]>0)
	  {
	    chip[i] = (int *)calloc(numexposures[i], sizeof(int));
	    if (fits_read_col(fptr, TINT, 6, row, 1, numexposures[i], &anynull,
			      chip[i], &anynull, &status))
	      printerror(status);
            if (fits_read_col(fptr, TDOUBLE, 7, row, 1, numexposures[i], &anynull,
                              temppsfe1, &anynull, &status))
              printerror(status);
            if (fits_read_col(fptr, TDOUBLE, 8, row, 1, numexposures[i], &anynull,
                              temppsfe2, &anynull, &status))
              printerror(status);

            // identify which exposure this is                                           
	  
	    // OLD HERE
	     if (write_psfs==0)
	      {
		for (j = 0; j<numexposures[i]; j++)
		  { // check each exposure in turn against a master list                              
		    //printf(" here %d %d %d %s\n",nobj,i,j,imagename[i][j]);
		    
		    found = 0;
		    for (jj=0; jj<nfound; jj++)
		      { // check if this image name is already in the master list                      
			//printf(" here %d %d %s %f %f \n",i,jj,expname[jj],temppsfe1[j],temppsfe2[j]); 

			comparisonlen = strlen(expname[jj]);
			if (strncmp(imagename[i][j],expname[jj],comparisonlen)==0)
			  {
			    found = 1;
			    fflag[i][jj]=1;
			    psfe1[i][jj] = temppsfe1[j];
			    psfe2[i][jj] = temppsfe2[j];
			    // printf(" found, %d %d %s %f %f \n",i+1,jj,expname[jj],temppsfe1[j],temppsfe2[j]);  
			    break;
			  }
		      }
		    if (found==0)
		      { // imagename doesn't match one previously found, add to master list            
			// first check that number of images has not been exceeded
			if (nfound >= nmaximages[0])
			  {
			    fprintf(stderr," error, too many exposures, %d found at object %d\n",nfound+1, i+1);
			    exit(EXIT_FAILURE);
			  }
			comparisonlen = strlen(imagename[i][j]);
			delimiterposn = strrchr(imagename[i][j], *delimiter);
			if (delimiterposn != NULL)
			  comparisonlen -= strlen(delimiterposn);
			strncpy(expname[nfound],imagename[i][j],comparisonlen);
			psfe1[i][nfound] = temppsfe1[j];
			psfe2[i][nfound] = temppsfe2[j];
			fflag[i][nfound]=1;
			// printf(" new, %d %d %s %f %f \n",i+1,nfound,expname[nfound],temppsfe1[j],temppsfe2[j]); 
			nfound++;
		      }
		  }
		  } 
	     // HACK HERE TO GET SOMETHING OUT WHEN ALL 5 EXPOSURES USED
	     
	    // HACK here

	     /* if (write_psfs==0)
	      {
		//printf(" %d %d numexp\n",numexposures[i],nmaximages[0]);
		if (numexposures[i]==nmaximages[0])
		  {

		    // printf(" %d yep\n",numexposures[i],nmaximages[0]);
		    // doesn't matter order can just write out all

		    for (jj = 0; jj<numexposures[i]; jj++)
		      {
			fflag[i][jj]=1;
			psfe1[i][jj] = temppsfe1[jj];
			psfe2[i][jj] = temppsfe2[jj];
			//printf(" %d %d %s %f %f \n",i,jj,expname[jj],temppsfe1[jj],temppsfe2[jj]);  
		      }
		  }
	      }*/


	      }
	else
	  {
	    chip[i] = (int *)calloc(1, sizeof(int));
            psfe1[i] = (double *)calloc(1, sizeof(double));
            psfe2[i] = (double *)calloc(1, sizeof(double));

	  }

        if (fits_read_col(fptr, TDOUBLE, 9, row, 1, 2, &floatnull, psfe[i],
                          &anynull, &status))
            printerror(status);
	//printf(" got exposures\n");fflush(stdout);

    }

    //printf(" got exposures\n");fflush(stdout);

    fits_read_col(fptr, TDOUBLE, 10, frow, felem, nobj2, &floatnull, strehl,
                  &anynull, &status);
    fits_read_col(fptr, TINT, 11, frow, felem, nobj2, &intnull, fitclass,
                  &anynull, &status);
    fits_read_col(fptr, TDOUBLE, 12, frow, felem, nobj2, &floatnull, size,
                  &anynull, &status);
    fits_read_col(fptr, TDOUBLE, 13, frow, felem, nobj2, &floatnull, bfrac,
                  &anynull, &status);
    fits_read_col(fptr, TDOUBLE, 14, frow, felem, nobj2, &floatnull, flux,
                  &anynull, &status);
    fits_read_col(fptr, TDOUBLE, 15, frow, felem, nobj2, &floatnull, snratio,
                  &anynull, &status);
    fits_read_col(fptr, TDOUBLE, 16, frow, felem, nobj2, &floatnull, profileSN,
                  &anynull, &status);
    if (fits_read_col(fptr, TINT, 17, frow, felem, nobj2, &intnull, sampling,
		      &anynull, &status)) printerror(status);

    for(i = 0; i < nobj; i++) {
        row = i + 1;
        if (fits_read_col(fptr, TDOUBLE, 18, row, 1, 2, &floatnull, probvals[i],
                          &anynull, &status))
            printerror(status);
    }

    for(i = 0; i < nobj; i++) {
        row = i + 1;
        if (fits_read_descript(fptr, 19, row, &repeat, &offset, &status))
            printerror(status);
        if (repeat > 0) {
            if (fits_read_col
                (fptr, TINT, 19, row, 1, repeat, &intnull, evalues[i],
                 &anynull, &status))
                printerror(status);
            if (fits_read_col
                (fptr, TFLOAT, 20, row, 1, repeat, &floatnull, posterior[i],
                 &anynull, &status))
                printerror(status);
        }
        numvals[i] = repeat;
    }

    //HACK REMOVAL HERE
    /*   if (write_psfs==0)
      {
	for (j=0; j<nfound; j++)
	  {
	    printf (" %s ",expname[j]);
	  }
	printf("\n");
	if (nfound != nmaximages[0])
	  {
	    printf ("error: found %d individual exposures \n",nfound);
	    printf (" continuing without writing individual PSF values\n");
	    write_psfs = -1;
	    fflush(stdout);
	  }
	  }*/

    if (fits_close_file(fptr, &status))
        printerror(status);

    //printf (" table read OK\n");

    return status;
}

void
printerror(int status)
{
    /*****************************************************/
    /*
     * Print out cfitsio error messages and exit program 
     */
    /*****************************************************/

    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error report */
        exit(status);           /* terminate the program, returning error status */
    }
    return;
}


int
writeposterior_fits(char *fitsname, float *f, int nobj, int pwidth,
                    int pheight)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[3], fpixel[3] = { 1, 1, 1 };
    int bitpix;
    int psfsize, psfsize3d;

    /*
     * write out psf data as a fits file 
     */

    if (pwidth != pheight) {
        fprintf(stderr, " PSF subimages must be square! \n");
        return 1;
    }

    if (pwidth <= 0 || pheight <= 0 || nobj <= 0) {
        fprintf(stderr, " error in subimage dimensions or number \n");
        return 1;
    }

    bitpix = -32;
    anaxis = 3;
    anaxes[0] = pwidth;
    anaxes[1] = pheight;
    anaxes[2] = nobj;

    psfsize = anaxes[0] * anaxes[1];
    psfsize3d = psfsize * nobj;

    printf(" dimensions: %d %d %d \n", pwidth, pheight, nobj);

    fits_create_file(&afptr, fitsname, &status);
    fits_create_img(afptr, bitpix, anaxis, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * write all PSF data into image array 
     */

    if (fits_write_pix(afptr, TFLOAT, fpixel, psfsize3d, f, &status)) {
        printf(" error reading pixel data \n");
        exit(2);
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

double rfunc (double r, float mag)
{
  // prior function for scalelength as a function of magnitude, derived from
  //   Simard et al 2002, assumes quasi-rayleigh distribution in r.  
  
  // The parameters in this function are afit[] and refmag, which are defined and set 
  // as global variables in front of the main code
  
  // NB returns ln(prior) 

  // parameters assume r is measured in arcsec

  double aval, rprior,tmag;

  if (r<=0.)
    {
      fflush(stdout);
      fprintf(stderr," illegal value given to rfunc: r=%f \n",r);
      exit(EXIT_FAILURE);
    }

  // function of magnitude for parameter a
  tmag = mag-refmag;

  // linear fit to median, then rescale to get aval
  aval = exp(afit[0] + afit[1]*tmag)/0.83255461;
  
  // exponent (4/3), convert median to required scale 
  // normalised to infinity by factor 3.sqrt(pi)/8.
  aval /= 1.134;
  rprior = r*exp(-pow((r/aval),(4./3.)))/(3.*sqrt(pi)*aval*aval/8.);

  return log(rprior);

}
