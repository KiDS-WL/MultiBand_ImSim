#ifndef __LENSUTIL_H__
#define __LENSUTIL_H__

void *ml_calloc(size_t nmemb, size_t size, float *mem, const char *format, ...);
void get_poly_order(int order, int chiporder, int chipvariation, int nchip,
		    int *ncoeffs, int *crossterm);

#endif

