#include <fitsio.h>

int
writepsfstamp(int nx, int ny, double *data, char *name)
{
    int status = 0;
    long naxes[2];
    fitsfile *fptr;

    naxes[0] = nx;
    naxes[1] = ny;
    fits_create_file(&fptr, name, &status);
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
    fits_write_img(fptr, TDOUBLE, 1, naxes[0] * naxes[1], data, &status);
    fits_close_file(fptr, &status);
    if (status)
        fits_report_error(stderr, status);
    return status;
}
