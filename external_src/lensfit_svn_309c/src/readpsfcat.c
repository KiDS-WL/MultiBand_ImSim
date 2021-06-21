#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

int
readpsfcat(char *catname, float *objx, float *objy, float *mag)
{
    FILE *catdata;
    int n;
    char string[2000];
    char delims[] = " \t";
    char *item = NULL;


    catdata = fopen(catname, "r");
    if (catdata == NULL) {
        printf(" failed to open file %s \n", catname);
        exit(2);
    }

    n = 0;

    while (fgets(string, 2000, catdata) != NULL) {
        if (strncmp(string, "#", 1) != 0) {
            item = strtok(string, delims);
            if (item == NULL) {
                fprintf(stderr,
                        " error reading 1st column (x) from catalogue \n");
                exit(2);
            }
            objx[n] = atof(item);
            item = strtok(NULL, delims);
            if (item == NULL) {
                fprintf(stderr,
                        " error reading 2nd column (y) from catalogue \n");
                exit(2);
            }
            objy[n] = atof(item);
            /*
             * magnitude is an optional item in the catalogue 
             */
            item = strtok(NULL, delims);
            if (item == NULL) {
                mag[n] = 0.;
            } else {
                mag[n] = (float) atof(item);
            }
            n++;
        }
    }

    fclose(catdata);

    return n;

}
