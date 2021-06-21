#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

int
readpsfcatsize(char *catname)
{
    FILE *catdata;
    char string[2000];
    int nobj;

    catdata = fopen(catname, "r");
    if (catdata == NULL) {
        printf(" failed to open file %s \n", catname);
        exit(2);
    }

    nobj = 0;

    while (fgets(string, 2000, catdata) != NULL) {
        if (strncmp(string, "#", 1) != 0) {
            nobj++;
        }
    }

    fclose(catdata);

    return nobj;

}
