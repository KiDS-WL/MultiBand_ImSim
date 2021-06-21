#include <math.h>
#define NRANSI
#include "nrutil.h"

double
dpythag(double a, double b)
{
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
        return absa * sqrt(1.0 +
                           (absb ==
                            0.0 ? 0.0 : (absb / absa) * (absb / absa)));
    else
        return (absb ==
                0.0 ? 0.0 : absb * sqrt(1.0 +
                                        (absa ==
                                         0.0 ? 0.0 : (absa / absb) * (absa /
                                                                      absb))));
}

#undef NRANSI
