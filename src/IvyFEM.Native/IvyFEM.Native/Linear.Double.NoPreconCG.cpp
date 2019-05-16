#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"

bool DoubleSolveNoPreconCG(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    double* r = new double[n];
    double* p = new double[n];
    double* Ap = new double[n];
    double* dummy = new double[n];
#define __FINALIZE \
delete[] r; \
delete[] p; \
delete[] Ap; \
delete[] dummy; \
fflush(stdout);

    int iter = 0;

    memset(X, 0, n * sizeof(double));
    memset(r, 0, n * sizeof(double));
    memset(p, 0, n * sizeof(double));
    memset(dummy, 0, n * sizeof(double));

    memcpy(r, B, n * sizeof(double));

    double rr;
    double sqInvNormRes0;
    {
        double sqNormRes0 = DoubleDot(n, r, r);
        if (sqNormRes0 < PrecisionLowerLimit)
        {
            convRatio = 0;
            printf("iter = %d norm: %e\n", iter, convRatio);

            __FINALIZE
                return true;
        }
        rr = sqNormRes0;
        sqInvNormRes0 = 1.0 / sqNormRes0;
    }

    memcpy(p, r, n * sizeof(double));

    for (iter = 0; iter < maxIter; iter++)
    {
        // Ap = A * p;
        DoubleMV(Ap,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, p,
            0.0, dummy);

        double alpha;
        {
            double pAp = DoubleDot(n, p, Ap);
            alpha = rr / pAp;
        }
        DoubleAxpy(r, -alpha, n, Ap, r);
        DoubleAxpy(X, alpha, n, p, X);

        double rrPrev = rr;
        {
            double sqNormRes = DoubleDot(n, r, r);
            if (sqNormRes * sqInvNormRes0 < tolerance * tolerance)
            {
                convRatio = sqrt(sqNormRes * sqInvNormRes0);
                printf("iter = %d norm: %e\n", iter, convRatio);

                __FINALIZE
                    return true;
            }
            rr = sqNormRes;
        }
        double beta = rr / rrPrev;

        DoubleAxpy(p, beta, n, p, r);
    }

    {
        double sqNormRes = DoubleDot(n, r, r);
        convRatio = sqrt(sqNormRes * sqInvNormRes0);
        printf("iter = %d norm: %e\n", iter, convRatio);
    }
    printf("Not converged\n");

    __FINALIZE
#undef __FINALIZE
        return false;
}
