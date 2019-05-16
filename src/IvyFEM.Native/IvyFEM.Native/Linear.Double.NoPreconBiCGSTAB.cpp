#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"

bool DoubleSolveNoPreconBiCGSTAB(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    double* r = new double[n];
    double* r0 = new double[n];
    double* p = new double[n];
    double* s = new double[n];
    double* Ap = new double[n];
    double* As = new double[n];
    double* dummy = new double[n];
#define __FINALIZE \
delete[] r; \
delete[] r0; \
delete[] p; \
delete[] s; \
delete[] Ap; \
delete[] As; \
delete[] dummy; \
fflush(stdout);

    int iter = 0;

    memset(X, 0, n * sizeof(double));
    memset(r, 0, n * sizeof(double));
    memset(r0, 0, n * sizeof(double));
    memset(p, 0, n * sizeof(double));
    memset(s, 0, n * sizeof(double));
    memset(dummy, 0, n * sizeof(double));

    memcpy(r, B, n * sizeof(double));

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
        sqInvNormRes0 = 1.0 / sqNormRes0;
    }

    memcpy(r0, r, n * sizeof(double));
    memcpy(p, r, n * sizeof(double));
    double r0r = DoubleDot(n, r0, r);

    for (iter = 0; iter < maxIter; iter++)
    {
        // Ap = A * p;
        DoubleMV(Ap,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, p,
            0.0, dummy);
        double alpha;
        {
            double denominator = DoubleDot(n, r0, Ap);
            alpha = r0r / denominator;
        }
        DoubleAxpy(s, -alpha, n, Ap, r);

        // As = A * s;
        DoubleMV(As,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, s,
            0.0, dummy);
        double omega;
        {
            double denominator = DoubleDot(n, As, As);
            double numerator = DoubleDot(n, s, As);
            omega = numerator / denominator;
        }

        DoubleAxpy(X, alpha, n, p, X);
        DoubleAxpy(X, omega, n, s, X);
        DoubleAxpy(r, -omega, n, As, s);

        {
            double sqNormRes = DoubleDot(n, r, r);
            if (sqNormRes * sqInvNormRes0 < tolerance * tolerance)
            {
                convRatio = sqrt(sqNormRes * sqInvNormRes0);
                printf("iter = %d norm: %e\n", iter, convRatio);

                __FINALIZE
                    return true;
            }
        }

        double beta;
        {
            double r0rPrev = r0r;
            r0r = DoubleDot(n, r0, r);
            beta = (r0r * alpha) / (r0rPrev * omega);
        }
        DoubleAxpy(p, beta, n, p, r);
        DoubleAxpy(p, -beta * omega, n, Ap, p);
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
