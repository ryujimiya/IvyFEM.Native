#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"
#include "IvyFEM.Native.Internal.h"

bool ComplexSolveNoPreconBiCGSTAB(__complex* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    __complex* r = new __complex[n];
    __complex* r0 = new __complex[n];
    __complex* p = new __complex[n];
    __complex* s = new __complex[n];
    __complex* Ap = new __complex[n];
    __complex* As = new __complex[n];
    __complex* dummy = new __complex[n];
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

    memset(X, 0, n * sizeof(__complex));
    memset(r, 0, n * sizeof(__complex));
    memset(r0, 0, n * sizeof(__complex));
    memset(p, 0, n * sizeof(__complex));
    memset(s, 0, n * sizeof(__complex));
    memset(dummy, 0, n * sizeof(__complex));

    memcpy(r, B, n * sizeof(__complex));

    double sqInvNormRes0;
    {
        double sqNormRes0 = std::real(__ComplexDotc(n, r, r));
        if (sqNormRes0 < PrecisionLowerLimit)
        {
            convRatio = 0;
            printf("iter = %d norm: %e\n", iter, convRatio);

            __FINALIZE
                return true;
        }
        sqInvNormRes0 = 1.0 / sqNormRes0;
    }

    memcpy(r0, r, n * sizeof(__complex));
    memcpy(p, r, n * sizeof(__complex));
    __complex r0r = __ComplexDotc(n, r0, r);

    for (iter = 0; iter < maxIter; iter++)
    {
        // Ap = A * p;
        ComplexMV(Ap,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, p,
            0.0, dummy);
        __complex alpha;
        {
            __complex denominator = __ComplexDotc(n, r0, Ap);
            alpha = r0r / denominator;
        }
        ComplexAxpy(s, -alpha, n, Ap, r);

        // As = A * s;
        ComplexMV(As,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, s,
            0.0, dummy);
        __complex omega;
        {
            __complex denominator = __ComplexDotc(n, As, As);
            __complex numerator = __ComplexDotc(n, s, As);
            omega = numerator / denominator;
        }

        ComplexAxpy(X, alpha, n, p, X);
        ComplexAxpy(X, omega, n, s, X);
        ComplexAxpy(r, -omega, n, As, s);

        {
            double sqNormRes = std::real(__ComplexDotc(n, r, r));
            if (sqNormRes * sqInvNormRes0 < tolerance * tolerance)
            {
                convRatio = sqrt(sqNormRes * sqInvNormRes0);
                printf("iter = %d norm: %e\n", iter, convRatio);

                __FINALIZE
                    return true;
            }
        }

        __complex beta;
        {
            __complex r0rPrev = r0r;
            r0r = __ComplexDotc(n, r0, r);
            beta = (r0r * alpha) / (r0rPrev * omega);
        }
        ComplexAxpy(p, beta, n, p, r);
        ComplexAxpy(p, -beta * omega, n, Ap, p);
    }

    {
        double sqNormRes = std::real(__ComplexDotc(n, r, r));
        convRatio = sqrt(sqNormRes * sqInvNormRes0);
        printf("iter = %d norm: %e\n", iter, convRatio);
    }
    printf("Not converged\n");

    __FINALIZE
#undef __FINALIZE
        return false;
}
