#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"
#include "IvyFEM.Native.Internal.h"

bool ComplexSolvePreconditionedCOCG(__complex* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B,
    int LUIndexsLength, int* LUPtrs, int* LUIndexs, __complex* LUValues,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    __complex* r = new __complex[n];
    __complex* z = new __complex[n];
    __complex* p = new __complex[n];
    __complex* Ap = new __complex[n];
    __complex* dummy = new __complex[n];
#define __FINALIZE \
delete[] r; \
delete[] z; \
delete[] p; \
delete[] Ap; \
delete[] dummy;

    int iter = 0;

    memset(X, 0, n * sizeof(__complex));
    memset(r, 0, n * sizeof(__complex));
    memset(z, 0, n * sizeof(__complex));
    memset(p, 0, n * sizeof(__complex));
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

    ComplexSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, r);

    memcpy(p, z, n * sizeof(__complex));
    __complex rz = __ComplexDotu(n, r, z);

    for (iter = 0; iter < maxIter; iter++)
    {
        // Ap = A * p;
        ComplexMV(Ap,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, p,
            0.0, dummy);

        __complex alpha;
        {
            __complex pAp = __ComplexDotu(n, p, Ap);
            alpha = rz / pAp;
        }
        ComplexAxpy(r, -alpha, n, Ap, r);
        ComplexAxpy(X, alpha, n, p, X);

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

        ComplexSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, r);

        __complex rzPrev = rz;
        rz = __ComplexDotu(n, r, z);
        __complex beta = rz / rzPrev;

        ComplexAxpy(p, beta, n, p, z);
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

bool ComplexSolveCOCG(__complex* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B, int fillinLevel,
    double convRatioTolerance)
{
    ULONGLONG t;
    int LUIndexsLength = 0;
    int* LUPtrs = NULL;
    int* LUIndexs = NULL;
    __complex* LUValues = NULL;
    t = GetTickCount64();
    ComplexCalcILU(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues,
        n, AIndexsLength, APtrs, AIndexs, AValues, fillinLevel);
    printf("    1: t = %lld\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = ComplexSolvePreconditionedCOCG(
        X,
        n, AIndexsLength, APtrs, AIndexs, AValues, B,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    ComplexDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\n", GetTickCount64() - t);
    fflush(stdout);
    return success;
}
