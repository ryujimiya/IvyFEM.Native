#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"
#include "IvyFEM.Native.Internal.h"

bool ComplexSolvePreconditionedBiCGSTAB(__complex * X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex * AValues, __complex * B,
    int LUIndexsLength, int* LUPtrs, int* LUIndexs, __complex * LUValues,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    __complex* r = new __complex[n];
    __complex* r0 = new __complex[n];
    __complex* p = new __complex[n];
    __complex* s = new __complex[n];
    __complex* Mp = new __complex[n];
    __complex* Ms = new __complex[n];
    __complex* AMp = new __complex[n];
    __complex* AMs = new __complex[n];
    __complex* z = new __complex[n];
    __complex* dummy = new __complex[n];
#define __FINALIZE \
delete[] r; \
delete[] r0; \
delete[] p; \
delete[] s; \
delete[] Mp; \
delete[] Ms; \
delete[] AMp; \
delete[] AMs; \
delete[] z; \
delete[] dummy; \
fflush(stdout);

    int iter = 0;

    memset(X, 0, n * sizeof(__complex));
    memset(r, 0, n * sizeof(__complex));
    memset(r0, 0, n * sizeof(__complex));
    memset(p, 0, n * sizeof(__complex));
    memset(s, 0, n * sizeof(__complex));
    memset(Mp, 0, n * sizeof(__complex));
    memset(Ms, 0, n * sizeof(__complex));
    memset(z, 0, n * sizeof(__complex));
    memset(dummy, 0, n * sizeof(__complex));

    memcpy(r, B, n * sizeof(__complex));

    double sqInvNormRes0;
    {
        double sqNormRes0 = std::real(__ComplexDotc(n, r, r));
        if (sqNormRes0 < PrecisionLowerLimit)
        {
            convRatio = 0;
            printf("iter = %d norm: %e\r\n", iter, convRatio);

            __FINALIZE
                return true;
        }
        sqInvNormRes0 = 1.0 / sqNormRes0;
    }

    memcpy(r0, r, n * sizeof(__complex));
    memcpy(p, r, n * sizeof(__complex));

    for (iter = 0; iter < maxIter; iter++)
    {
        ComplexSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, p);
        memcpy(Mp, z, n * sizeof(double));

        __complex r0r = __ComplexDotc(n, r0, r);

        // AMp = A * Mp;
        ComplexMV(AMp,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, Mp,
            0.0, dummy);
        __complex alpha;
        {
            __complex denominator = __ComplexDotc(n, r0, AMp);
            alpha = r0r / denominator;
        }
        ComplexAxpy(s, -alpha, n, AMp, r);

        ComplexSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, s);
        memcpy(Ms, z, n * sizeof(double));

        // AMs = A * Ms;
        ComplexMV(AMs,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, Ms,
            0.0, dummy);
        __complex omega;
        {
            __complex denominator = __ComplexDotc(n, AMs, AMs);
            __complex numerator = __ComplexDotc(n, s, AMs);
            omega = numerator / denominator;
        }

        ComplexAxpy(X, alpha, n, Mp, X);
        ComplexAxpy(X, omega, n, Ms, X);
        ComplexAxpy(r, -omega, n, AMs, s);

        {
            double sqNormRes = std::real(__ComplexDotc(n, r, r));
            if (sqNormRes * sqInvNormRes0 < tolerance * tolerance)
            {
                convRatio = sqrt(sqNormRes * sqInvNormRes0);
                printf("iter = %d norm: %e\r\n", iter, convRatio);

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
        ComplexAxpy(p, -beta * omega, n, AMp, p);
    }

    {
        double sqNormRes = std::real(__ComplexDotc(n, r, r));
        convRatio = sqrt(sqNormRes * sqInvNormRes0);
        printf("iter = %d norm: %e\r\n", iter, convRatio);
    }
    printf("Not converged\r\n");

    __FINALIZE
#undef __FINALIZE
        return false;
}

#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"

bool ComplexSolveBiCGSTAB(__complex* X,
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
    printf("    1: t = %lld\r\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = ComplexSolvePreconditionedBiCGSTAB(
        X,
        n, AIndexsLength, APtrs, AIndexs, AValues, B,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    ComplexDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\r\n", GetTickCount64() - t);
    fflush(stdout);
    return success;
}

bool ComplexSolveBiCGSTABWithPivoting(__complex* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B, int fillinLevel,
    double convRatioTolerance)
{
    ULONGLONG t;
    int LUIndexsLength = 0;
    int* LUPtrs = NULL;
    int* LUIndexs = NULL;
    __complex* LUValues = NULL;
    int* pivot = new int[n];
    int pivotingAIndexsLength = 0;
    int* pivotingAPtrs = new int[n + 1];
    int* pivotingAIndexs = new int[AIndexsLength];
    __complex* pivotingAValues = new __complex[AIndexsLength];
    __complex* pivotingB = new __complex[n];

    t = GetTickCount64();
    ComplexCalcILUWithPivoting(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues, pivot,
        n, AIndexsLength, APtrs, AIndexs, AValues, fillinLevel);
    {
        int pos = 0;
        for (int row = 0; row < n; row++)
        {
            int oldRow = pivot[row];
            int sPtr = APtrs[oldRow];
            int ePtr = APtrs[oldRow + 1];
            pivotingAPtrs[row] = pos;
            for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
            {
                int col = AIndexs[iPtr];
                __complex value = AValues[iPtr];
                pivotingAIndexs[pos] = col;
                pivotingAValues[pos] = value;
                pos++;
            }
        }
        pivotingAPtrs[n] = pos;
        pivotingAIndexsLength = pos;
    }
    for (int i = 0; i < n; i++)
    {
        pivotingB[i] = B[pivot[i]];
    }
    printf("    1: t = %lld\r\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = ComplexSolvePreconditionedBiCGSTAB(
        X,
        n, pivotingAIndexsLength, pivotingAPtrs, pivotingAIndexs, pivotingAValues, pivotingB,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    ComplexDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\r\n", GetTickCount64() - t);
    fflush(stdout);

    delete[] pivot;
    delete[] pivotingAPtrs;
    delete[] pivotingAIndexs;
    delete[] pivotingAValues;
    delete[] pivotingB;
    return success;
}


