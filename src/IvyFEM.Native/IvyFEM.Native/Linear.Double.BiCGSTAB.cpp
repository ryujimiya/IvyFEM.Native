#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"

bool DoubleSolvePreconditionedBiCGSTAB(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
    int LUIndexsLength, int* LUPtrs, int* LUIndexs, double* LUValues,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    double* r = new double[n];
    double* r0 = new double[n];
    double* p = new double[n];
    double* s = new double[n];
    double* Mp = new double[n];
    double* Ms = new double[n];
    double* AMp = new double[n];
    double* AMs = new double[n];
    double* dummy = new double[n];
#define __FINALIZE \
delete[] r; \
delete[] r0; \
delete[] p; \
delete[] s; \
delete[] Mp; \
delete[] Ms; \
delete[] AMp; \
delete[] AMs; \
delete[] dummy; \
fflush(stdout);

    int iter = 0;

    memset(X, 0, n * sizeof(double));
    memset(r, 0, n * sizeof(double));
    memset(r0, 0, n * sizeof(double));
    memset(p, 0, n * sizeof(double));
    memset(s, 0, n * sizeof(double));
    memset(Mp, 0, n * sizeof(double));
    memset(Ms, 0, n * sizeof(double));
    memset(dummy, 0, n * sizeof(double));

    memcpy(r, B, n * sizeof(double));

    double sqInvNormRes0;
    {
        double sqNormRes0 = DoubleDot(n, r, r);
        if (sqNormRes0 < PrecisionLowerLimit)
        {
            convRatio = 0;
            printf("iter = %d norm: %e\r\n", iter, convRatio);

            __FINALIZE
                return true;
        }
        sqInvNormRes0 = 1.0 / sqNormRes0;
    }

    memcpy(r0, r, n * sizeof(double));
    memcpy(p, r, n * sizeof(double));

    for (iter = 0; iter < maxIter; iter++)
    {
        DoubleSolveLU(Mp, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, p);

        double r0r = DoubleDot(n, r0, r);
        // AMp = A * Mp;
        DoubleMV(AMp,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, Mp,
            0.0, dummy);
        double alpha;
        {
            double denominator = DoubleDot(n, r0, AMp);
            alpha = r0r / denominator;
        }
        DoubleAxpy(s, -alpha, n, AMp, r);

        DoubleSolveLU(Ms, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, s);

        // AMs = A * Ms;
        DoubleMV(AMs,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, Ms,
            0.0, dummy);
        double omega;
        {
            double denominator = DoubleDot(n, AMs, AMs);
            double numerator = DoubleDot(n, s, AMs);
            omega = numerator / denominator;
        }

        DoubleAxpy(X, alpha, n, Mp, X);
        DoubleAxpy(X, omega, n, Ms, X);
        DoubleAxpy(r, -omega, n, AMs, s);

        {
            double sqNormRes = DoubleDot(n, r, r);
            if (sqNormRes * sqInvNormRes0 < tolerance * tolerance)
            {
                convRatio = sqrt(sqNormRes * sqInvNormRes0);
                printf("iter = %d norm: %e\r\n", iter, convRatio);

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
        DoubleAxpy(p, -beta * omega, n, AMp, p);
    }

    {
        double sqNormRes = DoubleDot(n, r, r);
        convRatio = sqrt(sqNormRes * sqInvNormRes0);
        printf("iter = %d norm: %e\r\n", iter, convRatio);
    }
    printf("Not converged\r\n");

    __FINALIZE
#undef __FINALIZE
        return false;
}

bool DoubleSolveBiCGSTAB(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B, int fillinLevel,
    double convRatioTolerance)
{
    ULONGLONG t;
    int LUIndexsLength = 0;
    int* LUPtrs = NULL;
    int* LUIndexs = NULL;
    double* LUValues = NULL;
    t = GetTickCount64();
    DoubleCalcILU(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues,
        n, AIndexsLength, APtrs, AIndexs, AValues, fillinLevel);
    printf("    1: t = %lld\r\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = DoubleSolvePreconditionedBiCGSTAB(
        X,
        n, AIndexsLength, APtrs, AIndexs, AValues, B,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\r\n", GetTickCount64() - t);
    fflush(stdout);
    return success;
}

bool DoubleSolveBiCGSTABWithPivoting(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B, int fillinLevel,
    double convRatioTolerance)
{
    ULONGLONG t;
    int LUIndexsLength = 0;
    int* LUPtrs = NULL;
    int* LUIndexs = NULL;
    double* LUValues = NULL;
    int* pivot = new int[n];
    int pivotingAIndexsLength = 0;
    int* pivotingAPtrs = new int[n + 1];
    int* pivotingAIndexs = new int[AIndexsLength];
    double* pivotingAValues = new double[AIndexsLength];
    double* pivotingB = new double[n];

    t = GetTickCount64();
    DoubleCalcILUWithPivoting(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues, pivot,
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
                double value = AValues[iPtr];
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
    bool success = DoubleSolvePreconditionedBiCGSTAB(
        X,
        n, pivotingAIndexsLength, pivotingAPtrs, pivotingAIndexs, pivotingAValues, pivotingB,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\r\n", GetTickCount64() - t);
    fflush(stdout);

    delete[] pivot;
    delete[] pivotingAPtrs;
    delete[] pivotingAIndexs;
    delete[] pivotingAValues;
    delete[] pivotingB;
    return success;
}
