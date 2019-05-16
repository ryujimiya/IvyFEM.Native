#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"

bool DoubleSolvePreconditionedCG(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
    int LUIndexsLength, int* LUPtrs, int* LUIndexs, double* LUValues,
    double convRatioTolerance)
{
    double convRatio = convRatioTolerance;
    double tolerance = convRatio;
    int maxIter = MaxIter;
    double* r = new double[n];
    double* z = new double[n];
    double* p = new double[n];
    double* Ap = new double[n];
    double* dummy = new double[n];
#define __FINALIZE \
delete[] r; \
delete[] z; \
delete[] p; \
delete[] Ap; \
delete[] dummy;

    int iter = 0;

    memset(X, 0, n * sizeof(double));
    memset(r, 0, n * sizeof(double));
    memset(z, 0, n * sizeof(double));
    memset(p, 0, n * sizeof(double));
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

    DoubleSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, r);

    memcpy(p, z, n * sizeof(double));
    double rz = DoubleDot(n, r, z);

    for (iter = 0; iter < maxIter; iter++)
    {
        // Ap = A * p;
        DoubleMV(Ap,
            1.0, n, AIndexsLength, APtrs, AIndexs, AValues, p,
            0.0, dummy);

        double alpha;
        {
            double pAp = DoubleDot(n, p, Ap);
            alpha = rz / pAp;
        }
        DoubleAxpy(r, -alpha, n, Ap, r);
        DoubleAxpy(X, alpha, n, p, X);

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

        DoubleSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, r);

        double rzPrev = rz;
        rz = DoubleDot(n, r, z);
        double beta = rz / rzPrev;

        DoubleAxpy(p, beta, n, p, z);
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

bool DoubleSolveCG(double* X,
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
    printf("    1: t = %lld\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = DoubleSolvePreconditionedCG(
        X,
        n, AIndexsLength, APtrs, AIndexs, AValues, B,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\n", GetTickCount64() - t);
    fflush(stdout);
    return success;
}

bool DoubleSolveCGWithPivoting(double* X,
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
    printf("    1: t = %lld\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = DoubleSolvePreconditionedCG(
        X,
        n, pivotingAIndexsLength, pivotingAPtrs, pivotingAIndexs, pivotingAValues, pivotingB,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\n", GetTickCount64() - t);
    fflush(stdout);

    delete[] pivot;
    delete[] pivotingAPtrs;
    delete[] pivotingAIndexs;
    delete[] pivotingAValues;
    delete[] pivotingB;
    return success;
}

bool DoubleSolveICCG(double* X,
    int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
    double convRatioTolerance)
{
    ULONGLONG t;
    int LUIndexsLength = 0;
    int* LUPtrs = NULL;
    int* LUIndexs = NULL;
    double* LUValues = NULL;
    t = GetTickCount64();
    DoubleCalcIC(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues,
        n, AIndexsLength, APtrs, AIndexs, AValues);
    printf("    1: t = %lld\n", GetTickCount64() - t);
    t = GetTickCount64();
    bool success = DoubleSolvePreconditionedCG(
        X,
        n, AIndexsLength, APtrs, AIndexs, AValues, B,
        LUIndexsLength, LUPtrs, LUIndexs, LUValues,
        convRatioTolerance);
    DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
    printf("    2: t = %lld\n", GetTickCount64() - t);
    fflush(stdout);
    return success;
}
