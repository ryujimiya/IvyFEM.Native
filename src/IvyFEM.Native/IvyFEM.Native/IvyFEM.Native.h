#pragma once

#ifdef IVYFEMNATIVESTUB_EXPORTS
#define DLL_API __declspec(dllexport)
#else
#define DLL_API __declspec(dllimport)
#endif

#include <complex>
typedef std::complex<double> __complex;

extern "C"
{
	//////////////////////////////////////////////////////////////////////////
	// DoubleLinear
	DLL_API double DoubleDot(int n, double *X, double *Y);

	DLL_API void DoubleMV(double *Z,
		double alpha, int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *X,
		double beta, double *Y);

	DLL_API void DoubleAxpy(double *Z, double alpha, int n, double *X, double *Y);

	DLL_API bool DoubleSolveNoPreconCG(double *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
		double convRatioTolerance);

	DLL_API void DoubleDeleteCSR(int *APtrsP, int *AIndexsP, double *AValuesP);

	// Note: 内部でメモリ割り当てを行っているので使用後DeleteCSRを呼び出す必要がある
	DLL_API void DoubleCalcILU(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, int fillinLevel);

	DLL_API void DoubleSolveLU(
		double *X, int n, int LUIndexsLength, int *LUPtr, int *LUIndexs, double *LUValues, double *B);

	DLL_API void DoubleCalcILUWithPivoting(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
		int *pivot, int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, int fillinLevel);

	DLL_API bool DoubleSolveCGWithPivoting(double *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B, int fillinLevel,
		double convRatioTolerance);

	DLL_API void DoubleCalcIC(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues);

    ////////////////////////////////////////////////////
    // DoubleCG
    DLL_API bool DoubleSolvePreconditionedCG(double* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
        int LUIndexsLength, int* LUPtrs, int* LUIndexs, double* LUValues,
        double convRatioTolerance);

    DLL_API bool DoubleSolveCG(double* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B, int fillinLevel,
        double convRatioTolerance);

    DLL_API bool DoubleSolveICCG(double *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
		double convRatioTolerance);
		
    ////////////////////////////////////////////////////
    // DoubleBiCGSTAB
    DLL_API bool DoubleSolveNoPreconBiCGSTAB(double *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
		double convRatioTolerance);

    DLL_API bool DoubleSolvePreconditionedBiCGSTAB(double* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B,
        int LUIndexsLength, int* LUPtrs, int* LUIndexs, double* LUValues,
        double convRatioTolerance);

    DLL_API bool DoubleSolveBiCGSTAB(double* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B, int fillinLevel,
        double convRatioTolerance);

    DLL_API bool DoubleSolveBiCGSTABWithPivoting(double* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, double* AValues, double* B, int fillinLevel,
        double convRatioTolerance);

    //////////////////////////////////////////////////////////////////////////
	// ComplexLinear
	DLL_API void ComplexDotc(__complex *ret, int n, __complex *X, __complex *Y);

	DLL_API void ComplexDotu(__complex *ret, int n, __complex *X, __complex *Y);

	DLL_API void ComplexMV(__complex *Z,
		__complex alpha, int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *X,
		__complex beta, __complex *Y);

	DLL_API void ComplexAxpy(__complex *Z, __complex alpha, int n, __complex *X, __complex *Y);

	DLL_API bool ComplexSolveNoPreconCOCG(__complex *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *B,
		double convRatioTolerance);

	DLL_API void ComplexDeleteCSR(int *APtrsP, int *AIndexsP, __complex *AValuesP);

	// Note: 内部でメモリ割り当てを行っているので使用後DeleteCSRを呼び出す必要がある
	DLL_API void ComplexCalcILU(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, __complex **LUValuesP,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, int fillinLevel);

	DLL_API void ComplexSolveLU(
		__complex *X, int n, int LUIndexsLength, int *LUPtr, int *LUIndexs, __complex *LUValues, __complex *B);

    DLL_API void ComplexCalcILUWithPivoting(int* LUIndexsLengthP, int** LUPtrsP, int** LUIndexsP, __complex** LUValuesP,
        int* pivot, int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, int fillinLevel);

	DLL_API void ComplexCalcIC(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, __complex **LUValuesP,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues);

    ////////////////////////////////////////////////////
    // ComplexCOCG
    DLL_API bool ComplexSolvePreconditionedCOCG(__complex* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B,
        int LUIndexsLength, int* LUPtrs, int* LUIndexs, __complex* LUValues,
        double convRatioTolerance);

    DLL_API bool ComplexSolveCOCG(__complex* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B, int fillinLevel,
        double convRatioTolerance);

    DLL_API bool ComplexSolveICCOCG(__complex *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *B,
		double convRatioTolerance);

    ////////////////////////////////////////////////////
    // ComplexBiCGSTAB
    DLL_API bool ComplexSolveNoPreconBiCGSTAB(__complex *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *B,
		double convRatioTolerance);

    DLL_API bool ComplexSolvePreconditionedBiCGSTAB(__complex* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B,
        int LUIndexsLength, int* LUPtrs, int* LUIndexs, __complex* LUValues,
        double convRatioTolerance);

    DLL_API bool ComplexSolveBiCGSTAB(__complex* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B, int fillinLevel,
        double convRatioTolerance);

    DLL_API bool ComplexSolveBiCGSTABWithPivoting(__complex* X,
        int n, int AIndexsLength, int* APtrs, int* AIndexs, __complex* AValues, __complex* B, int fillinLevel,
        double convRatioTolerance);
}
