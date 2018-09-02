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

	DLL_API void DoubleDeleteCSR(int *APtrsP, int *AIndexsP, double *AValuesP);

	// Note: 内部でメモリ割り当てを行っているので使用後DeleteCSRを呼び出す必要がある
	DLL_API void DoubleCalcILU(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, int fillinLevel);

	DLL_API void DoubleSolveLU(
		double *X, int n, int LUIndexsLength, int *LUPtr, int *LUIndexs, double *LUValues, double *B);

	DLL_API bool DoubleSolvePreconditionedCG(double *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
		int LUIndexsLength, int *LUPtrs, int *LUIndexs, double *LUValues);

	DLL_API bool DoubleSolveCG(double *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B, int fillinLevel);

	//////////////////////////////////////////////////////////////////////////
	// ComplexLinear
	DLL_API void ComplexDotc(__complex *ret, int n, __complex *X, __complex *Y);

	DLL_API void ComplexDotu(__complex *ret, int n, __complex *X, __complex *Y);

	DLL_API void ComplexMV(__complex *Z,
		__complex alpha, int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *X,
		__complex beta, __complex *Y);

	DLL_API void ComplexAxpy(__complex *Z, __complex alpha, int n, __complex *X, __complex *Y);

	DLL_API void ComplexDeleteCSR(int *APtrsP, int *AIndexsP, __complex *AValuesP);

	// Note: 内部でメモリ割り当てを行っているので使用後DeleteCSRを呼び出す必要がある
	DLL_API void ComplexCalcILU(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, __complex **LUValuesP,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, int fillinLevel);

	DLL_API void ComplexSolveLU(
		__complex *X, int n, int LUIndexsLength, int *LUPtr, int *LUIndexs, __complex *LUValues, __complex *B);

	DLL_API bool ComplexSolvePreconditionedCOCG(__complex *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *B,
		int LUIndexsLength, int *LUPtrs, int *LUIndexs, __complex *LUValues);

	DLL_API bool ComplexSolveCOCG(__complex *X,
		int n, int AIndexsLength, int *APtrs, int *AIndexs, __complex *AValues, __complex *B, int fillinLevel);

}
