#include "stdafx.h"
#include <vector>
#include <map>
#include <assert.h>
#include "IvyFEM.Native.h"
#include "Constants.h"

/*
template<typename T>
inline static long long StdFind(std::vector<T> v, T value)
{
	auto findItr = std::find(v.begin(), v.end(), value);
	if (findItr == v.end())
	{
		return -1;
	}
	long long index = findItr - v.begin();
	return index;
}
*/

double DoubleDot(int n, double *X, double *Y)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += X[i] * Y[i];
	}
	return sum;
}

void DoubleMV(double *Z,
	double alpha, int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *X,
	double beta, double *Y)
{
	for (int row = 0; row < n; row++)
	{
		double tmp = 0;
		int sPtr = APtrs[row];
		int ePtr = APtrs[row + 1]; // Note: APtr size is n + 1, APtr[n + 1] is always AIndexsLength
		for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
		{
			int col = AIndexs[iPtr];
			double value = AValues[iPtr];
			tmp += alpha * value * X[col];
		}
		Z[row] = tmp + beta * Y[row];
	}
}

void DoubleAxpy(double *Z, double alpha, int n, double *X, double *Y)
{
	for (int i = 0; i < n; i++)
	{
		Z[i] = alpha * X[i] + Y[i];
	}
}

static void DoubleGetCSR(int n, int *AIndexsLengthP, int **APtrsP, int **AIndexsP, double **AValuesP,
	std::map<int, double> *AIndexValues)
{
	int AIndexsLength = 0;
	int *APtrs = new int[n + 1];
	int *AIndexs = NULL;
	double *AValues = NULL;

	std::vector<int> tmpIndexs;
	std::vector<double> tmpValues;
	for (int row = 0; row < n; row++)
	{
		APtrs[row] = (int)tmpIndexs.size();
		for (auto pair : AIndexValues[row])
		{
			int col = pair.first;
			double value = pair.second;
			tmpIndexs.push_back(col);
			tmpValues.push_back(value);
		}
	}
	// note: LUPtr size is n + 1
	APtrs[n] = (int)tmpIndexs.size();

	AIndexsLength = (int)tmpIndexs.size();
	AIndexs = new int[AIndexsLength];
	AValues = new double[AIndexsLength];
	{
		for (int index = 0; index < AIndexsLength; index++)
		{
			AIndexs[index] = tmpIndexs[index];
			AValues[index] = tmpValues[index];
		}
	}

	*AIndexsLengthP = AIndexsLength;
	*APtrsP = APtrs;
	*AIndexsP = AIndexs;
	*AValuesP = AValues;
}

void DoubleDeleteCSR(int *APtrs, int *AIndexs, double *AValues)
{
	delete[] APtrs;
	delete[] AIndexs;
	delete[] AValues;
}

static void __DoubleCalcILU(
	int n, std::map<int, double> *LUIndexValues, int fillinLevel)
{
	int *level = new int[n * n];
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			level[row + col * n] = LUIndexValues[row].find(col) != LUIndexValues[row].end() ? // LU[row, col]
				0 : (fillinLevel + 1);
		}
	}

	for (int i = 1; i < n; i++)
	{
		for (int k = 0; k <= (i - 1); k++)
		{
			double LUik;
			{
				if (level[i + k * n] > fillinLevel)
				{
					continue;
				}
				auto LUikItr = LUIndexValues[i].find(k);
				if (LUikItr == LUIndexValues[i].end()) // LU[i, k]
				{
					continue;
				}
				auto LUkkItr = LUIndexValues[k].find(k);
				if (LUkkItr == LUIndexValues[k].end()) // LU[k, k]
				{
					printf("__DoubleCalcILU err 3\r\n");
					fflush(stdout);
					throw;
				}
				(*LUikItr).second /= (*LUkkItr).second;  //LU[i, k] LU[k, k]
				LUik = (*LUikItr).second;
			}
			for (auto pair : LUIndexValues[k])
			{
				int j = pair.first;
				double LUkj = pair.second;
				if (j >= k + 1 && j < n)
				{
					//
				}
				else
				{
					continue;
				}

				level[i + j * n] = min(level[i + j * n], level[i + k * n] + level[k + j * n] + 1);
				if (level[i + j * n] <= fillinLevel)
				{
					// LU[i, j] : ƒL[‚ª‚È‚¯‚ê‚Î¶¬‚³‚ê‚é
					LUIndexValues[i][j] -= LUik * LUkj; // LU[i, k] LU[k, j]
				}
			}
		}
	}

	delete[] level;
}

void DoubleCalcILU(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, int fillinLevel)
{
	std::map<int, double> *LUIndexValues = new std::map<int, double>[n];
	for (int row = 0; row < n; row++)
	{
		int sPtr = APtrs[row];
		int ePtr = APtrs[row + 1]; // Note: APtr size is n + 1, APtr[n + 1] is always AIndexsLength
		for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
		{
			int col = AIndexs[iPtr];
			double value = AValues[iPtr];
			LUIndexValues[row].insert(std::make_pair(col, value)); // LU[row, col]
		}
	}

	__DoubleCalcILU(n, LUIndexValues, fillinLevel);

	int LUIndexsLength = 0;
	int *LUPtrs = NULL;
	int *LUIndexs = NULL;
	double *LUValues = NULL;
	DoubleGetCSR(n, &LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues, LUIndexValues);

	*LUIndexsLengthP = LUIndexsLength;
	*LUPtrsP = LUPtrs;
	*LUIndexsP = LUIndexs;
	*LUValuesP = LUValues;

	delete[] LUIndexValues;
}

void DoubleSolveLU(
	double *X, int n, int LUIndexsLength, int *LUPtrs, int *LUIndexs, double *LUValues, double *B)
{
	memcpy(X, B, n * sizeof(double));

	// Ly = b 
	for (int row = 1; row < n; row++)
	{
		int sPtr = LUPtrs[row];
		int ePtr = LUPtrs[row + 1]; // Note: LUPtr size is n + 1, LUPtr[n + 1] is always LUIndexsLength
		for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
		{
			int col = LUIndexs[iPtr];
			double value = LUValues[iPtr];
			if (col >= 0 && col < row)
			{
				X[row] -= value * X[col]; // LU[row, col]
			}
		}
	}

	// Ux = y
	for (int row = n - 2; row >= 0; row--)
	{
		double lurr = 0;
		int sPtr = LUPtrs[row];
		int ePtr = LUPtrs[row + 1]; // Note: LUPtr size is n + 1, LUPtr[n + 1] is always AIndexsLength
		for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
		{
			int col = LUIndexs[iPtr];
			double value = LUValues[iPtr];
			if (col >= row + 1 && col < n)
			{
				X[row] -= value * X[col]; // LU[row, col]

			}
			else if (row == col)
			{
				lurr = value;
			}
		}
		X[row] /= lurr; // LU[row, row];
	}
}

bool DoubleSolvePreconditionedCG(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
	int LUIndexsLength, int *LUPtrs, int *LUIndexs, double *LUValues)
{
	double convRatio = ConvRatioTolerance;
	double tolerance = convRatio;
	int maxCnt = 1000;
	double *r = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *Ap = new double[n];
	double *dummy = new double[n];
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
			printf("iter = %d norm: %e\r\n", iter, convRatio);

			__FINALIZE
			return true;
		}
		sqInvNormRes0 = 1.0 / sqNormRes0;
	}

	DoubleSolveLU(z, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, r);

	memcpy(p, z, n * sizeof(double));
	double rz = DoubleDot(n, r, z);

	for (iter = 0; iter < maxCnt; iter++)
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
				printf("iter = %d norm: %e\r\n", iter, convRatio);

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
		printf("iter = %d norm: %e\r\n", iter, convRatio);
	}
	printf("Not converged\r\n");

	__FINALIZE
	return false;
}

bool DoubleSolveCG(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B, int fillinLevel)
{
	int t;
	int LUIndexsLength = 0;
	int *LUPtrs = NULL;
	int *LUIndexs = NULL;
	double *LUValues = NULL;
	t = GetTickCount();
	DoubleCalcILU(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues,
		n, AIndexsLength, APtrs, AIndexs, AValues, fillinLevel);
	printf("    1: t = %d\r\n", GetTickCount() - t);
	t = GetTickCount();
	bool success = DoubleSolvePreconditionedCG(
		X,
		n, AIndexsLength, APtrs, AIndexs, AValues, B,
		LUIndexsLength, LUPtrs, LUIndexs, LUValues);
	DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
	printf("    2: t = %d\r\n", GetTickCount() - t);
	fflush(stdout);
	return success;
}