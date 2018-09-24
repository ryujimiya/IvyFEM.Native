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

bool DoubleSolveNoPreconCG(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
	double convRatioTolerance)
{
	double convRatio = convRatioTolerance;
	double tolerance = convRatio;
	int maxIter = MaxIter;
	double *r = new double[n];
	double *p = new double[n];
	double *Ap = new double[n];
	double *dummy = new double[n];
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
			printf("iter = %d norm: %e\r\n", iter, convRatio);

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
				printf("iter = %d norm: %e\r\n", iter, convRatio);

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
		printf("iter = %d norm: %e\r\n", iter, convRatio);
	}
	printf("Not converged\r\n");

	__FINALIZE
#undef __FINALIZE
	return false;
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
					printf("__DoubleCalcILU err divide by zero\r\n");
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
					// LU[i, j] : キーがなければ生成される
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
	int LUIndexsLength, int *LUPtrs, int *LUIndexs, double *LUValues,
	double convRatioTolerance)
{
	double convRatio = convRatioTolerance;
	double tolerance = convRatio;
	int maxIter = MaxIter;
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
#undef __FINALIZE
	return false;
}

bool DoubleSolveCG(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B, int fillinLevel,
	double convRatioTolerance)
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
		LUIndexsLength, LUPtrs, LUIndexs, LUValues,
		convRatioTolerance);
	DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
	printf("    2: t = %d\r\n", GetTickCount() - t);
	fflush(stdout);
	return success;
}

/*
template<typename KeyT, typename ValueT>
static void StdMapCopy(std::map<KeyT, ValueT> &dest, const std::map<KeyT, ValueT> &src)
{
	dest.clear();
	for (std::pair<KeyT, ValueT> pair : src)
	{
		dest.insert(std::make_pair(pair.first, pair.second));
	}
}
*/

static void __DoubleCalcILUWithPivoting(
	int n, std::map<int, double> *LUIndexValues, int *pivot, int fillinLevel)
{
	for (int row = 0; row < n; row++)
	{
		pivot[row] = row;
	}
	int *level = new int[n * n];
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			level[row + col * n] = LUIndexValues[row].find(col) != LUIndexValues[row].end() ? // LU[row, col]
				0 : (fillinLevel + 1);
		}
	}

	for (int k = 0; k < (n - 1); k++)
	{
		int p = k;
		{
			double max = 0;
			{
				double LUkk = 0;
				auto LUkkItr = LUIndexValues[k].find(k); // LU[k, k]
				if (LUkkItr != LUIndexValues[k].end())
				{
					LUkk = (*LUkkItr).second;
				}
				max = fabs(LUkk);
			}
			for (int i = k + 1; i < n; i++)
			{
				double LUik;
				{
					auto LUikItr = LUIndexValues[i].find(k);
					if (LUikItr == LUIndexValues[i].end())
					{
						continue;
					}
					LUik = (*LUikItr).second;
				}
				double abs = fabs(LUik); // LU[i, k]
				if (abs > max)
				{
					max = abs;
					p = i;
				}
			}
		}
		if (k != p)
		{
			{
				int tmp = pivot[k];
				pivot[k] = pivot[p];
				pivot[p] = tmp;
			}
			{
				std::map<int, double> tmp = LUIndexValues[k];
				LUIndexValues[k] = LUIndexValues[p];
				LUIndexValues[p] = tmp;
			}
			for (int j = 0; j < n; j++)
			{
				int tmp = level[k + j * n];
				level[k + j * n] = level[p + j * n];
				level[p + j * n] = tmp;
			}
		}
		for (int i = k + 1; i < n; i++)
		{
			if (level[i + k * n] > fillinLevel)
			{
				continue;
			}
			double LUik;
			{
				auto LUkkItr = LUIndexValues[k].find(k);
				if (LUkkItr == LUIndexValues[k].end())
				{
					printf("__DoubleCalcILUWithPivoting divide by zero\r\n");
					fflush(stdout);
					throw;
				}
				auto LUikItr = LUIndexValues[i].find(k);
				if (LUikItr == LUIndexValues[i].end())
				{
					continue;
				}
				(*LUikItr).second /= (*LUkkItr).second; // LU[i, k] LU[k, k]
				LUik = (*LUikItr).second;
			}
			for(auto pair : LUIndexValues[k])
			{
				int j = pair.first;
				double LUkj = pair.second;
				if (j >= k + 1 && j < n)
				{

				}
				else
				{
					continue;
				}
				level[i + j * n] = min(level[i + j * n], level[i + k * n] + level[k + j * n] + 1);
				if (level[i + j * n] <= fillinLevel)
				{
					// LU[i, j] : キーがなければ生成される
					LUIndexValues[i][j] -= LUik * LUkj; // LU[i, k] LU[k, j]
				}
			}
		}
	}

	delete[] level;
}

void DoubleCalcILUWithPivoting(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
	int *pivot, int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, int fillinLevel)
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

	__DoubleCalcILUWithPivoting(n, LUIndexValues, pivot, fillinLevel);

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

bool DoubleSolveCGWithPivoting(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B, int fillinLevel,
	double convRatioTolerance)
{
	int t;
	int LUIndexsLength = 0;
	int *LUPtrs = NULL;
	int *LUIndexs = NULL;
	double *LUValues = NULL;
	int *pivot = new int[n];
	int pivotingAIndexsLength = 0;
	int *pivotingAPtrs = new int[n + 1];
	int *pivotingAIndexs = new int[AIndexsLength];
	double *pivotingAValues = new double[AIndexsLength];
	double *pivotingB = new double[n];

	t = GetTickCount();
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
	printf("    1: t = %d\r\n", GetTickCount() - t);
	t = GetTickCount();
	bool success = DoubleSolvePreconditionedCG(
		X,
		n, pivotingAIndexsLength, pivotingAPtrs, pivotingAIndexs, pivotingAValues, pivotingB,
		LUIndexsLength, LUPtrs, LUIndexs, LUValues,
		convRatioTolerance);
	DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
	printf("    2: t = %d\r\n", GetTickCount() - t);
	fflush(stdout);

	delete[] pivot;
	delete[] pivotingAPtrs;
	delete[] pivotingAIndexs;
	delete[] pivotingAValues;
	delete[] pivotingB;
	return success;
}

static void __DoubleCalcIC(
	int n, std::map<int, double> *LUIndexValues, std::map<int, double> *AIndexValues)
{
	double *D = new double[n];
	memset(D, 0, n * sizeof(double));

	D[0] = AIndexValues[0][0]; // A[0, 0]
	LUIndexValues[0][0] = 1.0; // LU[0, 0]

	// L
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j <= (i - 1); j++)
		{
			double Aij;
			{
				auto AijItr = AIndexValues[i].find(j);
				if (AijItr == AIndexValues[i].end())
				{
					continue;
				}
				Aij = (*AijItr).second;
			}
			double tmp = Aij; // A[i, j]
			for (auto pair : LUIndexValues[i])
			{
				int k = pair.first;
				double LUik = pair.second;
				if (k >= 0 && k <= (j - 1))
				{
					//
				}
				else
				{
					continue;
				}
				double LUjk;
				{
					auto LUjkItr = LUIndexValues[j].find(k);
					if (LUjkItr == LUIndexValues[j].end())
					{
						continue;
					}
					LUjk = (*LUjkItr).second;
				}
				tmp -= LUik * LUjk * D[k];  // LU[i, k] LU[j, k]
			}
			LUIndexValues[i][j] = (1.0 / D[j]) * tmp; // LU[i, j]
		}

		// i == jのとき
		{
			double tmp = AIndexValues[i][i]; //A[i, i]
			for (auto pair : LUIndexValues[i])
			{
				int k = pair.first;
				double LUik = pair.second;
				if (k >= 0 && k <= (i - 1))
				{
					//
				}
				else
				{
					continue;
				}
				tmp -= LUik * LUik * D[k];  // LU[i, k]
			}
			D[i] = tmp;
			LUIndexValues[i][i] = 1.0; // LU[i, i] 
		}
	}

	// DL^T
	for (int i = 0; i < n; i++)
	{
		for (auto pair : LUIndexValues[i])
		{
			int j = pair.first;
			double LUij = pair.second;
			if (j >= 0 && j <= (i - 1))
			{
				//
			}
			else
			{
				continue;
			}
			LUIndexValues[j][i] = D[j] * LUij; // LU[j, i] LU[i, j]
		}
		LUIndexValues[i][i] = D[i]; // LU[i, i]
	}

	delete[] D;
}

void DoubleCalcIC(int *LUIndexsLengthP, int **LUPtrsP, int **LUIndexsP, double **LUValuesP,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues)
{
	std::map<int, double> *LUIndexValues = new std::map<int, double>[n];
	std::map<int, double> *AIndexValues = new std::map<int, double>[n];
	for (int row = 0; row < n; row++)
	{
		int sPtr = APtrs[row];
		int ePtr = APtrs[row + 1]; // Note: APtr size is n + 1, APtr[n + 1] is always AIndexsLength
		for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
		{
			int col = AIndexs[iPtr];
			double value = AValues[iPtr];
			AIndexValues[row].insert(std::make_pair(col, value)); // A[row, col]
		}
	}

	__DoubleCalcIC(n, LUIndexValues, AIndexValues);

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
	delete[] AIndexValues;
}

bool DoubleSolveICCG(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
	double convRatioTolerance)
{
	int t;
	int LUIndexsLength = 0;
	int *LUPtrs = NULL;
	int *LUIndexs = NULL;
	double *LUValues = NULL;
	t = GetTickCount();
	DoubleCalcIC(&LUIndexsLength, &LUPtrs, &LUIndexs, &LUValues,
		n, AIndexsLength, APtrs, AIndexs, AValues);
	printf("    1: t = %d\r\n", GetTickCount() - t);
	t = GetTickCount();
	bool success = DoubleSolvePreconditionedCG(
		X,
		n, AIndexsLength, APtrs, AIndexs, AValues, B,
		LUIndexsLength, LUPtrs, LUIndexs, LUValues,
		convRatioTolerance);
	DoubleDeleteCSR(LUPtrs, LUIndexs, LUValues);
	printf("    2: t = %d\r\n", GetTickCount() - t);
	fflush(stdout);
	return success;
}

bool DoubleSolveNoPreconBiCGSTAB(double *X,
	int n, int AIndexsLength, int *APtrs, int *AIndexs, double *AValues, double *B,
	double convRatioTolerance)
{
	double convRatio = convRatioTolerance;
	double tolerance = convRatio;
	int maxIter = MaxIter;
	double *r = new double[n];
	double *r0 = new double[n];
	double *p = new double[n];
	double *s = new double[n];
	double *Ap = new double[n];
	double *As = new double[n];
	double *dummy = new double[n];
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
			printf("iter = %d norm: %e\r\n", iter, convRatio);

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
		DoubleAxpy(p, -beta * omega, n, Ap, p);
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
