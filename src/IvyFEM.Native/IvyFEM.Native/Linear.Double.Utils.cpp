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
