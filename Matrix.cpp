#include "Matrix.h"
#include <corecrt_malloc.h>
#include<iostream>
using namespace std;

Matrix* CreateMatrix(int n, int m)
{
	
		Matrix * matr;
		matr = (Matrix *)malloc(sizeof(Matrix));
		matr->n = n;
		matr->m = m;
		MallocMatrix(matr);
		return matr;


}

void MallocMatrix(Matrix* matr)
{
	
	matr->data = (int **)malloc(matr->m * sizeof(int *));	
	for (int i = 0; i < matr->m; i++)
		matr->data[i] = (int *)malloc(matr->n * sizeof(int));	
}

void FreeMatrix(Matrix* A)
{
	for (int i = 0; i < A->m; i++)
		free(A->data[i]);
	free(A->data);
	A->n = -1;
	A->m = -1;
	free(A);
}


void FillMatrixZero(Matrix* matr)
{
	for (int i = 0; i < matr->m; i++)
	{
		for (int j = 0; j < matr->n; j++)
		{
			matr->data[i][j] = 0;
		}
	}
}

void FillMatrixNumbers(Matrix* matr)
{
	for (int i = 0; i < matr->m; i++)
	{
		for (int j = 0; j < matr->n; j++)
		{
			matr->data[i][j] = i*matr->n+j;
		}
	}
}

void FillMatrixDiagonal(Matrix* matr)
{
	for (int i = 0; i < matr->m; i++)
	{
		for (int j = 0; j < matr->n; j++)
		{
			if (i == j)
				matr->data[i][j] = 0;
			else
				matr->data[i][j] = 1;
		}
	}
}

void PrintMatrix(Matrix* matr)
{
	for (int i = 0; i < matr->m; i++)
	{
		for (int j = 0; j < matr->n; j++)
		{
			cout << matr->data[i][j] << " ";
		}
		cout << endl;
	}
}

Matrix* AddMatrixes(Matrix* A, Matrix* B)
{
	if (A->n != B->n || A->m != B->m)
		return nullptr;
	Matrix * C = CreateMatrix(A->n, A->m);
	for (int i = 0; i < C->m; i++)
	{
		for (int j = 0; j < C->n; j++)
		{
			C->data[i][j] = A->data[i][j] + B->data[i][j];
		}
	}
	return C;
	
}

Matrix* MultiplyMatrixes(Matrix* A, Matrix* B)
{
	if (A->n != B->m)
		return nullptr;
	Matrix * C = CreateMatrix(A->m, B->n);
	FillMatrixZero(C);
	for (int i = 0; i < C->m; i++)
	{
		for (int j = 0; j < C->n; j++)
		{
			for (int k = 0; k < A->n; k++)
			{
				C->data[i][j] += (A->data[i][k] * B->data[k][j]);
			}			
		}
	}
	return C;

}