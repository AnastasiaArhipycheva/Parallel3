#pragma once
typedef struct
{
	int ** data;
	int n; // ширина
	int m; // высота
} Matrix;
Matrix * CreateMatrix(int n, int m);
void MallocMatrix(Matrix* matr);
void FreeMatrix(Matrix * A);
void FillMatrixZero(Matrix * matr);
void FillMatrixNumbers(Matrix * matr);
void FillMatrixDiagonal(Matrix * matr);
void PrintMatrix(Matrix * matr);
Matrix * AddMatrixes(Matrix * A, Matrix * B);
Matrix * MultiplyMatrixes(Matrix * A, Matrix * B);


