#pragma once
#include "Matrix.h"
#include <mpi.h>

typedef struct 
{

	MPI_Comm mesh3d;
	MPI_Comm mesh_ik, ring_j;
	MPI_Comm mesh_jk, ring_i;
	MPI_Comm mesh_ij, ring_k;
	int myrank =0;
	int coords[3];
	int numprocs; 	
	int q;	
	int n; 
	
} MeshInfo;
Matrix *  MultiplyMatrixesDNS(Matrix * A, Matrix * B, MeshInfo * topology);
Matrix * DistributeLeftMatrix(Matrix * A, MeshInfo * topology);
Matrix * DistributeRightMatrix(Matrix * A, MeshInfo  * topology);
int* MatrixToArrByRows(Matrix* A);
int * MatrixToArrByCols(Matrix * B);
Matrix * ArrToMatrixByRows(int * a, int n, int m);
Matrix * ArrToMatrixByCols(int * a, int n, int m);
int * MatrixToArrBlocksRows(Matrix * A, int blockrow, int blockcol);
int * MatrixToArrBlocksCols(Matrix * A, int blockrow, int blockcol);
void CopyBlockToArr(Matrix * A, int * a, int row, int col, int height, int width, int position);
void CopyArrBlocktoMatrix(Matrix * A, int * a, int row, int col, int height, int width, int position);
Matrix * ArrBlocksToMatrix(int * a, int blockrow, int blockcol, int n, int m);
void create_topology(MeshInfo *info);
