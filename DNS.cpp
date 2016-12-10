#include "DNS.h"
#include <mpi.h>
#include <iostream>
#define iDIM 0
#define jDIM 1
#define kDIM 2



void create_topology(MeshInfo *info) {

	
	MPI_Comm_size(MPI_COMM_WORLD, &info->numprocs);
	
	
	for( info->q = 1; (info->q+1)*(info->q+1)*(info->q+1) <= info->numprocs; info->q++ );

	// Создаем топологию, которую будем исспользовать
	int dims[3]; //число процессов вдоль каждого измерения(3 измерения)
	int periods[3];  // граничные условия
	dims[iDIM] = info->q;
	dims[jDIM] = info->q;
	dims[kDIM] = info->q;
	periods[iDIM] = periods[jDIM] = periods[kDIM] = 1;
	info->mesh3d = MPI_COMM_WORLD;


	int err = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &info->mesh3d);   //создания коммуникатора с декартовой топологией


	MPI_Cart_coords(info->mesh3d, info->myrank, 3, info->coords);   //определения координат процесса по его идентификатору 

	// i-k плоскости
	dims[iDIM] = dims[kDIM] = 1;
	dims[jDIM] = 0;

	MPI_Cart_sub(info->mesh3d, dims, &info->mesh_ik); //выделение подпространства в декартовой топологии

	// j кольцо
	dims[iDIM] = dims[kDIM] = 0;
	dims[jDIM] = 1;

	MPI_Cart_sub(info->mesh3d, dims, &info->ring_j);

	// j-k плоскости
	dims[jDIM] = dims[kDIM] = 1;
	dims[iDIM] = 0;

	MPI_Cart_sub(info->mesh3d, dims, &info->mesh_jk);

	//i кольцо
	dims[jDIM] = dims[kDIM] = 0;
	dims[iDIM] = 1;

	MPI_Cart_sub(info->mesh3d, dims, &info->ring_i);

	// i-j плоскости 
	dims[iDIM] = dims[jDIM] = 1;
	dims[kDIM] = 0;

	MPI_Cart_sub(info->mesh3d, dims, &info->mesh_ij);

	// k кольцо
	dims[iDIM] = dims[jDIM] = 0;
	dims[kDIM] = 1;

	MPI_Cart_sub(info->mesh3d, dims, &info->ring_k);
	//printf("%d %d %d %d", info->myrank, info->coords[0], info->coords[1], info->coords[2]);
}

Matrix* MultiplyMatrixesDNS(Matrix* A, Matrix* B, MeshInfo * t)
{

	
	Matrix * Ablock = DistributeLeftMatrix(A, t);
	Matrix * Bblock = DistributeRightMatrix(B, t);
	Matrix * Cblock = MultiplyMatrixes(Ablock, Bblock);
	Matrix * C = nullptr;

	int * buf = nullptr;
//	PrintMatrix(Bblock);
	int blksz = t->n / t->q;



	
		 buf = (int*)malloc(sizeof(int)*Cblock->n*Cblock->n);
		 memset(buf, 0, sizeof(int)*Cblock->n*Cblock->n);


	if(t->myrank)
	MPI_Reduce(MatrixToArrByCols(Cblock), buf, Cblock->n*Cblock->n, MPI_INT, MPI_SUM, 0, t->ring_j); //с сохранением результата в адресном пространстве одного процесса

	int *Cblocks = nullptr;
	if (t->coords[1] == 0)
	{


		Cblocks = (int*)malloc(sizeof(int)* A->n * B->n);

		MPI_Gather(buf, blksz*blksz, MPI_INT, Cblocks, blksz*blksz, MPI_INT, 0, t->mesh_ik); //сбор блоков данных от всех процессов группы


		if(buf != nullptr)free(buf);
	}


	if (t->myrank == 0) {
		C = ArrBlocksToMatrix(Cblocks, blksz,blksz, A->n, A->n);		
		free(Cblocks);


	}


	FreeMatrix(Ablock);
	FreeMatrix(Bblock);
	FreeMatrix(Cblock);

	return C;
}

Matrix *  DistributeLeftMatrix(Matrix* A, MeshInfo * topology)
{
	Matrix * block;
	int blockSize = topology->n / topology->q;
	int * buf = (int*)malloc(sizeof(int)*blockSize*blockSize);
	int * send = nullptr;
	if (topology->myrank == 0)
	{
		
		send = MatrixToArrBlocksRows(A, blockSize, blockSize);
	}
	if (topology->coords[2] == 0)
	{
		MPI_Scatter(send, blockSize*blockSize, MPI_INT, buf, blockSize*blockSize, MPI_INT, 0, topology->mesh_ij); // распределение блоков данных по всем процессам группы
	}

	if (topology->myrank == 0)
	{
		free(send);
	}

	MPI_Bcast(buf, blockSize*blockSize, MPI_INT, 0, topology->ring_k);
	block = ArrToMatrixByRows(buf, blockSize, blockSize);
	free(buf);
	//PrintMatrix(block);
	return block;
	
}

Matrix * DistributeRightMatrix(Matrix* A, MeshInfo * topology)
{
	Matrix * block;
	int blockSize = topology->n / topology->q;
	int * buf = (int*)malloc(sizeof(int)*blockSize*blockSize);
	int * send = nullptr;
	if (topology->myrank == 0)
	{
		
		send = MatrixToArrBlocksRows(A, blockSize, blockSize);
	}
	if (topology->coords[0] == 0)
	{
		MPI_Scatter(send, blockSize*blockSize, MPI_INT, buf, blockSize*blockSize, MPI_INT, 0, topology->mesh_jk);
	}
	if (topology->myrank == 0)
	{
		free(send);
	}
	MPI_Bcast(buf, blockSize*blockSize, MPI_INT, 0, topology->ring_i);
	block = ArrToMatrixByRows(buf, blockSize, blockSize);
	free(buf);
	//PrintMatrix(block);
	return block;
}

int* MatrixToArrByRows(Matrix* A)
{
	int * a = (int*)malloc(sizeof(int) * A->n*A->m);
	for (int i = 0; i < A->m; i++)
	{
		for (int j = 0; j < A->n; j++)
		{
			a[i*A->n + j] = A->data[i][j];
		}
	}
	return a;
}

int* MatrixToArrByCols(Matrix* A)
{

	int * a = (int*)malloc(sizeof(int) * A->n*A->m);
	for (int i = 0; i < A->n; i++)
	{
		for (int j = 0; j < A->m; j++)
		{
			a[i*A->m + j] = A->data[j][i];
		}
	}


	return a;
}

Matrix* ArrToMatrixByRows(int* a, int n, int m)
{
	Matrix * A = CreateMatrix(n, m);
	for (int i = 0; i < A->m; i++)
	{
		for (int j = 0; j < A->n; j++)
		{
			A->data[i][j] = a[i*A->n + j];
		}
	}
	return A;
}

Matrix* ArrToMatrixByCols(int* a, int n, int m)
{
	Matrix * A = CreateMatrix(n, m);
	for (int i = 0; i < A->n; i++)
	{
		for (int j = 0; j < A->m; j++)
		{
			 A->data[j][i] = a[i*A->m + j];
		}
	}
	return A;
}

int* MatrixToArrBlocksRows(Matrix* A, int blockrow, int blockcol)
{ 
	int fullrow = A->m / blockrow;
	int fullcol = A->n / blockcol;
	int remainsrow = A->m % blockrow;
	int remainscol = A->n % blockcol;
	int * a = (int *)malloc(sizeof(int) * A->n*A->m);
	int count = 0;
	int i, j;
	for ( i = 0; i < fullrow; i++)
	{
		for ( j = 0; j < fullcol; j++)
		{
			CopyBlockToArr(A, a, i*blockrow, j*blockcol, blockrow, blockcol, count);
			count += blockrow*blockcol;
		}
		
			CopyBlockToArr(A, a, i*blockrow, blockcol*fullcol, blockrow, remainscol, count);
			count += blockrow*remainscol;
		
		
	}	
		for (j = 0; j < fullcol; j++)
		{
			CopyBlockToArr(A, a, i*blockrow, j*blockcol, remainsrow, blockcol, count);
			count += remainsrow*blockcol;
		}
		CopyBlockToArr(A, a, fullcol*blockrow, blockcol*fullcol, remainsrow, remainscol, count);
	
	
	return a;
}

int* MatrixToArrBlocksCols(Matrix* A, int blockrow, int blockcol)
{
	int fullrow = A->m / blockrow;
	int fullcol = A->n / blockcol;
	int remainsrow = A->m % blockrow;
	int remainscol = A->n % blockcol;
	int * a = (int *)malloc(sizeof(int) * A->n*A->m);
	int count = 0;
	int i, j;
	for (j = 0; j < fullcol; j++)
	
	{
		for (i = 0; i < fullrow; i++)
		{
			CopyBlockToArr(A, a, i*blockrow, j*blockcol, blockrow, blockcol, count);
			count += blockrow*blockcol;
		}

		CopyBlockToArr(A, a, i*blockrow, j*blockcol, remainsrow, blockcol, count);
		count += remainsrow*blockcol;
		CopyBlockToArr(A, a, i*blockrow, blockcol*fullcol, blockrow, remainscol, count);
		count += blockrow*remainscol;


	}
	for (i = 0; i < fullrow; i++)
	{
		CopyBlockToArr(A, a, i*blockrow, blockcol*fullcol, blockrow, remainscol, count);
		count += blockrow*remainscol;
	}
	CopyBlockToArr(A, a, fullcol*blockrow, blockcol*fullcol, remainsrow, remainscol, count);


	return a;
}

void CopyBlockToArr(Matrix* A, int* a, int row, int col, int height, int width, int position)
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			a[position + i * width + j] = A->data[row + i][col + j];
		}
	}
}

void CopyArrBlocktoMatrix(Matrix* A, int* a, int row, int col, int height, int width, int position)
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			A->data[col + j][row + i] = a[position + i * width + j];
		}
	}
}

Matrix* ArrBlocksToMatrix(int* a, int blockrow, int blockcol, int n, int m)
{
	Matrix * A = CreateMatrix(n, m);
	int fullrow = A->m / blockrow;
	int fullcol = A->n / blockcol;
	int remainsrow = A->m % blockrow;
	int remainscol = A->n % blockcol;

	int count = 0;
	int i, j;
	for (j = 0; j < fullcol; j++)

	{
		for (i = 0; i < fullrow; i++)
		{
			CopyArrBlocktoMatrix(A, a, i*blockrow, j*blockcol, blockrow, blockcol, count);
			count += blockrow*blockcol;
		}

		CopyArrBlocktoMatrix(A, a, i*blockrow, j*blockcol, remainsrow, blockcol, count);
		count += remainsrow*blockcol;
		CopyArrBlocktoMatrix(A, a, i*blockrow, blockcol*fullcol, blockrow, remainscol, count);
		count += blockrow*remainscol;


	}
	for (i = 0; i < fullrow; i++)
	{
		CopyArrBlocktoMatrix(A, a, i*blockrow, blockcol*fullcol, blockrow, remainscol, count);
		count += blockrow*remainscol;
	}
	CopyArrBlocktoMatrix(A, a, fullcol*blockrow, blockcol*fullcol, remainsrow, remainscol, count);


	return A;
}



