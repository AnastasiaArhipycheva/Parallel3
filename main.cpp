#include "mpi.h"
#include "Matrix.h"
#include "DNS.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>


int main( int argc, char *argv[] ) 
{

	MPI_Init(&argc, &argv);
	int rank;
	int n = 20;
	Matrix * A = nullptr;
	Matrix * B = nullptr;
	Matrix * C = nullptr;
	MeshInfo t;
	create_topology(&t);


	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double start;
	double end;	
	
		t.n = n;
	//	if (rank == 0)
	//	{

			A = CreateMatrix(n, n);
			FillMatrixNumbers(A);
			//PrintMatrix(A);
			B = CreateMatrix(n, n);
			FillMatrixNumbers(B);
			//PrintMatrix(B);
			start = MPI_Wtime();
	//	}

		std::cout << "I start multiply, heeey!" << std::endl;
		C = MultiplyMatrixesDNS(A, B, &t);

		if (rank == 0)
		{
			end = MPI_Wtime() - start;
		//	PrintMatrix(C);
			std::cout << "I end multy!" << n << ", time :" << end << std::endl;

		}
		if (rank == 0)
		{
			FreeMatrix(A);
			FreeMatrix(B);
			FreeMatrix(C);
		}

		
		MPI_Barrier(MPI_COMM_WORLD);
	
		
	//	std::cin;

	MPI_Finalize();
	
	return 0;

}

