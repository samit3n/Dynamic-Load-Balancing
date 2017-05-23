/**
*
* MPI implementation of scalar product of two vectors
*/

#include <mpi.h>
#include <iostream>
#include <random>

using namespace std;


int main(int argc, char ** argv)
{
	
	int rank, size;

	unsigned vecSize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//some global variables amongst processes

	int * vec1, * vec2; 
	int * cnts, * displs; 		//portion sizes and memory displacements

	int *rbuf1, *rbuf2; 		//receive buffers

	/*
		Parameter check
	*/


	if(rank == 0){				

		if(argc < 2){
			cout << "Vector size defautls to " << size << endl;
			vecSize = size;
		}else{

			vecSize = stoi(argv[1]);

			if(vecSize <= 0 ){
				cerr << "Invalid vector len " << vecSize << endl;
				return 2;
			}
		}
	
		//vectors initialization

		vec1 = new int[vecSize];
		vec2 = new int[vecSize];
 

		//fill vectors by ascending integer sequence 

		for(int i = 0; i < vecSize; i++){

			vec1[i] = i;
			vec2[i] = i;

		}

		cnts = new int [size];
		displs = new int [size];

		int sum = 0;

		//simple portion assignment
		//last cpu takes the rest

		int perProc = vecSize / size;

		if(vecSize % size != 0){
			cout << "perCPU " << perProc <<  " modulo " << vecSize % size <<  endl;


			for(int i = 0; i < size-1; i++){
			
				cnts[i]=perProc;
				displs[i] = sum;
				sum += 	perProc;
			}
			cnts[size-1] = perProc + (vecSize % size);
			cout << "last " << cnts[size-1] << endl;
			displs[size-1] = sum;
			//perProc + (vecSize % size);

		}else{

			for(int i = 0; i < size; i++){
			
				cnts[i]= perProc;
				displs[i] = sum;
				sum += perProc;
			}

		}

		

	}
		
	int portion;

	//scatter assigned portion to processes

	MPI_Scatter(cnts, 1, MPI_INT, &portion, 1,
				 MPI_INT, 0, MPI_COMM_WORLD );

	rbuf1 = new int[portion];
	rbuf2 = new int[portion];

	//scatter vector data

	MPI_Scatterv(vec1, cnts, displs, MPI_INT, rbuf1 ,
					 portion, MPI_INT, 0, MPI_COMM_WORLD
					);

	MPI_Scatterv(vec2, cnts, displs, MPI_INT, rbuf2,
					 portion, MPI_INT, 0, MPI_COMM_WORLD
					);

	cout << "Proc " << rank << " received " << portion << " elems ";

	for(int i = 0 ; i < portion; i++){
		cout << '(' << rbuf1[i] << '*'<< rbuf2[i] << ") + ";
	}
	cout << endl;

	//compute

	int sum = 0;
	for(int i = 0; i < portion; i++){
		sum += rbuf1[i] * rbuf2[i];
	}

	int prod;

	//reduce SUM to root
	
	MPI_Reduce(&sum, &prod, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

	if(rank == 0){
		cout << "Scalar vector product is " << prod << endl;
	}


	delete[] rbuf1, rbuf2;
	delete[] cnts, displs;

	if(rank == 0) delete[] vec1, vec2;

	MPI_Finalize();



}


 