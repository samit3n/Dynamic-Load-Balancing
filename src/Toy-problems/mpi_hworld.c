/** 
 * Simple MPI Hello world
 * 
 * Determines process rank and world size.
 *
 */

#include <mpi.h>
#include <stdio.h>

int main(int argc, char ** argv)
{

	MPI_Init(NULL, NULL);

	int worldSize;

	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	char procName[MPI_MAX_PROCESSOR_NAME];

	int nameLen;
	MPI_Get_processor_name(procName, &nameLen);

	printf("Hello, this is Salomon's HPC World. ");
	printf("I'm processor %s with ID %d out of %d\n", procName, myRank, worldSize);
	
	MPI_Finalize();

	return 0;
}


