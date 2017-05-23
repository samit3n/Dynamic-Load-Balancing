/*
 * =====================================================================================
 *
 *       Filename:  wall_time.c
 *
 *    Description:  Wall time measure test 
 *
 *        Version:  1.0
 *        Created:  11/04/2016 01:40:10 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include <mpi.h>
#include <stdio.h>

using namespace std;

int main(int argc, char ** argv)
{

    double el_time;

    el_time = MPI_Wtime();

    cout << "Start timestamp " << el_time << endl;

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

    cout << "Stop timestamp " << MPI_Wtime() << " - elapsed " << MPI_Wtime() - el_time << endl;
	MPI_Finalize();

	return 0;
}
