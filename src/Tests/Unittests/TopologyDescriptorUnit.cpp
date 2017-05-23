/***********************************************
*
* 	File Name:		TileDescriptorTestBench.cpp

*	Project: 		DIP - Dynamic Load Balancing in HPC Applications
*	
*	Description:	Simple testbench for tile descriptor.
*
*	Author: 		
* 	Email:			
* 	Date:  			
*
***********************************************/

#include "../Sources/TopologyDescriptor.h"


int main()
{

	int rank, size;
	
	// Initialize MPI
	MPI_Init(&argc, &argv);
	// Get MPI rank and size
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	TopologyDescriptor td;



	MPI_Finalize();

}
