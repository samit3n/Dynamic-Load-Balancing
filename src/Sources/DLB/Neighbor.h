/***********************************************
*
*  File Name:       HaloBuffers.cpp
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Neighbor describing class
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/

#ifndef __NEIGHBOR_H__
#define __NEIGHBOR_H__

#include <mpi.h>

namespace DLB {


/**
 * @brief Creates neighbor with all necessary attributes
 * 
 * @param rank of neighbor
 * @param sendTag - tag to be passed to send call
 * @param recvTag - tag to be passed to recv call
 * @param pos - informative - position from actual rank
 * @param send - halo zone to be send
 * @param recv - halo zone to be received to
 * 
*/


class Neighbor {

public:


	Neighbor() 
	{
		wRank = -1;
		root = -1;
		myRank = -1;

		comm = MPI_COMM_NULL;
		count = 0;
		displ= 0;

		scatterCnts = NULL;
		scatterDispls = NULL;
	}

	~Neighbor()
	{
		// if(scatterCnts != NULL){
			// delete[] scatterCnts;
		// }
// 
		// if(scatterDispls != NULL){
			// delete[] scatterDispls;
		// }

	}	

	void clear(void)
	{
		wRank = -1;
		root = -1;
		myRank = -1;

		comm = MPI_COMM_NULL;
		count = 0;
		displ= 0;

		if(scatterCnts != NULL){
			delete[] scatterCnts;
			scatterCnts = NULL;
		}

		if(scatterDispls != NULL){
			delete[] scatterDispls;
			scatterDispls = NULL;
		}

	}


	int wRank; 		// his rank in COMM_WORLD
	int root;			// his rank HIS communicator (eg. root)

	MPI_Comm comm;  	// his communicator

	int myRank; 		//myRank in his communicator

	unsigned count; 	// count of items which I share with him
	unsigned displ;	// displacement in my buffer

	// valus passed directly to Iscatterv
	// aray of comm_size(comm)
	// must be freed manually
	// is Bcaster from neighbor, he passes the same

	int * scatterCnts;	
	int * scatterDispls;




};


} //DLB nspace end

#endif
