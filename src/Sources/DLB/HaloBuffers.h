/***********************************************
*
*  File Name:       HaloBuffers.cpp
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Halo array container
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/

#ifndef __DLB_HALLO_BUFFERS_H
#define __DLB_HALLO_BUFFERS_H

#include <cstddef>

namespace DLB{


class HaloBuffers {

public:

	HaloBuffers(void):
	sendTemp(NULL),
	sendParams(NULL),
	sendMap(NULL),
	recvTemp(NULL),
	recvParams(NULL),
	recvMap(NULL)
	{}

	HaloBuffers(unsigned size)
	{
		sendTemp = new float[size];
		sendParams = new float[size];
		sendMap = new int[size];

		recvTemp = new float[size];
		recvParams = new float[size];
		recvMap = new int[size];
	}

	~HaloBuffers(void)
	{
		delete[] sendTemp;
		delete[] sendParams;
		delete[] sendMap;

		delete[] recvTemp;
		delete[] recvParams;
		delete[] recvMap;
	}

	void resize(unsigned newSize)
	{

		delete[] sendTemp;
		delete[] sendParams;
		delete[] sendMap;

		delete[] recvTemp;
		delete[] recvParams;
		delete[] recvMap;

		sendTemp = new float[newSize];
		sendParams = new float[newSize];
		sendMap = new int[newSize];

		recvTemp = new float[newSize];
		recvParams = new float[newSize];
		recvMap = new int[newSize];
		
	}



	//send
	float * sendTemp;
	float * sendParams;
	int * sendMap;

	//receive
	float * recvTemp;
	float * recvParams;
	int * recvMap;







};


}


#endif
