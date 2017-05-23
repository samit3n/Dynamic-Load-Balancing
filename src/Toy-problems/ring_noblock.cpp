/**
* MPI rign communication
* 
* Ring communication is defined as set of processes
* in ring topology, where every process want's to send
* message to both his neighbors.
* 
* This variant implements non-blocking version of MPI communication, 
* so no deadlock should occur.
* 
* Author: Vojtech Dvoracek
* email: xdvora0y@stud.fit.vutbr.cz
*/

#include <mpi.h>
#include <iostream>

using namespace std;

const int MYTAG = 123; 	//generic message tag

int main(int argc, char ** argv)
{

	if(argc < 2){
		cerr << "Give me msg len";
		return 0;
	}
	
	int rank, left, right, size;		//rank holding variables

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// computer neighbors based on my rank

	if(rank == 0){
		left = size-1;
		right = 1;
	}else if(rank == size-1){
		left = rank -1;
		right = 0;
	}else{
		
		left = rank -1;
		right = rank +1;
	}

	int msgsize = std::stoi(argv[1]);		//parsing argument, message size

	if(msgsize < 1){						// 0-byte message makes no sense
		cerr << "Message too small";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	int *Srmsg = new int[msgsize];			//message holding buffers
	int *Slmsg = new int[msgsize];
	int *Rrmsg = new int[msgsize];
	int *Rlmsg = new int[msgsize];


	for(int i = 0; i < msgsize;i++){		//fill buffers with some data
											//simple sequence in this case
		Slmsg[i] = i+1;
		Srmsg[i] = i+2;
		Rlmsg[i] = i+3;
		Rrmsg[i] = i+4;

	}

	MPI_Status stat[4];
	MPI_Request req[4];

	int r;

	// initiate asynchronous communication

	MPI_Isend(Slmsg, msgsize, MPI_INT, left, 1, MPI_COMM_WORLD, &req[0]);
	
	MPI_Isend(Srmsg, msgsize, MPI_INT, right, 2, MPI_COMM_WORLD, &req[1]);

	MPI_Irecv(Rlmsg, msgsize, MPI_INT, left, 2, MPI_COMM_WORLD, &req[2]);

	MPI_Irecv(Rrmsg, msgsize, MPI_INT, right, 1, MPI_COMM_WORLD, &req[3]);
	
	MPI_Waitall(4, req, stat);			//wait till everyone is done
	


	cout << "Received " << msgsize << "*" << Rlmsg[0] << " from rank " << left <<  endl; 
	cout << "Received " << msgsize << "*" << Rrmsg[0] << " from rank " << right << endl;

	MPI_Finalize();

	delete[] Slmsg;
	delete[] Srmsg;
	delete[] Rrmsg;
	delete[] Rlmsg;

 	return 0;
}

