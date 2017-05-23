/**
* MPI ring communication
* 
* Ring communication is defined as set of processes
* in ring topology, where every process want's to send
* message to both his neighbors.
* 
* This variant implenets synchronous variant, which leads
* to deadlock when MPI switches from buffered to synch. transfer.
* 
* Author: Vojtech Dvoracek
* email: xdvora0y@stud.fit.vutbr.cz
*/



#include <mpi.h>
#include <iostream>

using namespace std;

const int MYTAG = 123;				//generic message tag

int main(int argc, char ** argv)
{

	if(argc < 2){
		cerr << "Give msg len";
		return 0;
	}
	
	int rank, left, right, size;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//determine neighbor rank based on mine

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

	// process args, produce error if 0-byte message is given -> nonsense

	int msgsize = std::stoi(argv[1]);

	if(size < 1){
		cerr << "Message too small";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	int *rmsg = new int[msgsize];
	int *lmsg = new int[msgsize];

	for(int i = 0; i < msgsize;i++){

		rmsg[i] = i;
		lmsg[i] = i;
	}

	MPI_Status stat;
	int r;

	// perform synchronous transfer
	
	r = MPI_Send(lmsg, msgsize, MPI_INT, left, 1, MPI_COMM_WORLD);

	r = MPI_Send(rmsg, msgsize, MPI_INT, right, 2, MPI_COMM_WORLD);

	r = MPI_Recv(rmsg, msgsize, MPI_INT, left, 2, MPI_COMM_WORLD, &stat);
	r = MPI_Recv(lmsg, msgsize, MPI_INT, right, 1, MPI_COMM_WORLD, &stat);



	cout << "Received " << msgsize << "*" << lmsg[0] << " from rank " << left <<  endl; 
	cout << "Received " << msgsize << "*" << rmsg[0] << " from rank " << right << endl;

	MPI_Finalize();

	delete[] rmsg;
	delete[] lmsg;

 	return 0;
}

