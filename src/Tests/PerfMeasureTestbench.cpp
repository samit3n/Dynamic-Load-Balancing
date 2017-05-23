
#include "../Sources/PerfMeasure.h"

#include <mpi.h>
#include <chrono>
#include <thread>
#include <random>

const int COMPUTE_LEN = 500;
const int PERIOD = 4;
const int DELAY_RANK = 2;

using namespace std;

int main(int argc, char *argv[])
{
	

	int rank, size;
	
	// Initialize MPI
	MPI_Init(&argc, &argv);
	// Get MPI rank and size
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// random generator, seeded by rank
	// to generate pseudo random delay
	default_random_engine generator(rank);
	uniform_real_distribution<double> isDelay(0.0, 1.0);

	// possible delays in miliseconds
	uniform_int_distribution<int> delay(50, 200);

	PerfMeasure pmes(rank, size, PERIOD);


	for(int  i = 0; i < 50; i++){

		pmes.start();
		
		//computation length - fixed
	    this_thread::sleep_for(chrono::milliseconds(COMPUTE_LEN));

	    // random delay
	    if(rank == DELAY_RANK){

	    	if(isDelay(generator) > 0.5){
				this_thread::sleep_for(chrono::milliseconds(delay(generator)));
			}
	    }


	    pmes.stop();
	    
	    if(rank == 0 && i % PERIOD == 1)
	    	pmes.print();


	}	

	if(rank != 0)
    	pmes.print();

    MPI_Finalize();



	return 0;
}