#include <iostream>
#include <mpi.h>

#include "DynamicBlockDescriptor.h"
#include "Asserts.h"
#include "Dims.h"

#include <vector>

using namespace std;
using DLB::DynamicBlockDescriptor;
using DLB::Dims;

class DBDTest : public DynamicBlockDescriptor {

	public:

		DBDTest( int rank, int worldSize, size_t edgeSize, Dims objSize):
		DynamicBlockDescriptor(rank, worldSize, edgeSize, objSize)
		{}

		void test(void){

			vector<Dims> v;
			Dims tmp;

			cout << "cols "<< cols << " rows " << rows << endl;


			for(unsigned i = 0; i < 256;i++){
				tmp = getCoordsByGID(i);
				cout << i << " => " <<  tmp << endl;
				v.push_back(tmp);

			}
			cout << "gids" << endl;

			for(auto t : v){
				cout << getGIDByCoords(t) << endl;
			}


		}


};

int main(int argc, char * argv[])
{

	MPI_Init(&argc, &argv);

	int rank;

	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	DBDTest dbd(rank, 2 , 128, Dims(8,8));

	if(rank == 0 ){
		cout << "edgeSize " << 128 << endl;
		cout << dbd.toString();
		dbd.test();

	} 




	MPI_Finalize();

	return 0;
}

