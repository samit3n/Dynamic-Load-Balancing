#include <mpi.h>
#include <zoltan_cpp.h>
#include <iostream>

#include <vector>

using namespace std;



int RootBehav(Zoltan *zz)
{
	int * arr;

	arr = new int[100];

	for(int i = 0 ; i < 100;i++) arr[i] = i;

	int changes, num_imp, num_exp, *imp_procs, *exp_procs;
	int *imp_to_part, *exp_to_part;
	int num_gid_entries, num_lid_entries;
	ZOLTAN_ID_PTR imp_global_ids, exp_global_ids;
	ZOLTAN_ID_PTR imp_local_ids, exp_local_ids;

	zz->LB_Partition(
		changes,
		num_gid_entries,
		num_lid_entries,
    	num_imp,
    	imp_global_ids,
    	imp_local_ids,
    	imp_procs,
    	imp_to_part,
    	num_exp,
    	exp_global_ids,
    	exp_local_ids,
    	exp_procs,
    	exp_to_part
    	); 


	return 0;

}

int OthersBehav(Zoltan * zz)
{
	int *arr;

	arr = new int[25];
	for(int i = 0 ; i < 25;i++) arr[i] = i;




	return 0;
}


int main(int argc, char * argv[])
{




	MPI_Init(&argc, &argv);
	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	cout << rank << ": " << size << endl;

	float version;

	Zoltan_Initialize(argc, argv, &version);
	//Dynamically create Zoltan object.

	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);
	zz->Set_Param("LB_METHOD", "RCB");
	zz->Set_Param("RBC_Reuse", "TRUE");


	if(rank == 0){

		
		RootBehav(zz);
		
	}else{

		OthersBehav(zz);
	}


	
	delete zz;
	MPI_Finalize();

}
