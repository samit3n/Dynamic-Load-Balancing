
#include <iostream>
#include <unistd.h>
#include <cstddef>
#include "../Sources/DLB/TileDescriptor.h"

#include "../Sources/DLB/TileMsg.h"

using namespace DLB;

using std::cout;


int main()
{


	// TileDescriptor neigh(0,10,5,5);

	// cout << td.getOverlapOffset(neigh) << endl;

    /*
    cout << ( (int *) (&tm.posx) -  (int *) (&tm)) * sizeof(int) << endl;
    cout << ( (int *) (&tm.posy) -  (int *) (&tm)) * sizeof(int) << endl;
    cout << ( (int *) (&tm.dimx) -  (int *) (&tm)) * sizeof(int) << endl;
    cout << ( (int *) (&tm.dimy) -  (int *) (&tm)) * sizeof(int) << endl;
    cout << ( (int *) (&tm.rank) -  (int *) (&tm)) * sizeof(int) << endl;
    cout << ( (int *) (&tm.hostNumber) -  (int *)(&tm)) * sizeof(int)  << endl;
    */

    TileDescriptor t0(0,0, 48, 64);
    TileDescriptor t1(48,0, 80, 64);
    TileDescriptor t2(0,64, 64,64 );
    TileDescriptor t3(64,64,64,64);

    cout << t0.getSharedEdge(t1) << endl;
    cout << t0.getSharedEdge(t2) << endl;
    cout << t0.getSharedEdge(t3) << endl;

    cout << t1.getSharedEdge(t0) << endl;
    cout << t1.getSharedEdge(t2) << endl;
    cout << t1.getSharedEdge(t3) << endl;

    cout << t1.isNeighbor(t0) << endl;
    cout << t1.isNeighbor(t2) << endl;
    cout << t1.isNeighbor(t3) << endl;
    // cout << td0.getSharedEdge(td1) << endl;








}