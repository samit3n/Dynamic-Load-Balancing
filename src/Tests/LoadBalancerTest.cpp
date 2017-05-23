#include <iostream>

#include "DynamicBlockDescriptor.h"
#include "LoadBalancer.h"
#include "Asserts.h"
#include "Dims.h"
#include "TileDescriptor.h"

#include <iomanip>
#include <vector>

using namespace DLB;

class LBTest : public LoadBalancer
{

public:

	LBTest(int rank, int edge, int world, Dims o):
	LoadBalancer(rank, edge, world, o)
	{
	}
};


int main()
{

	
	LoadBalancer lb(0, 128, 16, Dims(8,8) );

	vector<TileDescriptor> * t = lb.regularTiles();

	vector<float> times= {1.0,1.0,1.0, 2.0};

	auto vect = lb.splitByPerform(times);

	cout << "sizes: ";
	for(auto v: vect)
		cout << v << " ";

	cout << endl;

}