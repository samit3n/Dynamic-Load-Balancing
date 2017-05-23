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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE TileDescriptorTestBench

#include "../../Sources/DLB/TileDescriptor.h"
#include <boost/test/unit_test.hpp>

using DLB::TileDescriptor;
using  TD = DLB::TileDescriptor;

BOOST_AUTO_TEST_CASE(neighborhood)
{
	// simple 4 tiles
	TileDescriptor t0(0,0 , 64, 64);
	TileDescriptor t1(64,0,64,64);
	TileDescriptor t2(0,64, 64,64);
	TileDescriptor t3(64,64,64,64);

	BOOST_CHECK(t0.isNeighbor(t1));
	BOOST_CHECK(! t1.isNeighbor(t2));
	BOOST_CHECK(! t0.isNeighbor(t3));
	BOOST_CHECK(t2.isNeighbor(t3));
	BOOST_CHECK(t1.isNeighbor(t3));


	BOOST_CHECK_EQUAL(t0.getSharedEdge(t1), TD::RIGHT);
	BOOST_CHECK_EQUAL(t0.getSharedEdge(t2), TD::BOTTOM);

	BOOST_CHECK_EQUAL(t0.getOverlapOffset(t2), (unsigned) 128 );
	BOOST_CHECK_EQUAL(t0.getOverlapOffset(t1), (unsigned) 64 );
	BOOST_CHECK_EQUAL(t3.getOverlapOffset(t2), (unsigned) 192 );



}

BOOST_AUTO_TEST_CASE(size_mix)
{
	// different sizes
	TileDescriptor t0(0,0 , 64, 64);
	TileDescriptor t1(64,0,32,32);
	TileDescriptor t2(64,32, 10, 10);
	TileDescriptor t3(10,64, 20,20);


	BOOST_CHECK_EQUAL(t0.getSharedEdge(t1), TD::RIGHT);
	BOOST_CHECK_EQUAL(t0.getSharedEdge(t2), TD::RIGHT);
	BOOST_CHECK_EQUAL(t0.getSharedEdge(t3), TD::BOTTOM);

	BOOST_CHECK(t0.isNeighbor(t1));
	BOOST_CHECK( t0.isNeighbor(t2));
	BOOST_CHECK( t0.isNeighbor(t3));



	BOOST_CHECK_EQUAL(t0.getOverlapOffset(t1), (unsigned) 64 );
	BOOST_CHECK_EQUAL(t0.getOverlapOffset(t2), (unsigned) 96 );
	BOOST_CHECK_EQUAL(t0.getOverlapOffset(t3), (unsigned) 162  );


}


BOOST_AUTO_TEST_CASE(equalop)
{
	TileDescriptor td1(0,0 , 10, 10);
	TileDescriptor td2(0,0 , 10, 10);
	TileDescriptor td3(0,0, 5,5);

	BOOST_CHECK(td1 == td2);
	BOOST_CHECK(!(td1 != td2));

	BOOST_CHECK(td1 != td3);
}


BOOST_AUTO_TEST_CASE(dynamic)
{
	TileDescriptor td1(0,0 , 48, 64);
	TileDescriptor td2(0,0 , 80, 64);

	BOOST_CHECK(td1 == td2);
	BOOST_CHECK(!(td1 != td2));

	BOOST_CHECK_EQUAL(td1.getSharedEdge(td2), TD::RIGHT);
}

