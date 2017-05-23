/***********************************************
*
*  File Name:       TopologyDescriptor.h
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Describes topology made of irregular mesh
*  					of computation blocks.
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/


#ifndef __DLB_TOPOLOGY_DESCRIPTOR_H__
#define __DLB_TOPOLOGY_DESCRIPTOR_H__

#include <mpi.h>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <cstddef>
#include <algorithm>
#include <map>

#include <TileDescriptor.h>
#include <Asserts.h>
#include <BlockData.h>
#include <Logger/Logger.h>
#include <Neighbor.h>

using std::vector;
using std::map;
using std::find;
using std::pair;

namespace DLB {



/**
 * @brief Describes topology of problem domain, holds metadata of all
 * 		  domain blocks in the problem domain
 * 
 * @details Handles all domain metadata and ensures consistency
 * 			in case of dynamic load balancing.
 * 			
 */

class TopologyDescriptor {

public:

	TopologyDescriptor(int rank, int worldSize, size_t edgeSize);	


	~TopologyDescriptor(void);

	/**
	 * @brief Set tiles
	 * @details Sets actual buffer pointer, containing TileDesciptor
	 * 			objects
	 * 
	 * @param tiles - size == worldSize
	 */

	void setTiles(TileMsg * tiles);
	/**
	 * @brief Sets tiles from given vector
	 */

	void setTiles(const vector<TileDescriptor> & tds );
	/**
	 * @brief Sets tile assigned to process to t
	 */
	void setTile(const TileDescriptor & t);


	/**
	 * @brief Getters
	 */

	Dims getPosition(void) const;
	Dims getSize(void) const;

	TileDescriptor & tile(void) {return myTile; }

	unsigned getDomainArea(void) const;

	int getRank(const Dims & position) const;

	bool getTopChanged(void) { return topologyChanged; }

	BlockData getBlockData(void);

	vector<TileDescriptor> & getTiles(void);


	/**
 	* @brief Updates all topology related metadata 
 	* @details [long description]
 	* 
 	*/

	void updateTopology(void);

	/**
	 * @brief Updates middle column metadata
	 */
	void midUpdate(void);


	// TileMessage MPI datatype
	MPI_Datatype TileMsg_t;

	/**
	 * @brief Creates string representation of topology
	 */
	string toString(void);
	string commsToString(void);


protected:

	const int NO_RANK = -1;


	/**
	 * @brief Find all neighbors in topology and store their ranks into neighbors
	 */

	void countNeighborRanks(void);
	/**
	 * @brief Counts halo zone displacements for all neighbors
	 */
	void countDispls(void);

	/**
	 * @brief Initialize neighbor communicators
	 * @detailed Iterates over COMM_WORLD
	 * 		tile participate in MPI_Comm_split if actual rank 
	 * 		neighbor or my own.
	 * 		Every tile has its communicator to neighbors
	 * 		and participated in N other, where N is neighbor count.
	 */
	
	void initComms(void);

	/**
	 * @brief Initialize TileMsg as MPI datatype
	 */

	void initDtypes(void);
    /**
     * @brief Set block data border for halo zone
     *        2 - standard halo zone
     *        4 - edge shared with domain edge -> no halo zone
     *        
     * 
     * @param BlockData beeing generated
     */

    void getBorders(BlockData & bd);
    void setBorders(void);


	bool topologyChanged;

	vector<TileDescriptor> tiles;
	TileMsg * tileBuffer;


	// my TileDescriptor instance picked from array
	// duplicated for pointer safety
	TileDescriptor myTile;

	// true if edge is border of domain
	// false if there are neighbors
	// needed for halo offset counting
	bool borders[4];

	// tiles adjacent to my tile
	// this tiles communicator for scattering
	MPI_Comm myComm;
	// rank in myComm
	int myCommRank;

  	// COMM_WORLD ranks of neighbors
  	vector<int> neighbors;

  	map<int , Neighbor> nData;

  	// displacement for MPI scatter call
  	vector<int> displs;
  	// how many items, should I send / receive with neighbor
  	// sorted in order of neighbors 
  	vector<int>  counts;

	int rank, worldSize; //in COMM_WORLD
	size_t edgeSize;

	 // optional middle
  	bool middle;
  	MPI_Comm COMM_MIDDLE;
  	int midRank;
  	int midSize;



};

} //DLB nspace end

#endif