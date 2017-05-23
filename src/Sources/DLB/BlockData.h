  /***********************************************
*
*   File Name:    BlockData.cpp

*  Project:     
*  
*  Description:  Class holds all necessary information
*                to perform computation.
*
*  Author:      Vojtech Dvoracek
*   Email:      xdvora0y@stu.fit.vutbr.cz
*   Date:       29.3.2017
*
***********************************************/

#ifndef __BLOCK_DATA_H__
#define __BLOCK_DATA_H__

#include <vector>
#include <map>
#include <mpi.h>

#include <TileDescriptor.h>
#include <Neighbor.h>

using std::vector;
using std::map;

namespace DLB {

/**
 * @brief   Class holding data related to actual block
 * 
 * @details This class is updated by balancing methods
 *          Holds data arrays and sizes.
 *          Contains neighbor metadata as well.
 * 
 * 
 */

class BlockData
{

public:

    BlockData() {}

    BlockData(  float * oldTemp, float * newTemp, float * domParams, int * domMap,    // data arrays
                int top, int bottom, int left, int right,  // block bounds
                bool middle // is in middle column?
            ):

    oldTemp(oldTemp),
    newTemp(newTemp),
    domParams(domParams),
    domMap(domMap),
    top(top),
    bottom(bottom),
    left(left),
    right(right),
    middle(middle)
    {
        neighbors = NULL;
        nData = NULL;
        displs = NULL;
        counts = NULL;
    }

    ~BlockData(){}

    // data pointers

    float * oldTemp, * newTemp, *domParams;
    int * domMap;

    // data bounds

    unsigned top, bottom, left, right;

    // halo flag
    bool topF, bottomF, leftF, rightF;

    // each tile will have 1 communicator per neighbor
    // and one communicator containing all its neighbors
    // it is N+1 communicators, where N is amount of neighbors
    MPI_Comm myComm;
    int myCommRank;

    vector<int> * neighbors;

    map<int, Neighbor> * nData;

    vector<int> * displs;    // MPI scatter
    vector<int> * counts;

    // optional middle
    bool middle;
    MPI_Comm COMM_MIDDLE;
    int midRank;
    int midSize;


    /**
     * @brief Check pointer validity
     * @details [long description]
     * @return [description]
     */
    bool isValid(void) const
    {
        if( ! (nData && neighbors && displs && counts && oldTemp && newTemp && domParams && domMap))
            return false;
    
        return true;
    }

};

inline ostream& operator<<(ostream& os, const BlockData & b)
{
    return os << "borders (T,R,B,L): " << b.top << " " << b.right << " " << b.bottom << " " << b.left << endl;

}

/**
 * @brief Class describing temporary block and related memory 
 *        during data migration.
 * @details tile holds position and size of new block
 */

class TempBlock {

public:

    TempBlock(void){

        temp = NULL;
        params =  NULL;
        map = NULL;
    
    }

    // new temporary will become
    // actual tile afte data migration
  
    TileDescriptor tile;
    // temporary block data arrays
    // with size of new assigned block

    float * temp;
    float * params;
    int * map;

};


} //DLB nspace end

#endif