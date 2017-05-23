/***********************************************
*
*     File Name:        TileDescriptor.h

*    Project:         DIP - Dynamic Load Balancing in HPC Applications
*    
*    Description:    
*
*    Author:         Vojtech Dvoracek
*     Email:            
*     Date:              7.3.2017
*
***********************************************/

#ifndef __DLB_TILE_DESCRIPTOR_H__
#define __DLB_TILE_DESCRIPTOR_H__

#include <iostream>
#include <stdexcept>
#include <sstream>

#include <unistd.h>
#include <functional>
#include <cstring>

#include <Dims.h>
#include <TileMsg.h>

// using namespace std;

using std::string;
using std::endl;
using std::cout;
using std::endl;
using std::stringstream;
using std::hash;

namespace DLB {


const unsigned HALO_WIDTH = 2;


/**
 * @brief Describes one rectangular tile and its position in topology
 * 
 * @details Position si given as number of domain points from upper-left corner.
 * 
 */

class TileDescriptor {

public:

    // constants for rectangle edges
    // !! array indices, check TopologyDescriptor::(get/set)Borders()
    typedef enum ed { TOP = 0, RIGHT, BOTTOM, LEFT, NONE } TEdge;

    TileDescriptor(int rank, unsigned posx, unsigned posy, unsigned dimx, unsigned dimy);
    TileDescriptor(unsigned posx, unsigned posy, unsigned dimx, unsigned dimy);
    TileDescriptor(const TileMsg & data);
    TileDescriptor(void) {}

    //output overload
    // friend ostream& operator<<(ostream& os, const TileDescriptor& obj)
    // {
        // return os << "Tile:" << endl << "Position " <<  obj.position << endl << "Size " << obj.size << endl;
    // }

    
    void setData(const TileMsg &data);
    TileMsg getData(void);

    Dims getPosition(void) const;
    Dims getSize(void) const ;

    void setPosition(const Dims & pos);
    void setSize(const Dims & size);
    // width * height
    unsigned getArea(void) const;

    unsigned getExtArea(void) const;
    Dims getExtSize(int haloSize = 2) const;

    bool isMiddle(unsigned mid);

    //offset of middle column
    unsigned middleCol(void);

    



    int getRank(void) const {return rank; }
    void setRank(int rank) {this->rank = rank; }

    bool isNeighbor(const TileDescriptor & tile) const;
    TEdge getSharedEdge(const TileDescriptor & tile) const;
    unsigned getOverlapOffset(const TileDescriptor & tile,  unsigned & cnt) const;


    void setHostNumber(int hostNumber) {this->hostNumber = hostNumber; }

    // sets and stores hostNumber computed from hostname
    // implemented as hostname std::hash (C++11)
    int setHostNumber(string hostName);

    int getHostNumber(void) const {return hostNumber; }



    static string getHostname(void);

    string toString(void);

protected:

    Dims position, size;

    //MPI rank
    int rank;

    // number of host which is mapped to this process
    // based on hostname, creates hostname <> MPI rank mapping
    // useful for optimizing communications

    int hostNumber;

    // simple interval check
    inline bool isInRange(int start, int stop, int x) const;


};

/**
 * @brief Overloaded equality/inequality operators
 * @details compares two tile descriptors based on positon and dimmension
 */
bool operator==(const TileDescriptor & left, const TileDescriptor & right);

//find tile by rank
inline bool operator==(const TileDescriptor & left, const int rank)
{
    return left.getRank() == rank;
}

inline bool operator!=(const TileDescriptor & left, const TileDescriptor & right)
{
    return !(left == right);
}

ostream& operator<<(ostream& os, const TileDescriptor& obj);


} //DLB nspace end



#endif