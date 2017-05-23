/***********************************************
*
*     File Name:        TileDescriptor.cpp

*    Project:         DIP - Dynamic Load Balancing in HPC Applications
*    
*    Description:    
*
*    Author:         Vojtech Dvoracek
*     Email:            
*     Date:              7.3.2017
*
***********************************************/


#include "TileDescriptor.h"

using TileDescriptor = DLB::TileDescriptor;
using Dims = DLB::Dims;
using TileMsg = DLB::TileMsg;


TileDescriptor::TileDescriptor(int rank,unsigned posx, unsigned posy, unsigned dimx, unsigned dimy) :
position(posx, posy),
size(dimx, dimy),
rank(rank)
{
    setHostNumber(getHostname());
}

TileDescriptor::TileDescriptor(unsigned posx, unsigned posy, unsigned dimx, unsigned dimy) :
position(posx, posy),
size(dimx, dimy)
{
    setHostNumber(getHostname());
}


TileDescriptor::TileDescriptor(const TileMsg & data) :
position(data.posx, data.posy),
size(data.dimx, data.dimy),
rank(data.rank),
hostNumber(data.hostNumber)
{
    setHostNumber(getHostname());
}


ostream& DLB::operator<<(ostream& os, const TileDescriptor& obj)
{
    os << "rank: " << obj.getRank() << " pos: " << obj.getPosition();
    os << " size: " << obj.getSize() << " hostNumber: " << obj.getHostNumber() << endl;

    return  os;
}


/**
 * @brief Overloaded equality operator
 * @details compares two tile descriptors based on positon and dimmension
 * 
 * @param t1 left operand
 * @param t2 right operand
 * 
 * @return bool
 */

bool DLB::operator==(const TileDescriptor & left, const TileDescriptor & right)
{
    Dims pos1 = left.getPosition();
    Dims dim1 = left.getSize();

    Dims pos2 = right.getPosition();
    Dims dim2 = right.getSize();

    if(pos1.x == pos2.x && pos1.y == pos2.y){
        if(dim1.x == dim2.x && dim1.y == dim2.y){
            return true;
        }
    }
    return false;
}

/**
 * @brief Set data from received msg
 * 
 * @param data TileMsg object
 */
void TileDescriptor::setData(const TileMsg & data)
{
    position.x =  data.posx;
    position.y = data.posy;

    size.x = data.dimx;
    size.y = data.dimy;
}

TileMsg TileDescriptor::getData(void)
{
    TileMsg data;

    data.posx = position.x;
    data.posy = position.y;

    data.dimx = size.x;
    data.dimy = size.y;

    data.rank = rank;
    data.hostNumber = hostNumber;

    return data;
}

/**
 * @brief Tile position getter
 * @return actual position in domain points
 */

Dims TileDescriptor::getPosition(void) const
{
    return position;
}
/**
 * @brief Tile size getter
 * @return size in domain points
 */
Dims TileDescriptor::getSize(void) const
{
    return size;
}


void TileDescriptor::setSize(const Dims & size)
{
    this->size.x  = size.x;
    this->size.y = size.y;
}


void TileDescriptor::setPosition(const Dims & pos)
{
    this->position.x = pos.x;
    this->position.y = pos.y;
}

bool TileDescriptor::isMiddle(unsigned mid)
{   
    return isInRange(position.x, (position.x + size.x) - 1, mid); 
}


/**
 * @brief Returns tile dimensions extended by haloSize
 * @details [long description]
 * 
 * @param haloSize - size of halo zone
 */
Dims TileDescriptor::getExtSize(int haloSize) const
{
    return Dims(size.x + 2*haloSize, size.y + 2*haloSize);
}

/**
 * @brief Block area getter
 * @return [description]
 */

unsigned TileDescriptor::getArea(void) const
{
    return size.x * size.y;
}

unsigned TileDescriptor::getExtArea(void) const
{

    return (size.x + 2*HALO_WIDTH) * (size.y + 2*HALO_WIDTH);

}

/**
 * @brief Checks whether x in range of <start, stop>
 * 
 * @param start - beginninge
 * @param stop - end
 * @param x - tested value
 * @return true if x >= start and <= stop
 */

inline bool TileDescriptor::isInRange(int start, int stop, int x) const
{
    return x >= start && x <= stop;
}

/**
 * @brief Computes shared edge if there is any
 * @details Checks whether actual tile with given tile
 *             share some edge.
 * 
 * @param tile 
 * @return TEdge - some edge or NONE
 */

TileDescriptor::TEdge TileDescriptor::getSharedEdge(const TileDescriptor & tile) const
{

    Dims tpos = tile.getPosition();
    Dims tsize = tile.getSize();
    Dims tend(tpos.x + tsize.x -1, tpos.y + tsize.y -1);

    Dims end(position.x + size.x -1 , position.y + size.y -1);


    if( isInRange(position.x, end.x, tpos.x) || isInRange(position.x, end.x, tend.x) ||
        isInRange(tpos.x, tend.x, position.x) || isInRange(tpos.x, tend.x, end.x ) ){
        // tile starts or ends in range of x-axis

        if((end.y + 1) == tpos.y){
            return TEdge::BOTTOM;

        }else if(position.y == (tend.y + 1) ){
            return TEdge::TOP;


        }else{
            return TEdge::NONE;
        }

    }else if(isInRange(position.y, end.y, tpos.y) || isInRange(position.y, end.y, tend.y ) ||
            isInRange(tpos.y, tend.y, position.y) || isInRange(tpos.y, tend.y, end.y)){
        // tile starts or ends in range of y-axis

        if( (end.x + 1) == tpos.x){
            return TEdge::RIGHT;
        }else if(position.x == (tend.x + 1) ){
            return TEdge::LEFT;


        }else{
            return TEdge::NONE;
        }

    }else{
        return TEdge::NONE; //digonal or completely distant
    }
}

/**
 * @brief Checks whether two tiles share some edge
 * 
 * @return true if they do, false otherwise
 */

bool TileDescriptor::isNeighbor(const TileDescriptor & tile) const
{
    if(getSharedEdge(tile) != TEdge::NONE){
        return true;
    }
    return false;
}

/**
 * @brief Computes halo offset for given tile
 * @details Halo zones are linearized into array of size HALO_WIDHT * (2*SIZE_X + 2*SIZE_Y).
 *             Data will be send to neighbors using scatterv.
 *             Offset is counted as number of points from topleft corner up to position
 *             where shared edge with given tile begins. It can be used directly to scatterv.
 * 
 * @param tile - neighbor tile sharing halo zone
 * @param cnt - number of shared points
 * @return offset
 */
unsigned TileDescriptor::getOverlapOffset(const TileDescriptor & tile, unsigned & cnt ) const
{
    Dims tpos = tile.getPosition(); // absolute position
    Dims tsize = tile.getSize();    
    Dims tend(tpos.x + tsize.x , tpos.y + tsize.y ); // absolute tile end

    Dims end(position.x + size.x, position.y + size.y); //absolute end

    TEdge shared = getSharedEdge(tile);
    unsigned offset = 100;

    // offset is computed based on which edge is shared and
    // whether the neighbor tile begins of ends on the edge.

    // cout << "shared " << shared << endl;

    if(shared == TOP){
        // cout << "top" << endl;

        // neighbor begins and ends inside
        if(tpos.x >= position.x && tend.x <= end.x){

            offset = tpos.x - position.x;  // (end.x - tend.x);
            cnt = tend.x - tpos.x;

        // begins inside, ends somewheres
        }else if(tpos.x >= position.x && tend.x > end.x){

            offset = tpos.x - position.x;
            cnt = end.x - tpos.x;

        // begins somewhere, end inside
        }else if(tpos.x < position.x && tend.x <= end.x){

            offset = 0;
            cnt = tend.x - position.x;

        }else{ // whole border

            offset = 0;
            cnt = size.x;
        }

    }else if(shared == RIGHT){
        // cout << "right" << endl;

         // neighbor begins and ends inside
        if(tpos.y >= position.y && tend.y <= end.y){

            offset = size.x + ( tpos.y - position.y);
            cnt = tend.y - tpos.y;

         // begins inside, ends somewhere
        }else if(tpos.y >= position.y && tend.y > end.y){

            offset = size.x + (tpos.y - position.y);
            cnt = end.y - tpos.y;

        // begins somewhere, end inside
        }else if(tpos.y < position.y && tend.y <= end.y ){

            offset = size.x;
            cnt = tend.y - position.y;

        // whole side
        }else{

            offset = size.x;
            cnt = size.y;
        }


    }else if(shared == BOTTOM){
        // cout << "bottom" << endl;

        
        if(tpos.x >= position.x && tend.x <= end.x){
            offset = size.y + (tpos.x - position.x); // (end.x - tend.x);
            cnt = tend.x - tpos.x;

        }else if(tpos.x >= position.x && tend.x > end.x){

            offset = size.y + (tpos.x - position.x);
            cnt = end.x - tpos.x;

        }else if(tpos.x < position.x && tend.x <= end.x){

            offset = size.y;
            cnt = tend.x - position.x;
        }else{ // whole border

            offset = size.y;
            cnt = size.x;
        }

        // add offset of TOP and BOTTOM in buffer
        offset += size.x  + size.y;

    }else if(shared == LEFT){
        // cout << "left" << endl;

        if(tpos.y >= position.y && tend.y <= end.y){

            offset = tpos.y - position.y;
            cnt = tend.y - tpos.y;

        }else if(tpos.y >= position.y && tend.y > end.y){

            offset = tpos.y - position.y;
            cnt = end.y - tpos.y;

        }else if(tpos.y < position.y && tend.y <= end.y ){

            offset = 0;
            cnt = tend.y - position.y;

        }else{

            offset = 0;
            cnt = size.y;
        }

        // add offset of TOP and BOTTOM in buffer
        offset += size.x + size.y;

    }else{

        throw std::runtime_error("TileDescriptor::getOverlapOffset error, tile is not neighbor");
    }

    return offset;

}

/**
 * @brief Standard UNIX gethostname() wrapper
 * @return hostname or throws runtime_error
 */

string TileDescriptor::getHostname(void)
{
    const int STRING_SIZE = 50;

    char temp[STRING_SIZE];
    
    for(int i = 0; i < STRING_SIZE; i++){
        temp[i] = 0;
    }


    int err = gethostname(temp, STRING_SIZE);

    if(err){
        stringstream ss;
        ss << "TopologyDescriptor::getHostname - hostname resolution failed, errno = " << err << endl;
        throw std::runtime_error(ss.str()) ;

    }

    string hostname(temp, strlen(temp));
    hostname.erase(hostname.begin() + hostname.size(), hostname.end());


    // delete[] temp;

    return hostname;
}


/**
 * @brief Resolves integer hostname identifier from hostname
 * @details ID is constructe as all decimal numbers in hostname.
 *             Not general, you're welcome to reimplement.
 * 
 * @param hostname [description]
 * @return host ID
 */

int TileDescriptor::setHostNumber(string hostname)
{

    hostNumber = (int) hash<std::string>{}(hostname);
    return hostNumber;
} 


string TileDescriptor::toString(void)
{
    stringstream ss;
    ss << "Position: " << position << endl;
    ss << "Size: " << size << endl;
    ss << "Rank (world): " << rank << endl;
    ss << "Host number: " << hostNumber << endl;

    return ss.str();
}
