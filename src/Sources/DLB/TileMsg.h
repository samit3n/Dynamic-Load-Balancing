/***********************************************
*
*  File Name:       TopologyDescriptor.cpp
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     MPI message for TileDescriptor class 
*                   used to exchange topology between processes
*                   (DTO design pattern)
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            7.3.2017
*
*******************************************/

#ifndef __TILE_MSG_H__
#define __TILE_MSG_H__

/**
 * @brief Simple structure holding tile data to be sent over MPI
 */

#include <Dims.h>

namespace DLB {

typedef struct tmsg 
{
    tmsg() {}
    tmsg(unsigned posx, unsigned posy, unsigned dimx, unsigned dimy, int rank, int hostNumber):
        posx(posx), posy(posy), dimx(dimx), dimy(dimy), rank(rank), hostNumber(hostNumber) {}

    tmsg(unsigned posx, unsigned posy, unsigned dimx, unsigned dimy) : posx(posx), posy(posy), dimx(dimx), dimy(dimy) {}

    tmsg(const Dims &pos, const Dims & dim, int rank, int hN) :
    	posx(pos.x), posy(pos.y), dimx(dim.x), dimy(dim.y), rank(rank), hostNumber(hN)
    {}

    unsigned posx, posy;
    unsigned dimx, dimy;

    int rank;
    int hostNumber;

} TileMsg;


/**
* @brief Single model point representation for data migration using Zoltan
*/

typedef struct point {

    point() {}
    point(float temp, float param, int map) : temp(temp), param(param), map(map) {}

    float temp;
    float param;
    int map;

} Point_t;

} //DLB nspace end

#endif