/***********************************************
*
*  File Name:       Dims.cpp
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Dimensions holding class
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/

#ifndef __DIMS_H__
#define __DIMS_H__

#include <iostream>
#include <ostream>

using std::ostream;

namespace DLB {

/**
 * @brief Simple class holding dimensions (unsigned int)
 * @details Stores x and y dimension of type unsigned int
 * 
 * Useful for storing coordinates.
 * 
 */
	
class Dims {

public:

    Dims(unsigned x, unsigned y): x(x), y(y) {}
    Dims() {}

    unsigned x,y;

};


//related operators

ostream& operator<<(ostream& os, const Dims& obj);

bool operator==(const Dims & d1, const Dims & d2);

bool operator!=(const Dims & d1, const Dims & d2);

Dims operator-(const Dims & d1, const Dims & d2);
Dims operator+(const Dims & d1, const Dims & d2);


} //DLB nspace end

#endif