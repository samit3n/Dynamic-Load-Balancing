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

#include "Dims.h"


ostream& DLB::operator<<(ostream& os, const Dims& obj)
{
    return os << "[x,y] : ["<< obj.x << "," <<  obj.y << "]";
}


bool DLB::operator==(const Dims & d1, const Dims & d2)
{
    if(d1.x == d2.x){
        if(d1.y == d2.y)
            return true;
    }

    return false;
}
    
bool DLB::operator!=(const Dims & d1, const Dims & d2)
{
    return ! (d1 == d2);
}

DLB::Dims DLB::operator-(const Dims & d1, const Dims & d2)
{
    Dims d;
    d.x = d1.x - d2.x;
    d.y = d1.y - d2.y;

    return d;
}

DLB::Dims DLB::operator+(const Dims & d1, const Dims & d2)
{

    Dims d;
    d.x = d1.x + d2.x;
    d.y = d1.y + d2.y;

    return d;

}