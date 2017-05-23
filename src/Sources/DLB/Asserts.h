/***********************************************
*
*  File Name:       Asserts.h
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Class holds all necessary information
*                   to perform computation.
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/


#ifndef __DLB_ASSERTS_H__
#define __DLB_ASSERTS_H__

#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <mpi.h>

using std::runtime_error;
using std::stringstream;
using std::string;
using std::endl;
using std::cout;
using std::flush;

namespace DLB {



/**
 * @brief Macro generates location info.
 * May be passed to assert.
 * 
 * Uncomment empty version if these macros unsupported
 */

#define LOCATION ,__func__, __FILE__, __LINE__
#define COUTLOC  __func__ << " : " <<  __FILE__ << " : " << __LINE__ << endl

// debug flag
const bool DBG = false;


/**
 * Custom assert functions
 * 
 * Two basic variants of, with or without custom message.
 * 
 */


inline void MPI_assert(int result, const string  & msg)
{
    if(result){
        stringstream ss;
        ss << "MPI call returned " << result << endl << "message => " << msg << endl;
    
        throw runtime_error( ss.str() );
    }
}

inline void MPI_assert(int result)
{
    if(result){

        stringstream ss;
        ss << "MPI call returned " << result << endl;

        throw runtime_error(ss.str());
  }
}

/**
 * @brief Assert versions extended by LOCATION macro
 * @details [long description]
 * 
 * @param result - returned by MPI call
 * @param msg - custom message
 * @param func - __func__ value
 * @param file - __FILE__ value
 * @param line - __LINE__ value
 * 
 * Usage:
 * 
 *  MPI_assert(MPI_Some_call(), "Some call failed" LOCATION);
 *  
 *  Last comma is added by macro for compatibility reasons.
 */

inline void MPI_assert(int result, const string & msg, const char * func , const char * file, int line)
{
    if(result){

        stringstream ss;
        ss << "MPI call returned " << result << endl;
        ss << " in function " << func << " in file " << file  << " on line " << line << endl << "message => " << msg << endl;
     
        throw runtime_error(ss.str());
    }    
}

inline void MPI_assert(int result, const char * func , const char * file, int line)
{
    if(result){

        stringstream ss;
        ss << "MPI call returned " << result << endl;
        ss << " in function " << func << " in file " << file  << " on line " << line << endl ;

        throw runtime_error(ss.str());
    }    
}


inline void synCout(const string & msg, int rank, int worldSize, MPI_Comm comm = MPI_COMM_WORLD )
{

    for(int i = 0; i < worldSize;i++){
        if(rank == i){
           cout << msg << flush;
        }
    }
    MPI_Barrier(comm);
}


} //DLB nspace end


/**
 * @brief Suppress unused variable warnign
 * @details Kind of hacky template to selectively suppress unused variable warnings.
 * 
 * Don't wanna sabotate compiler's optimization efforts, it is especially for callbacks, with
 * given function type where we really don't need all params and don't wanna
 * shut up all warnings.
 * 
 * Found on stackoverflow.com
 */

template<class T> inline void unused( const T& ) { }



#endif


