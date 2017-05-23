/***********************************************
*
*  File Name:       PerfMeasure.h
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Class for process performance measurement
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            24.2.2017
*
***********************************************/

#ifndef __DLB_PERF_MEASURE_H__
#define __DLB_PERF_MEASURE_H__

#include <mpi.h>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <cmath>

#include <chrono>
#include <thread>

#include <Logger/Logger.h>

using std::cout;
using std::endl;
using std::vector;

namespace DLB {

typedef enum method {AVERAGE } AGR_METHOD;

/**
 * @brief Class for process performance measurement
 *
 */

class PerfMeasure
{

public:

    // time spend in actual iteration
    double iter;
    double io;
    double balance;

    double iterTotal;
    double iterAvg;
    double ioTotal;
    double balTotal;
    unsigned sleepTotal;
    double last;

    bool once;
    bool onceMult;

    // base time for delay
    double sleepfor;

    /**
     * @brief Period of sending perfomrance report
     *           expressed in number of iterations
     */
    unsigned period;

    // defines, whether measurement was started
    // in past or not

    unsigned iterCounter;

    // true if iterStart called
    bool running;

    // if true, period of iterations elapsed
    // imbalance detection should be performed
    bool periodEl; 


    // history of measurement for basic statistics
    vector<double> history;

    int rank;
    int worldSize;

    
    PerfMeasure(int rank, int worldSize, unsigned period = 10):
    balance(0.0),
    iterTotal(0.0),
    iterAvg(0.0),
    ioTotal(0.0),
    balTotal(0.0),
    sleepTotal(0.0),
    once(true),
    onceMult(true),
    period(period),
    iterCounter(0),
    periodEl(false),
    rank(rank),
    worldSize(worldSize)
    {

    }

    ~PerfMeasure() {}

    void iterStart(void)
    {
        iter = MPI_Wtime();
        running = true;
    }

    void iterStop(void)
    {
        // if not started, measurement will fail
        if(!running)
            throw std::runtime_error("PerfMeasure.stop called before start");
    
        running = false;

    
        //time elapsed from start() call
        iter = MPI_Wtime() - iter;
        iterTotal += iter;
        last = iter;
        history.push_back(iter);
    
        //number of iterations
        iterCounter++; 

        if(iterCounter == 5){
            sleepfor = agregate(AVERAGE);
            once = false;
            // cout << "once" << endl;
        }
    
        if(iterCounter %  period == 0){
            iterAvg = agregate(AVERAGE);
            periodEl = true;
        }
    }

    void ioStart(void)
    {
        io = MPI_Wtime();
    }
    void ioEnd()
    {
        io = MPI_Wtime() - io;
        ioTotal += io;
    }

    
    void balStart(void)
    {
        balance = MPI_Wtime();
    }

    void balStop(void)
    {
        balance = MPI_Wtime() - balance;
        balTotal += balance;
    }

    double total(void) const
    {
        return io + balance + iter;
    }

    bool periodElapsed(void) const {return periodEl; }


    void sleep(unsigned ms)
    {
        std::chrono::milliseconds msec(ms);
        // double now = MPI_Wtime();
        // cout <<"Wtime: " << now << " " ;
        std::this_thread::sleep_for(msec);
        // cout << MPI_Wtime() << " -> " <<  MPI_Wtime() - now << endl;
    }

    /**
     * Causes delay long enough, to be detected
     */

    void imbalDelay(bool delay, unsigned iter, unsigned nIterations, float multiply)
    {        
        if(once)
            return;

        // unsigned len = nIterations / 4; // length of delay
        unsigned stime = sleepfor * 1000 * multiply;

        if(delay && (iter > nIterations / 4)){
        // if(delay && (iter > nIterations / 10 && iter < (nIterations / 10  + len ) )){

            // if(onceMult){
                // cout << "stime " << stime << endl;
                // cout << "len " << len << endl;
// 
                // onceMult = false;
            // }
// 
            sleep(stime);

            sleepTotal += stime;
        }
    }

    void singleDelay(unsigned mult = 20)
    {
        if(iterCounter == 1){
            double i = history.back();
            sleep(i * mult * 1000);
        }

    }



    /**
    * @brief print measured results
    * @details [long description]
    */

    void print(void) const
    {
        cout << endl;
        // print measures from history vector
        cout << "PerfMeasure: process rank " << rank <<  " history len " << history.size() << endl;
    
        bool once = true;
        for(auto i: history){
    
            if(once){
                cout << i;
                once = false;
            }else{
                cout << ", " << i;
            }
        }
        cout << endl;
    }

    double getAgreg(void)
    {
        return agregate(AVERAGE);
    }
    
    double agregate(AGR_METHOD meth) const
    {
        double result = 0;

        if(history.size() == 0)
            throw std::runtime_error("PerfMeasure::agregate : history empty");

        
        if(meth == AVERAGE){

            result = accumulate(history.begin(), history.end(), 0.0) / history.size();

        }

        return result;
    }

    void reset(void)
    {
        periodEl = false;
        // once = true; // store new iteration time

        history.clear();


    }
};

} // DLB nspace end

#endif