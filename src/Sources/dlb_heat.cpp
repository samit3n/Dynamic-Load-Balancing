/**
 * @file        dlb_heat.cpp
 * @author      Vojtech Dvoracek, Master theesis
 * 
 *              Original sequential version and I/O routines
 *              created by:
 *              
 *              Jiri Jaros, Radek Hrbacek and Filip Vaverka\n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       Parallel heat distribution simulation model
 *              with dynamic load balancing (aka dlb) ability.
 *              Developed as master thesis at FIT VUT Brno.
 *
 * @version     2017
 * @date        10 April 2015, 10:22 (created) \n
 * @date        30 April 2017, 11:08 (last revised) \n
 *
 * @detail
 * This is the main file of the project.
 * 
 * 
 */


#include <mpi.h>
#include <hdf5.h>

#include <string.h>
#include <string>
#include <cmath>
#include <sstream>
#include <iostream>
#include <new>


#include "MaterialProperties.h"
#include "BasicRoutines.h"

// Dynamic Load Balancing files
#include <Asserts.h>
#include <DynamicBlockDescriptor.h>
#include <Logger/Logger.h>
#include <BlockData.h>
#include <Dims.h>
#include <HaloBuffers.h>



const bool DEBUG = true;
const int RANK = 1;


// dynamic load balancing namespace
using namespace DLB;

// constant Zoltan object size dimensions
// may be passed to DynamicBlockDescriptor ctor
// for testing only - value from user params should be used
const Dims OBJSIZE(16,16);

// some std stuff
using std::cout;
using std::endl;
using std::cerr;
using std::stringstream;

//----------------------------------------------------------------------------//
//---------------------------- Global variables ------------------------------//
//----------------------------------------------------------------------------//

/// Temperature data for sequential version.
float *seqResult = NULL;
/// Temperature data for parallel method.
float *parResult = NULL;

/// Parameters of the simulation
TParameters parameters;

/// Material properties
TMaterialProperties materialProperties;



//----------------------------------------------------------------------------//
//------------------------- Function declarations ----------------------------//
//----------------------------------------------------------------------------//

/// Sequential implementation of the Heat distribution
void SequentialHeatDistribution(float                     *seqResult,
                                const TMaterialProperties &materialProperties,
                                const TParameters         &parameters,
                                string                     outputFileName);

/// Parallel Implementation of the Heat distribution (Non-overlapped file output)
void ParallelHeatDistribution(float                     *& parResult,
                              const TMaterialProperties &materialProperties,
                              const TParameters         &parameters,
                              string                     outputFileName);

/// Store time step into output file
void StoreDataIntoFile(hid_t         h5fileId,
                       const float * data,
                       const size_t  edgeSize,
                       const size_t  snapshotId,
                       const size_t  iteration);

/// Store time step into output file using parallel HDF5
void StoreDataIntoFileParallel(hid_t h5fileId,
                               const float *data,
                               const size_t edgeSize,
                               const size_t tileWidth, const size_t tileHeight,
                               const size_t tilePosX, const size_t tilePosY,
                               const size_t snapshotId,
                               const size_t iteration);
// ----------------------
// Some useful functions
// ----------------------


/**
 * @brief Template function copying halo zone
 *        to linear buffer, to be send over MPI.
 *        
 * @details Expects data block extended by halo zone of size = 2
 *          on every side.
 * 
 * @param block - original data block
 * @param buff - allocated buffer with sufficient size ( = 2 * block perimeter)
 * @param ext - size of block including halo zones
 */

template <class T>
void HaloToBuff(T * block, T * buff, Dims ext)
{
    unsigned offset = 0;
    bool once = true;

    //top
    for(unsigned i = 2; i < 4;i++){
      for(unsigned j = 2; j < ext.x-2; j++){

        buff[offset] = block[i * ext.x + j];
        offset += 2;
      }

      if(once){
        offset = 1;
        once = false;
      }
    }
    offset -= 1;

    // right
    for(unsigned j = 2; j < ext.y-2;j++){
      for(unsigned i = ext.x - 4;i < ext.x -2;i++){
        buff[offset] = block[j * ext.x + i];
        offset++;
      }
    }

    //left
    for(unsigned j = 2; j < ext.y-2;j++){
      for(unsigned i = 2 ;i < 4 ;i++){
        buff[offset] = block[j * ext.x + i];
        offset++;
      }
    }

    //bottom
    unsigned tmp = offset;

    for(unsigned i = ext.y - 4; i < ext.y - 2;i++){
      for(unsigned j = 2; j < ext.x-2;j++){

        buff[offset] = block[i * ext.x + j];
        offset += 2;
      }
      offset = tmp + 1;
    }

  
}

/**
 * @brief Copies data from linear buffer to halo zones.
 * 
 * @details Expects data block extended by halo zone of size = 2
 *          on every side.
 * 
 * @param block - original data block
 * @param buff - allocated buffer with sufficient size ( = 2 * block perimeter)
 * @param ext - size of block including halo zones
 */

template <class T>
void BuffToHalo(T * block, T * buff, Dims ext)
{
    unsigned offset = 0;
    bool once = true;

    //top
    for(unsigned i = 0; i < 2;i++){
      for(unsigned j = 2; j < ext.x-2; j++){

        block[i * ext.x + j] = buff[offset];
        offset += 2;
      }
      if(once){
        offset = 1;
        once = false;
      }
    }
    offset -= 1;

    // right
    for(unsigned j = 2; j < ext.y-2;j++){
      for(unsigned i = ext.x - 2;i < ext.x ;i++){
        block[j * ext.x + i] = buff[offset];
        offset++;
      }
    }

    //left
    for(unsigned j = 2; j < ext.y-2;j++){
      for(unsigned i = 0 ;i < 2 ;i++){
        block[j * ext.x + i] = buff[offset];
        offset++;
      }
    }

    // bottom
    unsigned tmp = offset;
    for(unsigned i = ext.y - 2; i < ext.y;i++){
      for(unsigned j = 2; j < ext.x-2;j++){

        block[i * ext.x + j] = buff[offset];
        offset += 2;
      }
      offset = tmp + 1;
    }
}

/**
 * @brief Printable block representation
 * 
 * @details Template function for array block printing
 * 
 * @param block - data
 * @param extSize - block size including halo zones
 * @param rank - which ranks prints it
 * @param extended:
 *         true - inlcuded halo zones
 *         false - core only
 *         
 * @return [description]
 */

template <class T>
string BlockToString(T * block, Dims extSize, int rank, bool extended)
{


    stringstream  ss;

    if(block == NULL)
        throw std::runtime_error("PrintBlock (float *): NULL pointer passed in!");


    if(extended){

        ss << " ========= EXTENDED ARRAY - RANK " << rank << " ========" << endl;
  
        for(int i = 0; i < (int) extSize.y; i++){
            for(int j = 0; j < (int) extSize.x ; j++){

                ss <<  block[i*(extSize.x)+ j] << " ";
            }
            ss << endl;
        }

    }else{
        ss << " ========= DATA ONLY ========" << endl;
  
        for(int i = 2; i < (int) extSize.y - 2; i++){
            for(int j = 2; j < (int) extSize.x - 2 ; j++){

                  ss << block[i*(extSize.x)+ j] << " ";
            }
            ss << endl;
        }
    }

    return ss.str();
}

template<class T>
void PrintBlock(T * block, Dims extSize, int rank, bool extended)
{
  cout << BlockToString<T>(block, extSize, rank, extended);
}



//----------------------------------------------------------------------------//
//------------------------- Function implementation  -------------------------//
//----------------------------------------------------------------------------//


void ComputePoint(float  *oldTemp,
                  float  *newTemp,
                  float  *params,
                  int    *map,
                  size_t  i,
                  size_t  j,
                  size_t  edgeSize,
                  float   airFlowRate,
                  float   coolerTemp)
{
    // [i] Calculate neighbor indices
    const int center    = i * edgeSize + j;
    const int top[2]    = { center - (int)edgeSize, center - 2*(int)edgeSize };
    const int bottom[2] = { center + (int)edgeSize, center + 2*(int)edgeSize };
    const int left[2]   = { center - 1, center - 2};
    const int right[2]  = { center + 1, center + 2};

    // [ii] The reciprocal value of the sum of domain parameters for normalization
    const float frac = 1.0f / (params[top[0]]    + params[top[1]]    +
                            params[bottom[0]] + params[bottom[1]] +
                            params[left[0]]   + params[left[1]]   +
                            params[right[0]]  + params[right[1]]  +
                            params[center]);

    // [iii] Calculate new temperature in the grid point
    float pointTemp = 
        oldTemp[top[0]]    * params[top[0]]    * frac +
        oldTemp[top[1]]    * params[top[1]]    * frac +
        oldTemp[bottom[0]] * params[bottom[0]] * frac +
        oldTemp[bottom[1]] * params[bottom[1]] * frac +
        oldTemp[left[0]]   * params[left[0]]   * frac +
        oldTemp[left[1]]   * params[left[1]]   * frac +
        oldTemp[right[0]]  * params[right[0]]  * frac +
        oldTemp[right[1]]  * params[right[1]]  * frac +
        oldTemp[center]    * params[center]    * frac;

    // [iv] Remove some of the heat due to air flow (5% of the new air)
    pointTemp = (map[center] == 0)
              ? (airFlowRate * coolerTemp) + ((1.0f - airFlowRate) * pointTemp)
              : pointTemp;

    newTemp[center] = pointTemp;
}

/**
 * Sequential version of the Heat distribution in heterogenous 2D medium
 * @param [out] seqResult          - Final heat distribution
 * @param [in]  materialProperties - Material properties
 * @param [in]  parameters         - parameters of the simulation
 * @param [in]  outputFileName     - Output file name (if NULL string, do not store)
 *
 */
void SequentialHeatDistribution(float                      *seqResult,
                                const TMaterialProperties &materialProperties,
                                const TParameters         &parameters,
                                string                     outputFileName)
{

    // [1] Create a new output hdf5 file
    hid_t file_id = H5I_INVALID_HID;
  
    if (outputFileName != "")
    {
      if (outputFileName.find(".h5") == string::npos)
        outputFileName.append("_seq.h5");
      else
        outputFileName.insert(outputFileName.find_last_of("."), "_seq");
      
      file_id = H5Fcreate(outputFileName.c_str(),
                          H5F_ACC_TRUNC,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
      if (file_id < 0) ios::failure("Cannot create output file");
    }


    // [2] A temporary array is needed to prevent mixing of data form step t and t+1
    float * tempArray = new float[materialProperties.nGridPoints];
    // float *tempArray = (float *)_mm_malloc(materialProperties.nGridPoints * 
                                           // sizeof(float), DATA_ALIGNMENT);
       // [3] Init arrays
    for (size_t i = 0; i < materialProperties.nGridPoints; i++)
    {
      tempArray[i] = materialProperties.initTemp[i];
      seqResult[i] = materialProperties.initTemp[i];
    }

    // [4] t+1 values, t values
    float *newTemp = seqResult;
    float *oldTemp = tempArray;

    if (!parameters.batchMode) 
      printf("Starting sequential simulation... \n");
  
    //---------------------- [5] start the stop watch ------------------------------//
    double elapsedTime = MPI_Wtime();
    size_t i, j;
    size_t iteration;
    float middleColAvgTemp = 0.0f;

    // [6] Start the iterative simulation
    for (iteration = 0; iteration < parameters.nIterations; iteration++)
    {
        // [a] calculate one iteration of the heat distribution (skip the grid points at the edges)
        for (i = 2; i < materialProperties.edgeSize - 2; i++)
            for (j = 2; j < materialProperties.edgeSize - 2; j++)
                ComputePoint(oldTemp,
                        newTemp,
                        materialProperties.domainParams,
                        materialProperties.domainMap,
                        i, j,
                        materialProperties.edgeSize, 
                        parameters.airFlowRate,
                        materialProperties.coolerTemp);

        // [b] Compute the average temperature in the middle column
        middleColAvgTemp = 0.0f;
        for (i = 0; i < materialProperties.edgeSize; i++)
            middleColAvgTemp += newTemp[i*materialProperties.edgeSize +
                                materialProperties.edgeSize/2];
        middleColAvgTemp /= materialProperties.edgeSize;

        // [c] Store time step in the output file if necessary
        if ((file_id != H5I_INVALID_HID)  && ((iteration % parameters.diskWriteIntensity) == 0))
        {
            StoreDataIntoFile(file_id,
                        newTemp,
                        materialProperties.edgeSize,
                        iteration / parameters.diskWriteIntensity,
                        iteration);
        }

        // [d] Swap new and old values
        swap(newTemp, oldTemp);

        // [e] Print progress and average temperature of the middle column
        if ((iteration % (parameters.nIterations / 10l)) == 
            ((parameters.nIterations / 10l) - 1l) && !parameters.batchMode)
        {
          printf("Progress %ld%% (Average Temperature %.2f degrees)\n", 
                 iteration / (parameters.nIterations / 100) + 1, 
                 middleColAvgTemp);
        }
    } // for iteration

    //-------------------- stop the stop watch  --------------------------------//  
    double totalTime = MPI_Wtime() - elapsedTime;

    // [7] Print final result
    if (!parameters.batchMode)
        printf("\nExecution time of sequential version %.5f\n", totalTime);
    else
        printf("%s;%s;%f;%e;%e\n", outputFileName.c_str(), "seq",
                                middleColAvgTemp, totalTime,
                                totalTime / parameters.nIterations);   

    // Close the output file
    if (file_id != H5I_INVALID_HID) H5Fclose(file_id);
    
    // [8] Return correct results in the correct array
    if (iteration & 1)
      memcpy(seqResult, tempArray, materialProperties.nGridPoints * sizeof(float));

    // _mm_free(tempArray);
    delete[] tempArray;
} // end of SequentialHeatDistribution
//------------------------------------------------------------------------------

//-------- Definition of structures for topology description ------ //

//serialize halo zone to buffer




/**
 * @brief Exchange all halo zones (temp, params, map)on topology change.
 * 
 */

void ExchangeHaloZones( BlockData &bd,
                        DynamicBlockDescriptor & dbd,
                        float * sendTemp,
                        float * recvTemp,
                        float * sendParams,
                        float * recvParams,
                        const int * cts, //scatter counts 
                        const int * dpsl // scatter displs
                    )
{
    // size for all requests
    // int size = bd.neighbors->size()*3 + 3;

    MPI_Request  * send = new MPI_Request[2];
    MPI_Status * sendSt = new MPI_Status[2];

    MPI_Request  * recv = new MPI_Request[bd.neighbors->size() * 2];
    MPI_Status * recvSt = new MPI_Status[bd.neighbors->size() * 2];

    unsigned nCount = 0;


    if(DBG){
        char  * name = new char[MPI_MAX_OBJECT_NAME];
        int l, size;
        MPI_Comm_get_name( bd.myComm, name, &l);
        MPI_Comm_size( bd.myComm, &size);
        stringstream ss;
        ss << dbd.rank << " sending on comm " << name << " of size " << size <<   endl;
        synCout(ss.str(), dbd.rank, dbd.worldSize, bd.myComm);
    }

    MPI_assert( MPI_Iscatterv(sendTemp, cts, dpsl, MPI_FLOAT,  NULL, 0, MPI_FLOAT, bd.myCommRank, bd.myComm, &(send[0]) ),
       "Iscatter temp send failed" LOCATION );
    MPI_assert( MPI_Iscatterv(sendParams, cts, dpsl, MPI_FLOAT, NULL, 0, MPI_FLOAT, bd.myCommRank, bd.myComm, &(send[1])), 
      "Iscatter params send failed" LOCATION );

    for(auto it = bd.neighbors->begin(); it < bd.neighbors->end();it++){

        if(DBG){
            char  * name = new char[MPI_MAX_OBJECT_NAME];
            int l, size;
            MPI_Comm_get_name( bd.nData->at(*it).comm, name, &l);
            MPI_Comm_size(bd.nData->at(*it).comm, &size); 
            stringstream ss;
            ss << dbd.rank << " receiving "  << " on " << name << " of size " << size<< endl;
            synCout(ss.str(), dbd.rank, dbd.worldSize, bd.nData->at(*it).comm );
        }
        // MPI_Barrier(bd.nData->at(*it).comm);
        MPI_assert( MPI_Iscatterv(NULL, bd.nData->at(*it).scatterCnts, bd.nData->at(*it).scatterDispls, MPI_FLOAT, 
            &(recvTemp[ bd.nData->at(*it).displ ]), bd.nData->at(*it).count, 
            MPI_FLOAT, bd.nData->at(*it).root, bd.nData->at(*it).comm, &(recv[nCount + 0]) ), 
            "Iscatter temp send failed" LOCATION );
        MPI_assert( MPI_Iscatterv(NULL, bd.nData->at(*it).scatterCnts, bd.nData->at(*it).scatterDispls, MPI_FLOAT, 
            &(recvParams[bd.nData->at(*it).displ]), bd.nData->at(*it).count, MPI_FLOAT, 
            bd.nData->at(*it).root, bd.nData->at(*it).comm, &(recv[nCount + 1])), 
            "Iscatter params send failed" LOCATION );
          
        nCount+= 2;
    }
        
    MPI_assert( MPI_Waitall( 2, send, recvSt), "Receive waitall failed" LOCATION);

   
    MPI_assert( MPI_Waitall(bd.neighbors->size()*2, recv, recvSt), "Receive waitall failed" LOCATION);

    delete[] send;
    delete[] sendSt;
    delete[] recv;
    delete[] recvSt;

    return;

}

/**
 * @brief Function computing single simulation step on halo zones only.
 * 
 * @details Increases code readability when communication overlap should be used.
 * 
 * @param bd  - BlockData instance
 * @param dbd - DynamicBlockDescriptor
 * @param airFlowRate - air flow param
 * @param coolerTemp - cooler temperature param
 */
void ComputeHalo(BlockData & bd, DynamicBlockDescriptor & dbd, float airFlowRate, float  coolerTemp)
{

    Dims ext = dbd.getExtSize();
    //top
    if(! bd.topF){
        for(unsigned i = bd.top; i < bd.top + 2;i++)  {
          for(unsigned j = bd.left; j < bd.right;j++){
    
            ComputePoint(bd.oldTemp,
                         bd.newTemp,
                         bd.domParams,
                         bd.domMap,
                         i, j,
                         ext.x, 
                         airFlowRate,
                         coolerTemp
                         );
          }
        }
    }

    // righ
    if(! bd.rightF){
        for(unsigned i = bd.right - 2;i < bd.right;i++){
          for(unsigned j = bd.top; j <bd.bottom;j++){
    
              ComputePoint(bd.oldTemp,
                         bd.newTemp,
                         bd.domParams,
                         bd.domMap,
                         j, i,
                         ext.x, 
                         airFlowRate,
                         coolerTemp
                         );
          }
        }
    }

    //bottom
    if(! bd.bottomF){
        for(unsigned i = bd.bottom - 2;i < bd.bottom ;i++){
          for(unsigned j = bd.left; j < bd.right;j++){
    
              ComputePoint(bd.oldTemp,
                         bd.newTemp,
                         bd.domParams,
                         bd.domMap,
                         i, j,
                         ext.x, 
                         airFlowRate,
                         coolerTemp
                         );
          }
        }
    }

    //left
    if(! bd.leftF){
        for(unsigned i = bd.left ;i < bd.left+2 ;i++){
          for(unsigned j = bd.top; j < bd.bottom;j++){
    
             ComputePoint(bd.oldTemp,
                         bd.newTemp,
                         bd.domParams,
                         bd.domMap,
                         j, i,
                         ext.x, 
                         airFlowRate,
                         coolerTemp
                         );
            }
        }
    }
}


/**
* @function Parallel process behavior
* 
* @details Main simulation code same for all parallel processes
*          A lot of rank-based decision making is made inside DBD etc.
* 
* @param [in] &materialProperties - initial model data,
*             will be scattered amongst processes
* @param [in] parameters  - model parameters obatined from input data
* 
* @param [in] outputFileName - file to print simulation results
* @param [in] rank - MPI rank
* @param [in] size - MPI world size
*/


float * Behavior(   const TMaterialProperties  &materialProperties,
                    const TParameters          &parameters,
                    hid_t                      file_id,
                    int                        rank,
                    int                        size
                    
                    )
{
    if( DBG && rank == 0){

        cout << "Running config: " << endl;
        cout << "Object size: " << Dims(parameters.objDim, parameters.objDim) << endl;
        cout << "Balancing: " << (parameters.balance ? "ON" : "OFF") << endl;
        cout << "Balance period: " << parameters.balancePeriod << endl;
        cout << "World Size: " << size << endl << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(DBG){

      stringstream ss;
      ss << "rank " << rank << " working" << endl;
  
      synCout(ss.str(), rank, size);
    }

    //dimensions per block, all the same at this situation
    DynamicBlockDescriptor dbd(rank,
                               size,
                               parameters.edgeSize,
                               Dims(parameters.objDim, parameters.objDim),
                               parameters.threshold
                               );

    BlockData bd;

    PerfMeasure pm(rank, size, parameters.balancePeriod);

    dbd.zoltanInit();

    // loadInit distinguish between root and others
    // material properties may be empty in others
    bd = dbd.loadInit(materialProperties);

    // float * tempArray = bd.oldTemp;
    // hallo send and receive buffers
    HaloBuffers hb(2*dbd.getHaloLen());

    // store simulation timestamp
    double totalTime = MPI_Wtime();

    float middleColAvgTemp = 0.0f;

    // pack halo zones to buffers
    HaloToBuff<float>(bd.oldTemp, hb.sendTemp, dbd.getExtSize());
    HaloToBuff<float>(bd.domParams, hb.sendParams, dbd.getExtSize());
    HaloToBuff<int>(bd.domMap, hb.sendMap, dbd.getExtSize());

    int * cnts = new int[bd.counts->size()];
    int * displs = new int[bd.displs->size()];

    for(unsigned i = 0; i < bd.counts->size();i++){
        cnts[i] = bd.counts->at(i);
        displs[i] = bd.displs->at(i);
    }

    ExchangeHaloZones(bd, dbd, hb.sendTemp, hb.recvTemp, hb.sendParams, hb.recvParams, cnts, displs);

    BuffToHalo<float>(bd.oldTemp, hb.recvTemp, dbd.getExtSize());
    BuffToHalo<float>(bd.domParams, hb.recvParams, dbd.getExtSize());
    BuffToHalo<int>(bd.domMap, hb.recvMap, dbd.getExtSize());

    // Asynchronous requests

    MPI_Request *req = new MPI_Request[bd.neighbors->size()+1];
    MPI_Status * stat = new MPI_Status[bd.neighbors->size()+1];

    float  dummyF[1];

    // main simulatilson loop
    bool once = true;
    
    for(unsigned iter = 0; iter < parameters.nIterations; iter++){


        if(parameters.balance && pm.periodElapsed()){
            if(DBG && rank == 0) cout << "detecting" << endl;
  
            pm.balStart();

            if( dbd.loadBalance(pm, bd) ){

      
                if(DBG && rank == 0){
                  cout << "balancing" << endl;
                }

                hb.resize(2*dbd.getHaloLen());

                HaloToBuff<float>(bd.oldTemp, hb.sendTemp, dbd.getExtSize());
                HaloToBuff<float>(bd.domParams, hb.sendParams, dbd.getExtSize());
                HaloToBuff<int>(bd.domMap, hb.sendMap, dbd.getExtSize());

                delete[] cnts;
                delete[] displs;
                delete[] req;
                delete[] stat;
            
                cnts = new int[bd.counts->size()];
                displs = new int[bd.displs->size()];
                req = new MPI_Request[bd.neighbors->size() + 1];
                stat = new MPI_Status[bd.neighbors->size() + 1];

                for(unsigned i = 0; i < bd.counts->size();i++){
                    cnts[i] = bd.counts->at(i);
                    displs[i] = bd.displs->at(i);
                }
            
                ExchangeHaloZones(bd, dbd, hb.sendTemp, hb.recvTemp, hb.sendParams, hb.recvParams, cnts, displs);
            
                BuffToHalo<float>(bd.oldTemp, hb.recvTemp, dbd.getExtSize());
                BuffToHalo<float>(bd.domParams, hb.recvParams, dbd.getExtSize());
                BuffToHalo<int>(bd.domMap, hb.recvMap, dbd.getExtSize());

                MPI_assert(MPI_Barrier(MPI_COMM_WORLD));

            }
            pm.balStop();

            pm.reset();
        } //balancing end

        pm.iterStart(); //timestamp

        pm.imbalDelay(bd.middle, iter, parameters.nIterations, parameters.multiply);


        // compute halo zones
        ComputeHalo(bd, dbd, parameters.airFlowRate, materialProperties.coolerTemp);

        // init communications
        HaloToBuff<float>(bd.newTemp, hb.sendTemp, dbd.getExtSize());

        MPI_assert( MPI_Iscatterv(hb.sendTemp, cnts, displs, MPI_FLOAT, &dummyF, 0, MPI_FLOAT, bd.myCommRank, bd.myComm, &(req[0])),
                     "Temp scatter send failed" LOCATION );

        int idx = 1;
        for(auto n : *bd.neighbors){

            MPI_assert( MPI_Iscatterv(NULL, bd.nData->at(n).scatterCnts, bd.nData->at(n).scatterDispls,
                                     MPI_FLOAT,
                                     &(hb.recvTemp[bd.nData->at(n).displ]), bd.nData->at(n).count, MPI_FLOAT, 
                                     bd.nData->at(n).root, bd.nData->at(n).comm, &(req[idx])),
                    "Temp scatter receive failed" LOCATION );
            idx++;
        }

        if(DBG && once){
            cout << rank << " " << bd;
            once = false;
        }
        // compute the rest 
        unsigned top, right, bottom, left;

        top = bd.topF ? bd.top : bd.top + 2;
        right = bd.rightF ? bd.right : bd.right-2;
        bottom = bd.bottomF ? bd.bottom : bd.bottom - 2;
        left = bd.leftF ? bd.left : bd.left + 2;

        for(unsigned i = top; i < bottom ;i++){
            for(unsigned j = left; j < right ; j++){
            
                ComputePoint(bd.oldTemp,
                             bd.newTemp,
                             bd.domParams,
                             bd.domMap,
                             i, j,
                             dbd.getExtSize().x, 
                             parameters.airFlowRate,
                             materialProperties.coolerTemp
                             );
            }
        }

        // middle column output
        if ((iter % (parameters.nIterations / 10l)) == ((parameters.nIterations / 10l) - 1l))
        {

            if(bd.middle){

                float localSum;

                localSum = dbd.middleColAvg();

                if(bd.midRank == 0){

                    MPI_assert(MPI_Reduce(&localSum, &middleColAvgTemp, 1, MPI_FLOAT, MPI_SUM, 0, bd.COMM_MIDDLE),
                                "middleCol reduce failed" LOCATION
                                );

                    middleColAvgTemp /= materialProperties.edgeSize;
                    // cout << middleColAvgTemp << endl;

                    if(!parameters.batchMode){
                    printf("Progress %ld%% (Average Temperature %.2f degrees)\n", 
                            iter / (parameters.nIterations / 100) + 1, 
                            middleColAvgTemp);
                    }


                }else{

                    MPI_assert(MPI_Reduce(&localSum, NULL, 1, MPI_FLOAT, MPI_SUM, 0, bd.COMM_MIDDLE),
                                "middleCol reduce failed" LOCATION
                                );

                }
            }
        }

        // store to files
        if ( parameters.ioEnabled && (iter % parameters.diskWriteIntensity) == 0){

            stringstream ss;

            if(parameters.useParallelIO){
                ss << "Behavior: invalid HI, parallel " << rank << endl;
                if(file_id == H5I_INVALID_HID) throw runtime_error(ss.str());

                StoreDataIntoFileParallel(file_id,
                            bd.newTemp,
                            materialProperties.edgeSize,
                            dbd.getExtSize().x, dbd.getExtSize().y,
                            dbd.getPosition().x,   //offset in points
                            dbd.getPosition().y,   //offset in points
                            iter / parameters.diskWriteIntensity,
                            iter
                        );

            }else{

                pm.ioStart();

                // all processes
                float * out = dbd.collectData(false);
                // write sequentially
                
                pm.ioEnd();


                if(rank == 0){ // only 0
                    ss << "Behavior: invalid HID, sequential " << rank << endl;

                    if(file_id == H5I_INVALID_HID) throw runtime_error(ss.str());
                    // cout << "serial I/O writing" << endl;
                    
                    StoreDataIntoFile(file_id,
                                    out,
                                    materialProperties.edgeSize,
                                    iter / parameters.diskWriteIntensity,
                                    iter
                                );
                                
                    
                } 
            }


        } // I/O end

        // stop measuring before blcoking call
        pm.iterStop();
        // wait for communications completion
        MPI_assert( MPI_Waitall(bd.neighbors->size()+1, req, stat), "Waitall failed" LOCATION);
        BuffToHalo<float>(bd.newTemp, hb.recvTemp, dbd.getExtSize());

        // swap original pointers inside dbd as well
        dbd.swap(bd.newTemp, bd.oldTemp); 
        
    } //simulation loop

    totalTime = MPI_Wtime() - totalTime;

    if(bd.midRank == 0){

        // [7] Print final result
        if (!parameters.batchMode){
            printf("\nExecution time of parallel version %.5f\n", totalTime);
        }else{
          cout << "Outfile:" <<  parameters.outputFileName.c_str() << endl;
          cout << "Mode:" << (parameters.balance ? "parBal" : "par") << endl;
          cout << "ObjectSize:" << parameters.objDim << endl;
          cout << "MiddleCol:" << middleColAvgTemp << endl;
          cout << "TotalTime:" << totalTime << endl;
          cout << "IterTime:" << totalTime / parameters.nIterations << endl;
          cout << "IterTotal:" << pm.iterTotal << endl;
          cout << "SleepFor:" << pm.sleepfor << endl;
          cout << "SleepTotal[ms]:" << pm.sleepTotal << endl;
          cout << "IOTotal:" << pm.ioTotal << endl;
          cout << "BalanceTotal:" << pm.balTotal << endl;
          cout << "----" << endl;

          }
    } //simulation output end


    delete[] cnts;
    delete[] displs;

    // Asynchronous requests
    delete[] req ;
    delete[] stat;

    return dbd.collectData(false);
}

/**
 * string                     outputFileName)
 * Parallel version of the Heat distribution in heterogenous 2D medium
 * @param [out] parResult          - Final heat distribution
 * @param [in]  materialProperties - Material properties
 * @param [in]  parameters         - parameters of the simulation
 * @param [in]  outputFileName     - Output file name (if NULL string, do not store)
 *
 * @note This is the function that students should implement.                                                  
 * 
 */
 
void ParallelHeatDistribution(float                     *& parResult ,
                              const TMaterialProperties &materialProperties,
                              const TParameters         &parameters,
                              string                     outputFileName)
{
  // Get MPI rank and size
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  hid_t file_id = H5I_INVALID_HID;

  if(!parameters.useParallelIO)
  {
      // Serial I/O
      if(rank == 0 && outputFileName != "")
      {
          if(outputFileName.find(".h5") == string::npos)
              outputFileName.append("_par.h5");
          else
              outputFileName.insert(outputFileName.find_last_of("."), "_par");

          file_id = H5Fcreate(outputFileName.c_str(),
                              H5F_ACC_TRUNC,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
          if(file_id < 0) ios::failure("Cannot create output file");
      }
  }
  else
  {
      // Parallel I/O
      if(outputFileName != "")
      {
          if(outputFileName.find(".h5") == string::npos)
              outputFileName.append("_par.h5");
          else
              outputFileName.insert(outputFileName.find_last_of("."), "_par");

          hid_t hPropList = H5Pcreate(H5P_FILE_ACCESS);
          H5Pset_fapl_mpio(hPropList, MPI_COMM_WORLD, MPI_INFO_NULL);

          file_id = H5Fcreate(outputFileName.c_str(),
                              H5F_ACC_TRUNC,
                              H5P_DEFAULT,
                              hPropList);
          H5Pclose(hPropList);
          if(file_id < 0) ios::failure("Cannot create output file");
      }
  }

  /*
  *   Rank based behavior
  */
  if(rank == 0 && !parameters.batchMode) cout << "Starting parallel simulation..." << endl;

  //allocate and return parallel result
  parResult = Behavior(materialProperties, parameters, file_id, rank, size);
  
  // close the output file
  if (file_id != H5I_INVALID_HID) H5Fclose(file_id);
} // end of ParallelHeatDistribution
//------------------------------------------------------------------------------


/**
 * Store time step into output file (as a new dataset in Pixie format
 * @param [in] h5fileID   - handle to the output file
 * @param [in] data       - data to write
 * @param [in] edgeSize   - size of the domain
 * @param [in] snapshotId - snapshot id
 * @param [in] iteration  - id of iteration);
 */
void StoreDataIntoFile(hid_t         h5fileId,
                       const float  *data,
                       const size_t  edgeSize,
                       const size_t  snapshotId,
                       const size_t  iteration)
{
    hid_t   dataset_id, dataspace_id, group_id, attribute_id;
    hsize_t dims[2] = {edgeSize, edgeSize};

    string groupName = "Timestep_" + to_string((unsigned long long) snapshotId);

    // Create a group named "/Timestep_snapshotId" in the file.
    group_id = H5Gcreate(h5fileId,
                       groupName.c_str(),
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    // Create the data space. (2D matrix)
    dataspace_id = H5Screate_simple(2, dims, NULL);

    // create a dataset for temperature and write data
    string datasetName = "Temperature";
    dataset_id = H5Dcreate(group_id,
                         datasetName.c_str(),
                         H5T_NATIVE_FLOAT,
                         dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id,
           H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           data);

    // close dataset
    H5Sclose(dataspace_id);


    // write attribute
    string atributeName="Time";
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate2(group_id, atributeName.c_str(),
                            H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT);

    double snapshotTime = double(iteration);
    H5Awrite(attribute_id, H5T_IEEE_F64LE, &snapshotTime);
    H5Aclose(attribute_id);


    // Close the dataspace.
    H5Sclose(dataspace_id);

    // Close to the dataset.     H5Dclose(dataset_id);

    // Close the group.
    H5Gclose(group_id);
} // end of StoreDataIntoFile
//------------------------------------------------------------------------------

/**
 * Store time step into output file using parallel version of HDF5
 * @param [in] h5fileId   - handle to the output file
 * @param [in] data       - data to write
 * @param [in] edgeSize   - size of the domain
 * @param [in] tileWidth  - width of the tile
 * @param [in] tileHeight - height of the tile
 * @param [in] tilePosX   - position of the tile in the grid (X-dir)
 * @param [in] tilePosY   - position of the tile in the grid (Y-dir)
 * @param [in] snapshotId - snapshot id
 * @param [in] iteration  - id of iteration
 */
void StoreDataIntoFileParallel(hid_t h5fileId,
                               const float *data,
                               const size_t edgeSize,
                               const size_t tileWidth, const size_t tileHeight,
                               const size_t tilePosX, const size_t tilePosY,
                               const size_t snapshotId,
                               const size_t iteration)
{
    hid_t dataset_id, dataspace_id, group_id, attribute_id, memspace_id;
    const hsize_t dims[2] = { edgeSize, edgeSize };
    const hsize_t offset[2] = { tilePosY, tilePosX };
    const hsize_t tile_dims[2] = { tileHeight, tileWidth };
    const hsize_t core_dims[2] = { tileHeight - 4, tileWidth - 4 };
    const hsize_t core_offset[2] = { 2, 2 };

    string groupName = "Timestep_" + to_string((unsigned long)snapshotId);

    // Create a group named "/Timestep_snapshotId" in the file.
    group_id = H5Gcreate(h5fileId,
                         groupName.c_str(),
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create the data space in the output file. (2D matrix)
    dataspace_id = H5Screate_simple(2, dims, NULL);

    // create a dataset for temperature and write data
    string datasetName = "Temperature";
    dataset_id = H5Dcreate(group_id,
                           datasetName.c_str(),
                           H5T_NATIVE_FLOAT,
                           dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // create the data space in memory representing local tile. (2D matrix)
    memspace_id = H5Screate_simple(2, tile_dims, NULL);

    // select appropriate block of the local tile. (without halo zones)
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, core_offset, NULL, core_dims, NULL);

    // select appropriate block of the output file, where local tile will be placed.
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, core_dims, NULL);

    // setup collective write using MPI parallel I/O
    hid_t hPropList = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(hPropList, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, hPropList, data);

    // close memory spaces and property list
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Pclose(hPropList);

    // write attribute
    string attributeName = "Time";
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate2(group_id, attributeName.c_str(),
                              H5T_IEEE_F64LE, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT);

    double snapshotTime = double(iteration);
    H5Awrite(attribute_id, H5T_IEEE_F64LE, &snapshotTime);
    H5Aclose(attribute_id);

    // close the dataspace
    H5Sclose(dataspace_id);

    // close the dataset and the group
    H5Dclose(dataset_id);
    H5Gclose(group_id);
}
//------------------------------------------------------------------------------

/**
 * @brief Halo zones template function test
 * @details [long description]
 * 
 * @param rank [description]
 */
void BorderTest(int rank)
{

    Dims ext(36,36);

    const int bsize = 36*36;
    const int halosize = 256;  
    int * array = new int[bsize];
    int * halo = new  int[halosize];
    // int * test = new  int[halosize];

    //fill with values

    for(unsigned i = 0; i < bsize;i++){
        array[i] = 0;
    }
    

    for(unsigned i = 2; i < 4;i++)  {
        for(unsigned j = 2; j < ext.x-2;j++){

          if(i == 2){
            if(j < ext.x / 2)
              array[i * ext.x + j] = 1;
            else
              array[i * ext.x + j] = 2;
          }else{

            if(j < ext.x / 2)
              array[i * ext.x + j] = 3;
            else
              array[i * ext.x + j] = 4;
          }


      }
    }
    // right
    for(unsigned i = ext.x - 4;i < ext.x -2;i++){
      for(unsigned j = 2; j < ext.y-2;j++){

        if(i == ext.x - 4){
          if(j < ext.y / 2)
            array[j * ext.x + i] = 5;
          else
            array[j * ext.x + i] = 6;
        }else{

          if(j < ext.y / 2)
            array[j * ext.x + i] = 7;
          else
            array[j * ext.x + i] = 8;
        }
      }
    }

    //bottom
    for(unsigned i = ext.y - 4;i < ext.y - 2 ;i++){
      for(unsigned j = 2; j < ext.x-2;j++){

        if(i == ext.y - 4){
          if(j < ext.x / 2)
            array[i * ext.x + j] = 9;
          else
            array[i * ext.x + j] = 10;

        }else{
            if(j < ext.x / 2)
            array[i * ext.x + j] = 11;
          else
            array[i * ext.x + j] = 12;
        }
      }
    }

    //left
    for(unsigned i = 2 ;i < 4 ;i++){
      for(unsigned j = 2; j < ext.y-2;j++){
        if(i == 2){
          if(j < ext.y / 2)
            array[j * ext.x + i] = 13;
          else
            array[j * ext.x + i] = 14;

        }else{

          if(j < ext.y / 2)
            array[j * ext.x + i] = 15;
          else
            array[j * ext.x + i] = 16;

        }
      }
    }

    PrintBlock(array, Dims(36,36), rank, true);

    HaloToBuff(array, halo, Dims(36,36) );
    
    BuffToHalo(array, halo, Dims(36,36));

    cout << "Validate: " << endl;
    PrintBlock(array, Dims(36,36), rank, true);


}

/**
 * Main function of the project
 * @param [in] argc
 * @param [in] argv
 * @return
 */
int main(int argc, char *argv[])
{
    int rank, size;




    // Initialize MPI
    MPI_Init(&argc, &argv);

    ParseCommandline(argc, argv, parameters);

    // intiialize zoltan library
    float version;
    Zoltan_Initialize(argc, argv, &version);

    // Get MPI rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    TLogger::SetLevel(TLogger::BASIC);



    if (rank == 0)
    {
        // Create material properties and load from file
        materialProperties.LoadMaterialData(parameters.materialFileName, true);
        parameters.edgeSize = materialProperties.edgeSize;

        parameters.PrintParameters();
    }
    else
    {
        // Create material properties and load from file
        materialProperties.LoadMaterialData(parameters.materialFileName, false);
        parameters.edgeSize = materialProperties.edgeSize;
    }

    if (parameters.edgeSize % size)
    {
        if (rank == 0)
            printf("ERROR: number of MPI processes is not a divisor of N\n");

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (parameters.IsRunSequntial())
    {
        if (rank == 0)
        {
            // Memory allocation for output matrices.
            // seqResult = (float*)_mm_malloc(materialProperties.nGridPoints * sizeof(float), DATA_ALIGNMENT);
            seqResult = new float[materialProperties.nGridPoints];

            SequentialHeatDistribution(seqResult,
                                       materialProperties,
                                       parameters,
                                       parameters.outputFileName);
        }
    }

    if (parameters.IsRunParallel())
    {
        // Memory allocation for output matrix.

        ParallelHeatDistribution(parResult,
                                 materialProperties,
                                 parameters,
                                 parameters.outputFileName);

        // for(int i = 0; i < materialProperties.nGridPoints;i++){
            // cout << parResult[i] << " " ;
        // }
    
        TLogger::Log(TLogger::FULL, OUT_FMT_END_OF_SIMULATION);
    }

    // TLogger::Log(TLogger::FULL, OUT_FMT_END_OF_SIMULATION);

    // Validate the outputs
    if (parameters.IsValidation() && rank == 0)
    {
        if (parameters.debugFlag)
        {
          printf("---------------- Sequential results ---------------\n");
          PrintArray(seqResult, materialProperties.edgeSize);
    
          printf("----------------- Parallel results ----------------\n");
          PrintArray(parResult, materialProperties.edgeSize);
        }

        // PrintArray(seqResult, parameters.edgeSize);
        // printf("%s\n", "======================================" );
        // PrintArray(parResult, parameters.edgeSize);


        if (VerifyResults(seqResult, parResult, parameters, 0.001f))
          printf("Verification OK\n");
        else
          printf("Verification FAILED\n");
    } 

    /* Memory deallocation*/
    delete[] seqResult;
    delete[] parResult;

    MPI_Finalize();

    return EXIT_SUCCESS;
} // end of main
//------------------------------------------------------------------------------
