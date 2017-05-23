/***********************************************
*
*  File Name:       DynamicDynamicBlockDescriptor.h
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Core class for dynamic load balancing
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/

#ifndef __DYNAMIC_BLOCK_DESCRIPTOR_H__
#define __DYNAMIC_BLOCK_DESCRIPTOR_H__

#include <mpi.h>
#include <zoltan_cpp.h>


#include <string>
#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <iterator>

#include <new>
#include <stdexcept>
#include <immintrin.h>
#include <algorithm>


#include <MaterialProperties.h>
#include <BasicRoutines.h>
#include <Logger.h>
#include <Dims.h>
#include <TileMsg.h>
#include <Neighbor.h>
#include <TopologyDescriptor.h>
#include <LoadBalancer.h>
#include <PerfMeasure.h>
#include <BlockData.h>
#include <Asserts.h>



/*  !! DO NOT CHANGE !! */

namespace DLB {

const int NO_RANK = -1;   // constant for invalid rank, DO NOT CHANGE!

const int ZLT_ID_TOO_LONG = 1;
const int ZLT_UNSUP_DIMS = 2;

const int HALO_SIZE = 2;

/**
* @brief Maps given coordinates to array index respection halo zone size
* 
*/ 

inline unsigned halo(unsigned x, unsigned y, unsigned esize)
{
    return ((y + HALO_SIZE) * esize) + (x + HALO_SIZE);
}
              // CheckRank method relies on it

// using namespace std;




/**
 * @brief Class describing distributed topology.
 * @details Provides useful methods for topology related prolems.
 * Computes all necessary information, like block size and neighborhood
 * based MPI rank and world size.
 * 

 */

class DynamicBlockDescriptor{

public:


  /**
   * Describe appropriate side of 2D matrix
   */

    typedef enum pos {

        TOP = 0, 
        RIGHT, 
        BOTTOM, 
        LEFT

    } BlockPosition;


    /**
     * @brief DynamicBlockDescritor constructor
     *         
     * @details Constrcuts DBD and all attribute objects.
     * 
     * @param rank
     * @param edgeSize [in] size of 2D domain matrix edge 
     * @param rank [in] current process MPI rank
     * @param worldSize [in] size of MPI world
     */
    DynamicBlockDescriptor( int rank, int worldSize, size_t edgeSize, Dims objSize, double threshold);

    ~DynamicBlockDescriptor(void);


    /**
     * @brief Initial static load balance
     * @details Initialize data etc. 
     *          Must be called before loadBalance() by all
     *          processes.
     *          
     * @param props - initial model data from file
     *                significant only for root
                    
     * @return updated BlockData

     */

    BlockData loadInit(const TMaterialProperties & props);


    /**
     * @brief Dynamic load balance method
     * 
     * @details This method needs measured performance, which is send to master
     * and evaluated. Re-partition and data migration may occur within this call
     * based on measured values.
     *
     * Upon true returned, balancing performed, BlockData will be updated
     * Re-initialization of neighbor comms is necessary.
     * 
     * @param pm - performance measurement object with data about last iterations
     * @param block - BlockData reference, will be update if load balance occurs
     * @return true if balancing performed
     * 
     */

    bool loadBalance(PerfMeasure & pm, BlockData & block);

    /**
     * @brief   Compares old and new decomposition and setes appropriate
     *          migration arrays - import/export GIDs etc.
     *          
     * @details Computes sets of import and export objects.
     * 
     * @param oldTiles - old decomposition
     * @param newTiles - new decoposition
     * 
     * @return persistent GIDs (these will remain on currnet process)
     */

    list<unsigned> * resolveMigration(  const vector<TileDescriptor> & oldTiles,
                            const vector<TileDescriptor> & newTiles
                        );

    /**
     * @brief Data migration function
     * 
     * @details Initializes memory for assigned block and
     *          calls migrate function
     * 
     * @param persist - remaining GIDs
     */

   void migrate(const vector<TileDescriptor> & newTiles,
                            list<unsigned> & persist
                    );

   /**
    * @brief Collect data to master process for serial I/O purposes.
    * 
    * @details Collects Zoltan objects from all processes to master
    * 
    * @param old - if true, data are picked from oldTemp array
    * @return coollected data pointer, this can be passed to HDF5
    */

   float * collectData(bool old = true);

    /**
     * @brief Initialize Zoltan related params etc.
     * 
     * @details This may be part of constructor, but It was
     *          designed as part of API, when underlying library changes.
     */

    void zoltanInit(void);


    /**
     * @brief Custom pointer swapping
     * @details Swaps args as well as pointer inside dbd
     * 
     * @param p1 pointer
     * @param p2 pointer
     */

    template<class T>
    void swap(T& p1, T& p2)
    {
        // swap arguments
        T c( std::move(p1) ); 
        p1=std::move(p2);
        p2=std::move(c);
    
        // swap pointers in class for migration purposes
        std::swap(bdata.oldTemp, bdata.newTemp);
    }

    bool isTopChanged(void) {return tdesc.getTopChanged(); }


    /**
     * @brief Metadata update
     * 
     * @details [long description]
     */

    void blockUpdate(void);


    // compute and return middle column temp
    float middleColAvg(void);

    /**
     * @brief Calculates MPI ranks, which belong to neighbor blocks.
     * @return vector reference
     */
    const vector<int> & getAdjRanks(void);


    /**
     * @brief Creates instances of MPI datatypes describing block
     * @return non-zero on failure
     */

    int initDatatypes(int blockPadding = 2, int resizeValue = -1);

    /**
     * @brief Computes displacements vector for scattering
     */
    const    vector<int> & getDisplsVect(void);
    /**
     * @brief Computes displacement vector for scattering
     * @return int array pointer, which may be passed to scatter
     */

    int * getDisplsArr(void);

    /**
    * @brief Returns offset of block by given process rank
    * 
    */

    int * getBlockOffset(void);

    /**
    * @brief Return position of actual block or block belonging to rank.
    * @details Position is in range <0, blocks_in_row>, <0, blocks_in_column>
    * 
    * @param anyrank - any valid MPI rank
    * @return Dims object 
    */

    Dims getPosition(void)    {return tdesc.tile().getPosition(); }    //return block position


    /**
    * @brief Print given array of values in matrix format
    * @details block must be initialized pointer with size at least 
    * of blockArea
    * 
    * @param block - array to print
    * @param extended - include or exclude halo borders
    */


    /**
     * @brief Returns extended area of currently assigned lock
     * @return [description]
     */

    unsigned getExtArea(void)                { return tdesc.tile().getExtArea(); }
    /**
     * @brief Returns area fo whole computational domain
     */

    unsigned getDomainArea(void)             { return tdesc.getDomainArea(); }
    /**
     * @brief Returns area of current block
     */

    unsigned getBlockArea(void)              { return tdesc.tile().getArea(); }

    /**
     * @brief Returns stored MPI rank
     */
    int getRank(void)  const            { return rank; }

    int getRank(const Dims & blockPosition);

    unsigned getHaloLen(void);

    /**
     * @brief Returns dimensions of current block
     * @details Dims object
     */

    Dims getBlockSize()       { return tdesc.tile().getSize(); }
    /**
     * @brief Returs size of block extended by halo borders
     * @return Dims object
     */

    Dims getExtSize()         { return tdesc.tile().getExtSize(); }
    /**
     * @brief Position of current block in range <0, blocks_in_row>, <0, blocks_in_col>
     * @return [description]
     */
    // Dims getOffset()    const                         { return offset; }
    /**
     * @return stored MPI world size
     */
    int getWorldSize() const       { return worldSize; }

    // int * getBounds(void);

    /**
     * @brief Returns BlockData pointer
     * @return valid pointer or throws runtime_error
     */
    BlockData    getBlockData(void);

    void initBlockData(const TMaterialProperties & data);
    
    /**

     * @brief Print object metadata for debug purpose
     */

    string toString(void);
    string tilesToString(void);

    // vector<unsigned> * getAssignedObjs(void);

    int rank, worldSize;


protected:


    // balancing counter
    unsigned balanceSeq;
    
    size_t edgeSize;

    bool datatypesInitialized;


    // describes actual domain partitioning
    TopologyDescriptor tdesc;

    // load balancing object
    LoadBalancer    lb;

    // all necessary data for actual block
    BlockData bdata;

    // temporary data for data migration purpose
    TempBlock newBlock;


    // zoltan object size
    Dims objectSize;

    // object cols & rows
    unsigned rows;
    unsigned cols;

    //assigned - to this particular block at the moment
    // total - all object in the model domain
    unsigned assignedObjsCnt, totalObjsCnt;

    
    /**
     *  @details Collecting data from all processes to root.
     *  If collectDataFlag = true, zoltan unpack callbacks, unpacks
     *  data to resArray.
     *  
     *  resArray must be previously allocated to
     *  (nGridPoints + halo) * sizof(float)
     */

    bool collectDataFlag;
    bool oldArray;
    float  * resArray;

    /**
     * @brief Mapping from Zoltan object position to object GID
     *
     * @details Computed on regular mesh basis
     * 
     * @param objPos [description]
     * @return GID
     */

    unsigned getGIDByCoords(const Dims & objPos);

    /**
     * @brief Inversion of getGIDByCoords()
     * 
     * @details Computer position from object GID
     * 
     * @param gid  global ID
     * @return object position
     */
    Dims getCoordsByGID(unsigned gid);

    /**
     * @brief GIDs assigned to current block
     * @return GID list
     */

    list<unsigned> * getAssignGIDs(void);

    /**
     * @brief GIDs assigned to block t
     * @return GID list
     */
    list<unsigned> * getAssignGIDs(const TileDescriptor & t);

    /**
    * 
    * @brief       Sorting object to rank assignement
    * 
    * @detailes    Sorting objects to be assigned to ranks
    *              according to mesh numbering beginning at upper-left
    *              corner to bottom-right (0, N).
    *              Useful when external parititioner used, which may choose
    *              his own rules.
    * 
    * @param list [description]
    * @return [description]
    */

    void sortGIDs(unsigned * list, unsigned objsPerBlock);


    /**
    *  @brief Checks if assigned objects makes together valid block
    *         and compute position and size based on assigned objects.
    *         
    *  @return  TD with update params based on assigned objs
    */

    TileDescriptor resolvePartition(unsigned * gids, unsigned size);

    /**
     * @brief Computes assigned objects to all tiles defined by param
     *
     * @details Object lists are computed by given block position and size
     * 
     * @param tds - vector of tiles
     * @return std::map mapping rank to object list
     */

    map<int, list<unsigned>*>  getAssignedObjs(const vector<TileDescriptor> & tds);


    /**
     * @brief Check new partition validity and create relevant
     *        TileDescriptor
     *                
     * @details
     * 
     * @param assigned    - list of assigned objects
     * @return TileDescriptor describing new block
     */

    TileDescriptor resolvePartition(list<unsigned> & assigned);



    void initNewBlock(TileDescriptor t);

    /**
     * @brief Move objects from actual block to new block arrays
     * @details [long description]
     * 
     * @param persist - objects assigned to actual rank
     * @param init - first time migration  (begin of simulation)
     * @param collect - collecting results from ranks to root
     */

    void movePersistObj(const list<unsigned> & persist, bool init, bool collect = false);



/*===============================

    ZOLTAN CALLBACKS   

  ===============================
*/

    // if true, zoltan callbalcks produces debug output
    bool callbackDbg;

    static int zolt_num_obj_fn(void * data, int * err);
    static int zolt_num_geom_fn(void * data, int *err);

    static void zolt_geom_multi_fn( void * data, 
                              int num_gid_entries, 
                              int num_lid_entries, 
                              int num_obj,
                              ZOLTAN_ID_PTR global_ids, 
                              ZOLTAN_ID_PTR local_ids,
                              int num_dim,
                              double *geom_vec,
                              int * ierr
                            );

    static void zolt_geom_fn( void * data, 
                              int num_gid_entries, 
                              int num_lid_entries, 
                              ZOLTAN_ID_PTR global_id, 
                              ZOLTAN_ID_PTR local_id,
                              double *geom_vec,
                              int * ierr
                            );

    static int zolt_obj_size_fn(void * data,
                     int num_gid_entries,
                     int num_lid_entries,
                     ZOLTAN_ID_PTR global_id,
                     ZOLTAN_ID_PTR local_id,
                     int *ierr
                     );

    static void zolt_pack_obj_fn( void * data,
                            int num_gid_entries,
                            int num_lid_entries,
                            ZOLTAN_ID_PTR global_id,
                            ZOLTAN_ID_PTR local_id,
                            int dest, 
                            int size,
                            char * buf,
                            int *ierr
                            );

    static void zolt_unpack_obj_fn( void * data,
                                  int num_gid_entries,
                                  ZOLTAN_ID_PTR global_id,
                                  int size,
                                  char *buf,
                                  int *ierr
                                  );


    static void zolt_obj_list_fn( void *data,
                                int num_gid,
                                int num_lid, 
                                ZOLTAN_ID_PTR global_ids,
                                ZOLTAN_ID_PTR local_ids,
                                int wgt_dim,
                                float * obj_wgts,
                                int *ierr);
    
    /**
     * @brief Computes current block dimensions
     */

    void computeBlockSize(void);      
    /**
     * @brief Computes neighbor block ranks
     */
    void computeAdjacentRanks(void);    
    /**
     * @brief Computes MPI scatterv data displacements
     */
    void computeDisplacements(void);

    /**
     * @brief Inavlid rank assertion - throws runtime error
     * 
     * @details Checks if computed rank is valid in actual MPI world.
     * Throws runtime_error on inavlid rank
     *  
     * (May perform other checks in future)
     * 
     * @param rank - rank to be checked
     * 
     */

    void checkRank(int rank);



};

} //DLB nspace end

#endif 

