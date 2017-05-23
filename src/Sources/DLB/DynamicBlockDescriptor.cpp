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


#include "DynamicBlockDescriptor.h"

using DBD = DLB::DynamicBlockDescriptor;
using std::pair;
using std::list;

DBD::DynamicBlockDescriptor(int rank, int worldSize,  size_t edgeSize, Dims objSize, double threshold):
rank(rank),
worldSize(worldSize),
edgeSize(edgeSize),
tdesc(rank, worldSize, edgeSize),
lb(rank, edgeSize, worldSize, objSize, threshold),
objectSize(objSize)
{

    // bdata = tdesc.getBlockData();
   
    // regular matrix of zoltan objects
    cols = edgeSize / objSize.x;
    rows = edgeSize / objSize.y;

    if(edgeSize % objectSize.x != 0)
        throw runtime_error("DBD : edgeSize % objSize != 0 ");

    totalObjsCnt = (edgeSize / objectSize.x) * (edgeSize / objectSize.y);

    if(rank == 0){
        assignedObjsCnt = totalObjsCnt; 
    }else{
        assignedObjsCnt = 0;
    }

    callbackDbg = false;
    oldArray = true;

    collectDataFlag = false;

    if(rank == 0){
        resArray = new float[edgeSize * edgeSize];

        for(unsigned i = 0; i < edgeSize * edgeSize;i++)
            resArray[i] = 0.0;

    }

    balanceSeq = 0;

}

DBD::~DynamicBlockDescriptor(void)
{

  // oldTemp is returned as result
  // do not dealocate
  delete[] bdata.newTemp;
  delete[] bdata.domMap;
  delete[] bdata.domParams;

}



void DBD::zoltanInit(void)
{
    /**
     * Init zoltan object in LoadBalancer class
     */

    // communicator for Zoltan, duplicate of COMM_WORLD is created
    // as precausion for possible clashes

    MPI_assert(MPI_Comm_dup(MPI_COMM_WORLD, &(lb.zoltComm)), "zoltComm duplication failed" LOCATION);

    // create zoltan object
    lb.zz =  new Zoltan(lb.zoltComm);

    /*
    * Set params and callbacks
    */

    // lb.zz->Set_Param("LB_METHOD", "RIB");
    // lb.zz->Set_Param("RCB_RECTILINEAR_BLOCKS", "1");

    // lb.zz->Set_Param("RCB_REUSE", "1");
   
    // lb.zz->Set_Param("AVERAGE_CUTS", "1");
   
    // lb.zz->Set_Param("OBJ_WEIGHTS_COMPARABLE", "1");

    // GID - global object ID, LIDs are unused
    // number of UINT describing entry
    lb.zz->Set_Param("NUM_GID_ENTRIES", "1");
    lb.zz->Set_Param("NUM_LID_ENTRIES", "1");

    // remember domain decomposition
    // lb.zz->Set_Param("KEEP_CUTS", "1");
    // number of part = processes
    stringstream ss;
    ss << worldSize;
    lb.zz->Set_Param("NUM_GLOBAL_PARTS", ss.str() );
    lb.zz->Set_Param("NUM_LOCAL_PARTS", "1");

    // LB_partition return import, export lists
    lb.zz->Set_Param("RETURN_LISTS", "ALL");

    // amount of output
    lb.zz->Set_Param("DEBUG_LEVEL", "0");
    // if DEBUG_LEVEL > 5 (routine trace info) is printed by root
    lb.zz->Set_Param("DEBUG_PROCESSOR", "0");
    // debug Zoltan memory management
    // lb.zz->Set_Param("DEBUG_MEMORY", "3");

    //function callbacks
    // static private methods, pass "this" as user data

    lb.zz->Set_Fn(  ZOLTAN_NUM_OBJ_FN_TYPE,       (void (*)()) zolt_num_obj_fn,     this);
    lb.zz->Set_Fn(  ZOLTAN_GEOM_MULTI_FN_TYPE,    (void (*)()) zolt_geom_multi_fn,  this);
    lb.zz->Set_Fn(  ZOLTAN_GEOM_FN_TYPE,          (void (*)()) zolt_geom_fn,        this);
    lb.zz->Set_Fn(  ZOLTAN_NUM_GEOM_FN_TYPE,      (void (*)()) zolt_num_geom_fn,    this);
    lb.zz->Set_Fn(  ZOLTAN_OBJ_SIZE_FN_TYPE,      (void (*)()) zolt_obj_size_fn,    this);
    lb.zz->Set_Fn(  ZOLTAN_OBJ_LIST_FN_TYPE,      (void (*)()) zolt_obj_list_fn,    this);

    //migration callbacks
    lb.zz->Set_Fn(  ZOLTAN_PACK_OBJ_FN_TYPE,      (void (*)()) zolt_pack_obj_fn,    this);
    lb.zz->Set_Fn(  ZOLTAN_UNPACK_OBJ_FN_TYPE,    (void (*)()) zolt_unpack_obj_fn,  this);

}





/**
 * @brief Allocate necessary temporary arrays before data migration
 * @details [long description]
 * 
 * @param t - describing new tile geometry
 */
void DBD::initNewBlock(TileDescriptor t)
{

    newBlock.tile = t;
    newBlock.temp = new float[newBlock.tile.getExtArea()];
    newBlock.params = new float[newBlock.tile.getExtArea()];
    newBlock.map = new int[newBlock.tile.getExtArea()];

    for(unsigned i = 0; i < newBlock.tile.getExtArea(); i++ ){

        newBlock.temp[i] = 0.0;
        newBlock.params[i] = 0.0;
        newBlock.map[i] = 0;

    }

}

void DBD::initBlockData(const TMaterialProperties & data)
{
    unsigned earea = ( data.edgeSize + 4) * (data.edgeSize + 4);

    bdata.newTemp = new float[earea];
    bdata.oldTemp = new float[earea];
    bdata.domParams = new float[earea];
    bdata.domMap = new int[earea];

    Dims d(data.edgeSize, data.edgeSize);
    unsigned esize = data.edgeSize + 4;

    for(unsigned i = 0; i < d.y;i++){
        for(unsigned j = 0; j < d.x;j++){

            bdata.newTemp[halo(j,i,esize)] = data.initTemp[i*d.x + j];
            bdata.oldTemp[halo(j,i,esize)] = data.initTemp[i*d.x + j];
            bdata.domParams[halo(j,i,esize)] = data.domainParams[i*d.x + j];
            bdata.domMap[halo(j,i,esize)] = data.domainMap[i*d.x + j];
        }
    }

}



DLB::BlockData DBD::loadInit(const TMaterialProperties & props )
{

    if(lb.zz == NULL) throw runtime_error("loadInit: Zoltan not initialized");


    TileMsg * buf = new TileMsg[worldSize];

    unsigned objsPerBlock = totalObjsCnt / worldSize;

    unsigned * rootObjs = new unsigned[objsPerBlock];

    // compute staic decomposition - regular mesh
    vector<TileDescriptor> * vtd = lb.regularTiles();

    map<int, list<unsigned> *> gids = getAssignedObjs(*vtd);

    stringstream ss;

    if(rank == 0){

        // persist = getAssignGIDs(vtd->at(0));
        if(DBG){
            cout << COUTLOC;
            ss << vtd->at(rank);
            ss << "assigned: ";
            for(auto x : *(gids.at(rank)) )
                ss << x << " ";
            ss << endl;
    
            synCout(ss.str(),rank, worldSize);
        }

        // set actual tile
        tdesc.setTile(TileDescriptor(rank, 0,0, edgeSize, edgeSize));
        initNewBlock(vtd->at(0)); 

        initBlockData(props); //laod initial data

        lb.num_import = 0;
        lb.import_procs = NULL;
        lb.import_to_part = NULL;
        lb.import_global_ids = NULL;
        lb.import_local_ids = NULL;

        lb.num_export = totalObjsCnt - gids.at(0)->size();
        // cout << "r0 export" << lb.num_export << endl;
        lb.export_procs = new int[lb.num_export];
        lb.export_global_ids = new unsigned[lb.num_export];
        lb.export_local_ids = new unsigned[lb.num_export];
        lb.export_to_part = new int[lb.num_export];

        unsigned offset = 0;
        for(auto const& x : gids){
            auto it = x.second->begin();

            if(x.first == 0) //skip rank 0
                continue;

            for(unsigned i = 0; i < x.second->size();i++){
                lb.export_procs[offset + i] = x.first;
                lb.export_global_ids[offset + i] = *it;
                lb.export_local_ids[offset + i] = *it;
                lb.export_to_part[offset + i] = 0;
                advance(it, 1);
            }   
            offset += x.second->size();

        }


        // migrate data by Zoltan
        lb.migrateData();


        // move objects, which are not migrated
        // bdata will have actual data
        movePersistObj(*(gids.at(0)), false);

        tdesc.setTile(vtd->at(0));


        // new tile message
        TileMsg msg(vtd->at(0).getPosition(), vtd->at(0).getSize(), rank, tdesc.tile().getHostNumber());

        // receive tile info, rank<->hostNumber mapping at this point

        MPI_assert( MPI_Gather(&msg, 1, tdesc.TileMsg_t, buf, 1, tdesc.TileMsg_t, 0, MPI_COMM_WORLD ), "Gather failed" LOCATION );

        MPI_assert( MPI_Bcast(buf, worldSize, tdesc.TileMsg_t, 0, MPI_COMM_WORLD), "Bcast failed" LOCATION );
        // bcast new topology to all
        // udpate topol. related data

        tdesc.setTiles(buf);
        tdesc.updateTopology();

    }else{

        if(DBG){
            ss << vtd->at(rank);
            ss << "assigned: ";
            for(auto x : *(gids.at(rank)) )
                ss << x << " ";
            ss << endl;
            synCout(ss.str(),rank, worldSize);
        }

        // create msg with rank and hostNumber only
        // send msg to master
    
        tdesc.setTile(vtd->at(rank));

        initNewBlock(vtd->at(rank));

        lb.num_import = gids.at(rank)->size();

        lb.import_global_ids = new unsigned[lb.num_import];
        lb.import_local_ids = new unsigned[lb.num_import];
        lb.import_procs = new int[lb.num_import];
        lb.import_to_part = new int[lb.num_import];

        lb.num_export = 0;
        lb.export_procs = 0;
        lb.export_global_ids = NULL;
        lb.export_local_ids = NULL;
        lb.export_to_part = NULL;

        auto x = gids.at(rank);
        auto it = x->begin();

        for(unsigned i = 0; i < x->size();i++){

            lb.import_procs[i] = 0;
            lb.import_global_ids[i] = *it;
            lb.import_local_ids[i] = *it;
            lb.import_to_part[i] = 0;
            advance(it, 1);
        }   


        lb.migrateData();

        movePersistObj(*(gids.at(rank)), true);

        TileMsg msg(vtd->at(rank).getPosition(), vtd->at(rank).getSize(), rank, tdesc.tile().getHostNumber());

        // send my tile to master
        MPI_assert( MPI_Gather(&msg, 1, tdesc.TileMsg_t, NULL, 1, tdesc.TileMsg_t, 0, MPI_COMM_WORLD ), "Gather failed" LOCATION);

        // receive new topology from master
        MPI_assert( MPI_Bcast(buf, worldSize, tdesc.TileMsg_t, 0, MPI_COMM_WORLD), "Bcast failed" LOCATION );

        //set complete topology        
        tdesc.setTiles(buf); 
        tdesc.updateTopology();


    }

    lb.clearArrays();

    delete[] buf;
    
    for(auto i : gids)
        delete i.second;
    
    delete vtd;
    delete[] rootObjs;
    
    
    // calls TileDescriptor::getBlockData for topology info
    // then adds data array from DBD
    return getBlockData();


}



void DBD::sortGIDs(unsigned * lst, unsigned objsPerBlock)
{
    if(lst == NULL) throw runtime_error("sortGIDs null passed");

    map<unsigned, vector<unsigned> > tmp;

    vector<unsigned> l;
    list<unsigned> first;

    for(int i = 0; i < worldSize;i++){
        for(unsigned o = 0; o < objsPerBlock; o++ ){
            l.push_back(lst[i*objsPerBlock + o]);
            // map[i].push_back(i*objsPerBlock + o);
        }
        sort(l.begin(), l.end());
        first.push_back(l.front());
        tmp.insert(pair<unsigned, vector<unsigned>>(l.front(), l));
        l.clear();
    }
    first.sort();

    int of = 0;
    for(auto f: first){
        for(unsigned i =0; i < objsPerBlock;i++){
            lst[of * objsPerBlock + i] = tmp[f].at(i);
        }
        of++;
    }

}


bool DBD::loadBalance(PerfMeasure & pm, BlockData & block)
{
    int balancing = 0;
    bool restoreRegular = false;
    double * rbuf = new double[worldSize];

    TileMsg * tbuf = new TileMsg[worldSize];

    list<unsigned> * persist = NULL;
    vector<TileDescriptor> *tls = NULL;

    vector<unsigned> sizes;

    double time = pm.getAgreg();

    stringstream ss;

    // if(rank == 0 || rank == 3)
        // time *= 3;
    if(rank == 0){



        //collect times from other ranks
        MPI_assert( MPI_Gather( &time, 1, MPI_DOUBLE, 
                                rbuf, 1, MPI_DOUBLE,
                                0, MPI_COMM_WORLD
                              ), 
                                "loadBalance: root Gather failed" LOCATION );

       

        // detect imbalance if occurs
        // send balancing info to others
        for(int i = 0; i < worldSize;i++){
            auto it = find(tdesc.getTiles().begin(), tdesc.getTiles().end(), i);
            unsigned tmp = (*it).getSize().x / objectSize.x;
            sizes.push_back(tmp);
        }       

        // move to vector
        cout << "times_" << balanceSeq << ": ";

        vector<double>  times;
        for(int i =0; i < worldSize;i++){
            times.push_back(rbuf[i]);
            cout << rbuf[i] << " ";
        }
        cout << endl;

        balanceSeq++;

        if(lb.isBalanced(times)){

            if(lb.imbalance){ //actualy imbalance, recovering to regular

                restoreRegular = true;
                balancing = 1;
                MPI_assert( MPI_Bcast(&balancing, 1, MPI_INT, 0, MPI_COMM_WORLD), "loadBalance bcast failed" LOCATION  );
            }

        }else{

            balancing = 1;
            restoreRegular = false;
            MPI_assert( MPI_Bcast(&balancing, 1, MPI_INT, 0, MPI_COMM_WORLD), "loadBalance bcast failed" LOCATION );
        }

        if(!balancing){

            // send no balancing
            MPI_assert( MPI_Bcast(&balancing, 1 , MPI_INT, 0, MPI_COMM_WORLD), "loadBalance bcast failed" LOCATION  );

        }else{

            // balancing = 1;
            // lb.imbalance = true;

            // MPI_assert( MPI_Bcast(&balancing, 1, MPI_INT, 0, MPI_COMM_WORLD), "loadBalance bcast failed" LOCATION );

            // obtain new  topology
            if(restoreRegular){
                tls = lb.regularTiles();
                lb.imbalance = false;
            }else{
                tls = lb.getPartition(times, tdesc.getTiles());
                lb.imbalance = true;
            }
                

            if(DBG){
                cout << COUTLOC;
                ss << tls->at(rank);
                ss << "assigned: ";
                list<unsigned> * tmp = getAssignGIDs(tls->at(rank));
                for(auto x : *tmp )
                    ss << x << " ";
                ss << endl;
                cout << ss.str();
                delete tmp;
            }


            for(int i =0; i < worldSize;i++){  // new topology to send buffer
                tbuf[i] = tls->at(i).getData();
            }


            // Bcast new topology to others
            MPI_assert( MPI_Bcast(tbuf, worldSize, tdesc.TileMsg_t, 0, MPI_COMM_WORLD ), "topology Bcast failed" LOCATION);
        
            // set migration data            
            persist = resolveMigration(tdesc.getTiles(), *tls);

            // must be called before updateTopology
            // callbackDbg = true;
            migrate(*tls, *persist);
            // callbackDbg = false;

            tdesc.setTiles(*tls);


            

            tdesc.updateTopology();
            block = getBlockData();

            if(DBG){
                cout << COUTLOC;
                synCout(tdesc.commsToString(), rank, worldSize);
            }

        }

        cout << "sizes_" << balanceSeq << ": ";
        for(int i = 0; i < worldSize;i++){
            auto it = find(tdesc.getTiles().begin(), tdesc.getTiles().end(), i);
            unsigned tmp = (*it).getSize().x / objectSize.x;
            cout << tmp << " ";
        }       
        cout << endl;


        // MPI_assert(MPI_Barrier(MPI_COMM_WORLD), "loadBalance: sync before partition" LOCATION);
        // call load balance

    }else{


        // send measured performance
        MPI_assert( MPI_Gather( &time, 1, MPI_DOUBLE,
                                NULL, 0, MPI_DOUBLE, 
                                0, MPI_COMM_WORLD ),
                    
                                "loadBalance: Gather failed" LOCATION );

        MPI_assert( MPI_Bcast(&balancing, 1, MPI_INT, 0, MPI_COMM_WORLD), "loadBalance bcast failed" LOCATION );

        if(balancing == 1){
            // perform balancing

            MPI_assert( MPI_Bcast(tbuf, worldSize, tdesc.TileMsg_t, 0, MPI_COMM_WORLD ), "topology Bcast failed" LOCATION);
            vector<TileDescriptor> vtd;

            for(int i = 0; i < worldSize;i++){
                vtd.push_back(TileDescriptor(tbuf[i]));
            }

            if(DBG){
                ss << vtd.at(rank);
                ss << "assigned: ";
                list<unsigned> * tmp = getAssignGIDs(vtd.at(rank));
                for(auto x : *tmp )
                    ss << x << " ";
                ss << endl;
                cout << ss.str();
                delete tmp;
            }   


            persist = resolveMigration(tdesc.getTiles(), vtd);

            // callbackDbg = true;
            migrate(vtd, *persist);
            // callbackDbg = false;


            tdesc.setTiles(vtd);
            
            if(DBG && rank == 1){

                cout << COUTLOC;
                cout << tdesc.toString();
            }

            tdesc.updateTopology();

            block = getBlockData();

            if(DBG)  synCout(tdesc.commsToString(), rank, worldSize);

        } //end imbalance 1


    } // else end

    delete[] rbuf;
    delete[] tbuf;

    if(tls != NULL) delete tls;
    if(persist != NULL) delete persist;


    return (bool) balancing;
  
    // block = ();
}


/**
 * @brief Computes which GIDs should be imported 
 * and which are persist - remains assigned to same process
 * 
 * @details Sets appropriate arrays in  LoadBalancer for migration
 * 
 */

list<unsigned> * DBD::resolveMigration(   
        const vector<TileDescriptor> & oldTiles,    // tiles in actual decomposition
        const vector<TileDescriptor> & newTiles     // tiles is new decomposition
    )
{

    stringstream ss;

    map<int, list<unsigned> *>  oldGIDs,  newGIDs;

    // get GIDs mapped to apporpriate ranks
    // in old and new decomposition
    oldGIDs = getAssignedObjs(oldTiles);
    newGIDs = getAssignedObjs(newTiles);

    vector<int> importProcs;
    vector<unsigned> importGids;

    vector<int> exportProcs;
    vector<unsigned> exportGids;



    list<unsigned> import, exp;

    list<unsigned> *  persist = new list<unsigned> ();

    oldGIDs.at(rank)->sort();
    newGIDs.at(rank)->sort();

    if(DBG){
        ss << "newGIDs:";
        for(auto g: *(newGIDs.at(rank)))
            ss << g << " ";
        ss << endl;
    
        ss << "oldGIDs:";
        for(auto g: *(oldGIDs.at(rank)))
            ss << g << " ";
        ss << endl;
    }

    
    // computer difference of old and new assigned blocks
    set_difference(
                (newGIDs.at(rank))->begin(),
                (newGIDs.at(rank))->end(),
                (oldGIDs.at(rank))->begin(),
                (oldGIDs.at(rank))->end(),
                back_inserter(import)
                );

    if(DBG){
        ss << "rank " << rank << " import:";
        for(auto d : import)
            ss << d << " ";
        ss << endl;
    }

    set_difference(
                (oldGIDs.at(rank))->begin(),
                (oldGIDs.at(rank))->end(),
                (newGIDs.at(rank))->begin(),
                (newGIDs.at(rank))->end(),
                back_inserter(exp)
                );

    if(DBG){
        ss << "rank " << rank << " export:";
        for(auto d : exp)
            ss << d << " ";
        ss << endl;
    }

    set_intersection(
                (newGIDs.at(rank))->begin(),
                (newGIDs.at(rank))->end(),
                (oldGIDs.at(rank))->begin(),
                (oldGIDs.at(rank))->end(),
                back_inserter(*persist)
                );

    if(DBG){
        ss << "rank " << rank << " persist:";
        for(auto p : *persist)
            ss << p << " ";
        ss << endl << endl;
        synCout(ss.str(), rank, worldSize);

    }



    // find where GIDs are in actual decomposition
    // and set migration data
    // TODO: candidate to optimize complexity

    for(auto p : oldGIDs){
        for(auto gid: import){

            if( find( begin(*(p.second)), end(*(p.second)), gid) != end(*(p.second)) ){
                importProcs.push_back(p.first);
                importGids.push_back(gid);
            }
        }
    }

    for(auto p: newGIDs){
        for(auto gid: exp){

             if( find( begin(*(p.second)), end(*(p.second)), gid) != end(*(p.second)) ){
                exportProcs.push_back(p.first);
                exportGids.push_back(gid);
            }
        }
    }


    if(importGids.size() > 0){

        lb.num_import = importGids.size();
      
        lb.import_global_ids = new unsigned[lb.num_import];
        lb.import_local_ids = new unsigned[lb.num_import];
        lb.import_procs = new int[lb.num_import];
        lb.import_to_part = new int[lb.num_import];
    
    
        for(int i = 0; i < lb.num_import;i++){
    
            lb.import_global_ids[i] = importGids[i];
            lb.import_local_ids[i] = importGids[i];
            lb.import_procs[i] = importProcs[i];
            lb.import_to_part[i] = 0;
    
        }
    }else{
    
        lb.num_import = 0;
        lb.import_global_ids  = NULL;
        lb.import_local_ids  = NULL;
        lb.import_procs  = NULL;
        lb.import_to_part = NULL;

    }

    if(exportGids.size() > 0){

        lb.num_export = exportGids.size();
        lb.export_global_ids =  new unsigned[lb.num_export];
        lb.export_local_ids = new unsigned[lb.num_export];
        lb.export_procs = new int[lb.num_export];
        lb.export_to_part = new int[lb.num_export];
    
        for(int i = 0; i < lb.num_export;i++){
    
            lb.export_global_ids[i] = exportGids[i];
            lb.export_local_ids[i] = exportGids[i];
            lb.export_procs[i] = exportProcs[i];
            lb.export_to_part = 0;
        }
    }else{

        lb.num_export = 0;
        lb.export_global_ids = NULL;
        lb.export_local_ids = NULL;
        lb.export_procs = NULL;
        lb.export_to_part = NULL;
    }



    return persist;
}

/**
 * @brief Perform data migration if loadBalancing is applied
 * @details [long description]
 * q
 * @param gids [description] 
 * @return [description]
 */

void DBD::migrate(  const vector<TileDescriptor> & newTiles,
                    list<unsigned> & persist
                )
{


    // find tile describing my rank
    auto r = find(begin(newTiles), end(newTiles), rank);

    if(r == end(newTiles)) throw runtime_error("migrate: rank not found in newTiles");

    TileDescriptor td = *r;


    initNewBlock(td);

    lb.migrateData();

    movePersistObj(persist, false);

    lb.clearArrays();


    // return getBlockData();
}




/**
 * @brief Collect result data to root
 * @details [long description]
 */

float * DBD::collectData(bool old)
{

    // TileDescriptor tres(rank, 0,0, edgeSize, edgeSize);

    // initNewBlock(tres);

    map<int, list<unsigned> *> objs = getAssignedObjs(tdesc.getTiles());

   

    vector<int> procs;
    vector<int> gids;


    if(rank == 0){

        /*
        for(auto item : objs){
            cout << item.first << ":";
            for(auto o: *(item.second)){
                cout << o << " ";
            }
            cout << endl;
        }
        */

        for(auto item : objs){
            if(item.first == 0)
                continue;

            for(auto gid : *(item.second)){
                procs.push_back(item.first);
                gids.push_back(gid);
            }
        }

        lb.num_import = gids.size();
        lb.import_global_ids = new unsigned[lb.num_import];
        lb.import_local_ids = new unsigned[lb.num_import];
        lb.import_procs = new int[lb.num_import];
        lb.import_to_part = new int[lb.num_import];
    
        for(unsigned i = 0; i < gids.size();i++){
    
            lb.import_global_ids[i] = gids[i];
            lb.import_local_ids[i] = gids[i];
            lb.import_procs[i] = procs[i];
            lb.import_to_part[i] = 0;
        }

        lb.num_export = 0;
        lb.export_global_ids = NULL;
        lb.export_local_ids = NULL;
        lb.export_procs = NULL;
        lb.export_to_part = NULL;

    }else{

        lb.num_export = (objs.at(rank))->size();
        lb.export_global_ids = new unsigned[lb.num_export];
        lb.export_local_ids = new unsigned[lb.num_export];
        lb.export_procs = new int[lb.num_export];
        lb.export_to_part = new int[lb.num_export];

        unsigned i = 0;
        for(auto gid: *(objs.at(rank)) ){ 

            lb.export_global_ids[i] = gid;
            lb.export_local_ids[i] = gid;
            lb.export_procs[i] = 0;
            lb.export_to_part[i] = 0;
            i++;
        }

        lb.num_import = 0;
        lb.import_global_ids = NULL;
        lb.import_local_ids = NULL;
        lb.import_procs = NULL;
        lb.import_to_part = NULL;

    }

    if(rank == 0){
        for(unsigned i = 0; i < edgeSize * edgeSize;i++)
            resArray[i] = 0.0;
    }
    

    // swap(bdata.oldTemp, bdata.newTemp);

    oldArray = old; // bool arg 
    collectDataFlag = true; // callbacks write to resArray
    lb.migrateData();
    collectDataFlag = false;

    // move data on rank  0 to result
    if(rank == 0){
        movePersistObj(*objs.at(rank), false, true);
    }


    oldArray = false;

    // if(this->rank == 0){
// 
       // res = new float[edgeSize * edgeSize];
       // 
        // for(int i = 0; i < edgeSize;i++){
            // for(int j = 0; j < edgeSize;j++){
                // res[i*edgeSize + j] = bdata.oldTemp[halo(j,i, edgeSize+4)];
            // }
        // }
    // }

    lb.clearArrays();

    for(auto i : objs){
        delete i.second;
    }

    if(rank == 0)
        return resArray;
    else
        return NULL;

}

// int DBD::getRank(const Dims & blockPosition)
// {
    // 
// }

string DBD::toString(void)
{

    stringstream ss;
      
    ss << "========= PRINTING DBD INFO ===========" << endl;
    ss << "Rank: " << rank << endl;
    ss << "Size: " << worldSize << endl;
    ss << "Object size" << objectSize << endl;
    ss << "Cols: " << cols << " Rows: " << rows << endl;
    ss << "Neighbors: " << bdata.neighbors->size() << endl;
    ss << "blockSize " << tdesc.tile().getSize() << endl;
    ss << "blockPosition " << tdesc.tile().getPosition() << endl; 
    ss << "extSize " << tdesc.tile().getExtSize() << endl;
    
    ss << "Bounds: [T,R,B,L]" << "[" << bdata.top << ", " << bdata.right << ", ";
    ss << bdata.bottom << ", " << bdata.left << "]" << endl;
    
    ss << "isMiddle: " << bdata.middle << endl;
    ss << "===============================" << endl;
    
    return ss.str();
}

string DBD::tilesToString(void)
{
    return tdesc.toString();
}

void DBD::checkRank(int rank)
{
    if(rank < NO_RANK || rank > worldSize ) 
        throw runtime_error(string("CheckRank: Invalid rank") + to_string(rank));
}


/**
 * @brief Stores and return actual block data
 * @details Topology metadata are retrieved from
 *           TopologyDescriptor, array pointers are added
 *           and whole new block is stored and returned
 * 
 * @return Blockdata
 */

DLB::BlockData DBD::getBlockData(void)
{
   BlockData b = tdesc.getBlockData();

   b.oldTemp = bdata.oldTemp;
   b.newTemp = bdata.newTemp;
   b.domParams = bdata.domParams;
   b.domMap = bdata.domMap;

   //store bdata as attribute
   bdata = b; 

   return bdata;

}

/**
 * @brief Returns length of halo zones
 * 
 * @return [description]
 */

unsigned DBD::getHaloLen(void)
{ 
    Dims tmp = tdesc.tile().getSize();

    return 2*tmp.x + 2*tmp.y;
}



DLB::TileDescriptor DBD::resolvePartition(unsigned * gids, unsigned size)
{
    list<unsigned> tmp;
    for(unsigned i = 0; i < size;i++)
        tmp.push_back(gids[i]);

    return resolvePartition(tmp);
}

DLB::TileDescriptor DBD::resolvePartition(list<unsigned> & assigned)
{   
    TileDescriptor d;

    //sort GIDs
    assigned.sort();

    Dims first = getCoordsByGID(assigned.front());
    Dims last = getCoordsByGID(assigned.back());

    Dims size = (last + objectSize) - first; 

    d = tdesc.tile();
    d.setPosition(first);
    d.setSize(size);

    // how many objects
    unsigned objs = (size.x * size.y) / (objectSize.x * objectSize.y);

    if(assigned.size() != objs){
        stringstream ss ;
        ss << rank;
        throw runtime_error(ss.str() +  ": Object amount does not match assigned block");
    }

    return d;
}

/**
 * @brief [brief description]
 * @details Results in update bdata arrays
 * 
 * @param persist [description]
 * @param init [description]
 */

void DBD::movePersistObj(const list<unsigned> & persist, bool init, bool collect)
{
    // if(persist == NULL) throw runtime_error("persist = NULL");

    if(init){

        // dont have any previous data
        // eg. BlockData are not initialized
        bdata.oldTemp = newBlock.temp;
        bdata.newTemp = new float[newBlock.tile.getExtArea()];

        for(unsigned i = 0; i < newBlock.tile.getExtArea();i++){
            bdata.newTemp[i] = bdata.oldTemp[i];
        }

        bdata.domParams = newBlock.params;
        bdata.domMap = newBlock.map;

        newBlock.temp = NULL;
        newBlock.params = NULL;
        newBlock.map = NULL;


    }else{
        // have previous data

        // vector<unsigned> * objs = getAssignedObjs();

        Dims relOld, relNew;
        Dims esizeOld = tdesc.tile().getExtSize();
        Dims esizeNew;

        if(collect){

             for(auto obj : persist){

                Dims objPos = getCoordsByGID(obj);
                relOld = objPos - tdesc.tile().getPosition() ;
        
                for(unsigned i = 0; i < objectSize.y;i++){
                    for(unsigned j = 0 ; j < objectSize.x;j++){
                        // resolve coordinates relative to block icncluding halo zones
                        unsigned oldIdx = halo( relOld.x + j, relOld.y + i, esizeOld.x);

                        // resArray without halo zones
                        unsigned newIdx = (objPos.y + i) * edgeSize + (objPos.x + j);
                          // move data from old arrays to new ones
                        if(oldArray){
                            resArray[newIdx] = bdata.oldTemp[oldIdx];
                        }else{
                            resArray[newIdx] = bdata.newTemp[oldIdx];

                        }
                    }
                }
            }


        }else{

            Dims esizeNew = newBlock.tile.getExtSize();

            for(auto obj : persist){
    
                relNew = getCoordsByGID(obj) - newBlock.tile.getPosition() ;
                relOld = getCoordsByGID(obj) - tdesc.tile().getPosition() ;
        
                for(unsigned i = 0; i < objectSize.y;i++){
                    for(unsigned j = 0 ; j < objectSize.x;j++){
                        // resolve coordinates relative to block icncluding halo zones
                        unsigned oldIdx = halo( relOld.x + j, relOld.y + i, esizeOld.x);
                        unsigned newIdx = halo( relNew.x + j, relNew.y + i, esizeNew.x); 
                          // move data from old arrays to new ones
                        newBlock.temp[newIdx] = bdata.oldTemp[oldIdx];
                        newBlock.params[newIdx] = bdata.domParams[oldIdx];
                        newBlock.map[newIdx] = bdata.domMap[oldIdx];
                    }
                }
            }

            // copy persistent objects
           
            // delete old array, actual are in new block
            delete[] bdata.newTemp;
            delete[] bdata.oldTemp;
            delete[] bdata.domParams;
            delete[] bdata.domMap;
    
            bdata.oldTemp = newBlock.temp;
            // make new arrays actual
            bdata.newTemp = new float[newBlock.tile.getExtArea()];
            for(unsigned i = 0; i < newBlock.tile.getExtArea();i++){
                bdata.newTemp[i] = bdata.oldTemp[i];
            }
            bdata.domParams = newBlock.params;
            bdata.domMap = newBlock.map;
    
        }
    }

}

/**
 * @brief [brief description]
 * @details [long description]
 * @return [description]
 */

float DBD::middleColAvg(void)
{
    unsigned mid = edgeSize / 2;

    unsigned offset = mid -  tdesc.tile().getPosition().x;

    float avg = 0.0;
    for(unsigned i = 0; i < tdesc.tile().getSize().y; i++){

        avg += bdata.newTemp[ halo(offset, i, tdesc.tile().getExtSize().x) ];
    }

    // avg = avg / tdesc.tile().getSize().y;

    return avg;
}



/**
 * Zoltan callbacks and related functions
 * 
 */

/**
* returns number of objects assigned to current CPU
*/

int DBD::zolt_num_obj_fn(void * data, int * err)
{
    unused(err);
    return ((DynamicBlockDescriptor *) data)->assignedObjsCnt;
}


/**
* Returns number of dimensions needed to express geometry
*/

int DBD::zolt_num_geom_fn(void * data, int *err)
{
    unused(data);
    unused(err);

    return 2;
}

/**
 *  @brief Compute object location in mesh based on 
 *  given object ID.
 * 
 *  @detailed Object IDs are given in list
 *  
 */

void DBD::zolt_geom_multi_fn( void * data, 
    int num_gid_entries, 
    int num_lid_entries, 
    int num_obj,
    ZOLTAN_ID_PTR global_ids, 
    ZOLTAN_ID_PTR local_ids,
    int num_dim,
    double *geom_vec,
    int * ierr
    )
{

    unused(local_ids);

    DynamicBlockDescriptor * dbd =  (DynamicBlockDescriptor *) data;

    if(num_gid_entries != 1 || num_lid_entries != 1){
        *ierr  = ZLT_ID_TOO_LONG;
        return;
    }      

    if(num_dim != 2){
        *ierr = ZLT_UNSUP_DIMS;
        return;
    }

    Dims tmp;

    // geom_vect will contain object coordinates 
    // computed from regular mesh
    for(int i = 0; i < num_obj;i++){

        tmp = dbd->getCoordsByGID(global_ids[i]);

        geom_vec[i*num_dim] = (double) tmp.x;
        geom_vec[i*num_dim+1] = (double) tmp.y;
    }
}


/**
* @brief Returns position of single object.
*/

void DBD::zolt_geom_fn( void * data, 
    int num_gid_entries, 
    int num_lid_entries, 
    ZOLTAN_ID_PTR global_id, 
    ZOLTAN_ID_PTR local_id,
    double *geom_vec,
    int * ierr
    )
{

    unused(local_id);

    if(num_gid_entries != 1 || num_lid_entries != 1){
        *ierr  = ZLT_ID_TOO_LONG;
        return;
    }

    DynamicBlockDescriptor * dbd =  (DynamicBlockDescriptor *) data;


    Dims tmp;

    // geom_vect will contain object coordinates 
    // computed from regular mesh
    tmp = dbd->getCoordsByGID(*global_id);

    geom_vec[0] = (double) tmp.x;
    geom_vec[1] = (double) tmp.y;

}

/**
*   @brief Returns object size in bytes
*/
int DBD::zolt_obj_size_fn(void * data,
    int num_gid_entries,
    int num_lid_entries,
    ZOLTAN_ID_PTR global_id,
    ZOLTAN_ID_PTR local_id,
    int *ierr
   )
{

    unused(global_id);
    unused(local_id);

    if(num_gid_entries != 1 || num_lid_entries != 1){
        *ierr  = ZLT_ID_TOO_LONG;
        return 0;
    }

    DynamicBlockDescriptor * p = (DynamicBlockDescriptor *) data;

    unsigned objarea = p->objectSize.x * p->objectSize.y;

    if(p->collectDataFlag){

        return objarea * sizeof(float);

    }else{

        return 2*(objarea * sizeof(float)) + (objarea * sizeof(int));
    }

}

/**
 * @brief Extract data object from actual block
 * @detailed Reads data from oldTemp array
 * 
 */

void DBD::zolt_pack_obj_fn( void * data,
    int num_gid_entries,
    int num_lid_entries,
    ZOLTAN_ID_PTR global_id,
    ZOLTAN_ID_PTR local_id,
    int dest, 
    int size,
    char * buf,
    int *ierr
    )
{

    unused(size);
    unused(local_id);
    unused(dest);


    if(num_gid_entries != 1 || num_lid_entries != 1){
        *ierr  = ZLT_ID_TOO_LONG;
        return;
    }

    DynamicBlockDescriptor * dbd = (DynamicBlockDescriptor *) data;


    if(dbd->callbackDbg) cout << "rank "  << dbd->getRank() << " pack: " << *global_id << endl;

    Dims esize = dbd->tdesc.tile().getExtSize();
    Dims blockPos = dbd->tdesc.tile().getPosition();

    Dims objPos = dbd->getCoordsByGID(*global_id);

    Dims rel = objPos - blockPos; //relative object position

    // Point_t pt;
    // void *fbuf = NULL;
    // Point_t * pbuf = NULL;

    // if(dbd->collectDataFlag){
        // fbuf = (float *) buf;
    // }else{
        // pbuf = (Point_t *) buf;
    // }


    // extract object from block int Point_t structure
    // which will be passed to buffer

    float *fbuf, *fdata;
    int * ibuf, *idata;

    if(dbd->collectDataFlag){

        fbuf = (float *) buf;

         for(unsigned i = 0; i < dbd->objectSize.y;i++){
        // for(unsigned j = 0 ; j < dbd->objectSize.x;j++){

            unsigned di = halo(rel.x, rel.y + i, esize.x );
            //((rel.y + i + HALO_SIZE)*esize.x) + (rel.x + j + HALO_SIZE);
            unsigned bi = i* dbd->objectSize.x ;

            std::memcpy(&(fbuf[bi]),
                       &(dbd->bdata.newTemp[di]), 
                       dbd->objectSize.x*sizeof(float));
        }


    }else{

        unsigned size;
        unsigned offset;

        for(int x = 0; x < 3;x++){
            if(x == 0){
                fdata = dbd->oldArray ? dbd->bdata.oldTemp : dbd->bdata.newTemp;
                fbuf = reinterpret_cast<float *>(buf);
                size = sizeof(float);
                offset = 0;

            }else if(x == 1){
                offset = dbd->objectSize.x * dbd->objectSize.y;
                fdata = dbd->bdata.domParams;
                fbuf =  reinterpret_cast<float *>(buf);
                fbuf += offset;
                size = sizeof(float);

            }else if(x == 2){

                offset =  2*dbd->objectSize.x * dbd->objectSize.y;
                idata = dbd->bdata.domMap;
                fbuf = reinterpret_cast<float *>(buf); //second array offset
                fbuf += offset;
                ibuf = reinterpret_cast<int*>(fbuf);
                size = sizeof(int);
            }


             for(unsigned i = 0; i < dbd->objectSize.y;i++){
            // for(unsigned j = 0 ; j < dbd->objectSize.x;j++){
    
                unsigned di = halo(rel.x, rel.y + i, esize.x );
                //((rel.y + i + HALO_SIZE)*esize.x) + (rel.x + j + HALO_SIZE);
                unsigned bi = i* dbd->objectSize.x ;

                // if( (unsigned char *) &(ibuf[bi]) > ((unsigned char *) buf) + memsize || 
                    // (unsigned char *) &(fbuf[bi]) > ((unsigned char *) buf) + memsize){
                    // throw runtime_error("pack_obj memory hazard");
                // }

                if(x == 2){

                     std::memcpy( &(ibuf[bi]),
                                  &(idata[di]),
                                  dbd->objectSize.x*size
                                );

                }else{
                     std::memcpy(&(fbuf[bi]),
                                 &(fdata[di]),
                                 dbd->objectSize.x*size
                                );

                }
    
            }

        }

    }
}

/**
* @brief Unpack data from Zoltan to block
* @detailed Extracts data into allocated newblock
*/

void DBD::zolt_unpack_obj_fn( 
                            void * data,       
                            int num_gid_entries,  
                            ZOLTAN_ID_PTR global_id,
                            int size,
                            char *buf,
                            int *ierr
                            )

{

    unused(size);

    if(num_gid_entries != 1){
        *ierr  = ZLT_ID_TOO_LONG;
        return;
    }

    DynamicBlockDescriptor *dbd = (DynamicBlockDescriptor *) data;
    unsigned gid = * global_id;

    // if(rank == 12 ||  rank == 13 || rank == 14 ||  rank == 15 )
        // cout << rank << ": " << gid << endl;

    // if(gid >= dbd->totalObjsCnt){
        // cout << dbd->rank << ": " << gid << endl;
        // return;
    // }

   // unsigned objarea = dbd->objectSize.x * dbd->objectSize.y; 
    // unsigned memsize = 2*(objarea * sizeof(float)) + (objarea * sizeof(int));


    if(dbd->callbackDbg) cout << "rank "  << dbd->getRank() << " unpack: " << *global_id << endl;

    // Dims blockPos = dbd->tdesc.tile().getPosition();
    Dims objPos = dbd->getCoordsByGID(gid);
    
    Dims blockPos, esize;


    float * fdata, *fbuf;
    fdata = fbuf = NULL;
    int * idata, * ibuf;
    idata = ibuf = NULL;

    // collecting data to write
    if(dbd->collectDataFlag){

        blockPos = Dims(0,0);
        esize = Dims(dbd->edgeSize, dbd->edgeSize);

        Dims rel = objPos - blockPos; //relative object position


        unsigned size = sizeof(float);
        fbuf = (float *) buf;

        for(unsigned i = 0; i < dbd->objectSize.y;i++){
        // for(unsigned j = 0 ; j < dbd->objectSize.x;j++){

            // offset in object (source) buffer
            unsigned bi = i* dbd->objectSize.x; 
            unsigned di = (rel.y + i) * esize.x + (rel.x); 

            //copy object data from buffer to block eg. unpack
            if(dbd->rank != 0) throw runtime_error("unpack: possible segfault");
            // target array withou halo zone
            std::memcpy( &(dbd->resArray[di]), 
                         &(fbuf[bi]),
                        dbd->objectSize.x * size
                        );
            // dbd->resArray[di] = fbuf[bi];

        }

    }else{

        // unpacking to arbitraty block during domain mapping
        blockPos = dbd->newBlock.tile.getPosition();
        esize = dbd->newBlock.tile.getExtSize();

        Dims rel = objPos - blockPos; //relative object position


        unsigned size;
        unsigned offset;


        for(int x = 0; x < 3;x++){

            if(x == 0){
                offset = 0;
                fdata =  dbd->newBlock.temp;
                fbuf = reinterpret_cast<float *>(buf);
                size = sizeof(float);
            }else if(x == 1){

                offset = dbd->objectSize.x * dbd->objectSize.y ;
                fdata =  dbd->newBlock.params;
                // fbuf = (float *) (buf + offset);
                fbuf = reinterpret_cast<float*>(buf);
                fbuf += offset;
                size = sizeof(float);

            }else if(x == 2){
                offset = 2*dbd->objectSize.x * dbd->objectSize.y;
                idata = dbd->newBlock.map;
                fbuf = reinterpret_cast<float*>(buf);
                fbuf += offset;
                ibuf = reinterpret_cast<int*>(fbuf);
                size = sizeof(int);
            }


             for(unsigned i = 0; i < dbd->objectSize.y;i++){
            // for(unsigned j = 0 ; j < dbd->objectSize.x;j++){
    
                unsigned di = halo(rel.x, rel.y + i, esize.x );
                //((rel.y + i + HALO_SIZE)*esize.x) + (rel.x + j + HALO_SIZE);
                unsigned bi = i* dbd->objectSize.x ;

               // if( (unsigned char *) &(ibuf[bi]) > ((unsigned char *) buf) + memsize || 
                    // (unsigned char *) &(fbuf[bi]) > ((unsigned char *) buf) + memsize){
                    // throw runtime_error("pack_obj memory hazard");
                // }
    
                if(x == 2){
                    std::memcpy(  &(idata[di]), 
                                  &(ibuf[bi]) ,
                                  dbd->objectSize.x*size
                                );

                }else{

                     std::memcpy(  &(fdata[di]), 
                                  &(fbuf[ bi]) ,
                                  dbd->objectSize.x*size
                                );


                }
            }

        }
    }

}

/**
 * @brief 
 * @details [long description]
 * 
 * @param  [description]
 * @return [description]
 */

map<int, list<unsigned>* >  DBD::getAssignedObjs(const vector<TileDescriptor> & tds)
{
    map<int, list<unsigned> *>   lst;

    for(auto& tile : tds){

        lst.insert(pair<int, list<unsigned> *>( tile.getRank(), getAssignGIDs(tile) ) );
    } 

    return lst;
}

/**
* @brief Return list of assigned objects
* @return new list of assigned GIDs
*/


list<unsigned> * DBD::getAssignGIDs(const TileDescriptor & t)
{

    list<unsigned>  * lst = new list<unsigned>();

    Dims blockPos = t.getPosition();
    Dims esize = t.getSize();

    // int idx = 0;
    // increment by object size on y axis
    for(unsigned y = blockPos.y; y < blockPos.y + esize.y; y+= objectSize.y){
        for(unsigned x = blockPos.x; x < blockPos.x + esize.x; x+= objectSize.x){

            lst->push_back(getGIDByCoords(Dims(x, y)));
        }
    }

    return lst;
}


list<unsigned> * DBD::getAssignGIDs(void)
{

    list<unsigned> * lst = new list<unsigned>();

    Dims blockPos = tdesc.tile().getPosition();
    Dims esize = tdesc.tile().getSize();


    // int idx = 0;
    // increment by object size on y axis
    for(unsigned y = blockPos.y; y < blockPos.y + esize.y; y+= objectSize.y){
        for(unsigned x = blockPos.x; x < blockPos.x + esize.x; x+= objectSize.x){

            lst->push_back(getGIDByCoords(Dims(x, y)));
        }
    }

    return lst;
}



/**
* @brief Create list of object IDs associated with this process
*/

void DBD::zolt_obj_list_fn(     void *data,
    int num_gid,
    int num_lid, 
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int wgt_dim,
    float * obj_wgts,
    int *ierr)
{

    unused(num_gid);
    unused(num_lid);
    unused(local_ids);
    unused(wgt_dim);
    unused(obj_wgts);
    unused(ierr);

    DynamicBlockDescriptor * dbd = (DynamicBlockDescriptor *) data;

    list<unsigned>  *ids = dbd->getAssignGIDs();

    unsigned idx = 0;
    for(auto obj : *ids){

        global_ids[idx] = obj;
        idx++;
    }

    delete ids;
}


/**
 * @brief Returns GID computed from object position
 *  relative to object size
 
 * @details GID are assigned to every object from upper
 *          left corner ascendingly. Upper-left is 0.
 * 
 * @param objPos zoltan obj position
 * @return GID
 */
unsigned inline DBD::getGIDByCoords(const Dims & objPos)
{

    int col = objPos.x / objectSize.x ;
    int row = objPos.y / objectSize.y; 

    return (row * cols) + col;

}

DLB::Dims inline DBD::getCoordsByGID(unsigned gid)
{
    Dims p;

    p.x = (gid % cols)  * objectSize.x;
    p.y = ((gid - (gid % cols)) / cols ) * objectSize.y;

    return p;
}
