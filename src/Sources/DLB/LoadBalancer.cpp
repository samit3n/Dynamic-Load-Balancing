/***********************************************
*
*  File Name:       LoadBalancer.cpp
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Main class implementing dynamic load balancing
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            22.5.2017
*
***********************************************/

#include "LoadBalancer.h"

using LoadBalancer = DLB::LoadBalancer;
using TileDescriptor = DLB::TileDescriptor;
using std::runtime_error;
using std::pair;
using std::list;
using std::min_element;
using std::max_element;
using std::accumulate;
using std::for_each;


LoadBalancer::LoadBalancer(int rank, size_t edgeSize,  int worldSize, Dims objectSize, double threshold):
edgeSize(edgeSize),
rank(rank),
worldSize(worldSize),
objectSize(objectSize),
threshold(threshold)
{
    // init empty values
    zz = NULL; // Zoltan null until created from DBD::zoltanInit()
    
    changed = -1;
    gid_entries = lid_entries = -1;
    num_import = -1;
  
    import_global_ids =  import_local_ids = NULL;
    export_global_ids =  export_local_ids = NULL;
  
    import_procs = NULL;
    import_to_part = NULL;
    export_procs = NULL;
    export_to_part = NULL;
    num_export = -1;

    // object mesh dimensions

    objCols = edgeSize / objectSize.x;
    objRows = edgeSize / objectSize.y;

    rows = cols = 0; //will be updated by regularTiles()
    imbalance = false;

}

LoadBalancer::~LoadBalancer(void)
{
    if(zz != NULL) delete zz;
}


/**
 * @brief Load balance and returned topology
 * @details Sets new relative part sizes to Zoltan, based on run stats
 * 
 * @param stats rank <-> runtime pairs
 * @param tds - actual tile decomposition to be updated
 * @return new decomposition
 */

void inline LoadBalancer::check(double& x, bool max){

    // Zoltan returns +/- DBL_MAX for outer boundaries

    if(x < 0)
        x = 0.0;

    if(x > edgeSize)
    {
        if (max){
            x = (double) edgeSize;
        }else{
        x = (double) (edgeSize - objectSize.x);
        }
    }
}

/**
 * @brief Migrate object between procs based on zoltan migrate
 * @details Appropriate arrays must be set before call
 */

void LoadBalancer::migrateData(void)
{

    MPI_assert( MPI_Barrier(MPI_COMM_WORLD), "migrateData barrier" LOCATION);

    MPI_assert( zz->Migrate(
                num_import,
                import_global_ids,
                import_local_ids,
                import_procs,
                import_to_part,
                num_export,
                export_global_ids,
                export_local_ids,
                export_procs,
                export_to_part

              ), "Migrate failed" LOCATION );

    MPI_assert( MPI_Barrier(MPI_COMM_WORLD), "migrateData barrier" LOCATION);

}



/**
 * @brief Generate regular mesh
 * 
 */
vector<TileDescriptor> * LoadBalancer::regularTiles(void)
{

    vector<TileDescriptor> * tls = new vector<TileDescriptor>();


    Dims blockSize;

    // number of cpu's is Even power of 2
    // defines if work will be divied by squares or rectangles

    // since size is even or odd power of 2, conversion to int
    // is safe

    bool isCpuEven = ((int) log2(worldSize)) % 2 == 0;

    if(isCpuEven){

        // assigned block will be square
        blockSize.x = blockSize.y = edgeSize / (int) sqrt(worldSize);

    }else{


        // assigned block will be rectangle
        // we can get X-size of block the same as in previous case
        // makin number of cpus even power of 2
        // y size is then twice that long to make rectangle

        int tmp = worldSize*2;
        blockSize.x = edgeSize / (int) sqrt(tmp);
        blockSize.y = blockSize.x * 2;
    }


    // these will ramin fixed during simulation
    // only size of blocks in row changes
    cols = edgeSize / blockSize.x; // blocks in row = cols
    rows = edgeSize / blockSize.y; // blocks in col = rows

  // working directly on vector

    int tmpRank = 0;
    for(unsigned r = 0; r < rows;r++){
        for(unsigned c = 0 ; c < cols;c++){

            TileDescriptor tempTd;

            tempTd.setRank(tmpRank);
            tempTd.setPosition( Dims(c * blockSize.x, r * blockSize.y) );
            tempTd.setSize( Dims( blockSize.x, blockSize.y) );
            tls->push_back(tempTd);

            tmpRank++;

        }
    }

    if(tmpRank > worldSize + 1)
        throw runtime_error("regularTiles: rank out of bounds");


    return tls;
}



bool LoadBalancer::isBalanced(vector<double> & times)
{
    if(DBG){
        cout << COUTLOC << endl << "times:" << endl;
        for(auto t: times){
            cout << t << " "; 
            cout << endl;
        }
    }

    // normalize times, when imbalanced


    // double avg = accumulate(times.begin(), times.end(), 0.0 )/ times.size();    

    auto max = max_element(times.begin(), times.end());
    auto min = min_element(times.begin(), times.end());

    // cout << "isBalanced: min = " << *min << " max = "  << *max << endl;
    // if(*max - *min > avg * 0.2){
    if(*max >  (*min * threshold) ){

        // cout << COUTLOC << "Imbalance detected" << endl;
        return false;

    }else{

        return true;
    }
}

/**
 * @brief   Split number of object in one row by measured performance
 * @details 
 * 
 * @param times performance of processes
 * @return [description]
 */

vector<unsigned> LoadBalancer::splitByPerform(const vector<double> & times)
{   
    vector<unsigned> objects;
    vector<unsigned> sizes;
    vector<double> normalized = times;

    // cout << "cols" << objCols << endl;


    // cout << "times:";
    // for(float x: times) cout << std::fixed << std::setprecision(4) << x << " "; cout << endl;

    // cout << "normalized:";
    for_each(normalized.begin(), normalized.end(), [](double & t){ t = 1.0 / t;});

    // for(float x: normalized) cout << std::fixed << std::setprecision(2) << x << " "; cout << endl;



    // for(float x: normalized) cout << std::fixed << std::setprecision(2) << x << " "; cout << endl;

    double sum = accumulate(normalized.begin(), normalized.end(),0.0);
    double unit = objCols / sum; // unit amount of objects relative to norm. values

    for_each(normalized.begin(), normalized.end(), [unit](double &t){ t = t * unit; });

    // for(float x: normalized) cout << std::fixed << std::setprecision(2) << x << " "; cout << endl;


    for(auto x: normalized){
        
        int size = round(x);
        size = size == 0 ? 1 : size;
        sizes.push_back(size);
    }
    // for(float x: sizes) cout  << x << " "; cout << endl;

    unsigned objSum = accumulate(sizes.begin(), sizes.end(), 0);

    // if sum of assigned objects do not match, simpley spread the rest between
    // may be optimized by selecting parts

    if(objSum != objCols){
        int dif = objSum - objCols;
        unsigned absval = abs(dif);

        if(dif > 0){ 
            //more object assigned than exist
            //more likely because of ceiling

            for(unsigned d = 0; d < absval;d++){
                sizes[d % sizes.size()] -= 1;
            }
        }else{ //some objects not assigned

            for(unsigned d = 0; d < absval;d++){
               sizes[d % sizes.size()] += 1;
            }
        }
    }

    return sizes;

}

vector<TileDescriptor> * LoadBalancer::getPartition(  const vector<double> & times,
                                        const vector<TileDescriptor> & tiles
                                    )
{
    // objCols = edgeSize / objectSize.x;
    // objRows = edgeSize / objectSize.y;

    vector<TileDescriptor> byRow;
    vector<TileDescriptor> * newTiles = new vector<TileDescriptor>();
    // vector<TileDescriptor> row;
    vector<unsigned> sizeOnRow;;
    
    auto it = times.begin();
    auto tileIt = tiles.begin();
    unsigned pos = 0;

    // cout << COUTLOC << endl;
    // cout << "=Balance=" << endl;

    // cout << "getPart times: ";
    // for(auto t : times) cout << t << " "; cout << endl;

    for(unsigned i = 0; i < rows;i++){

        sizeOnRow = splitByPerform(vector<double>(it,it+cols));

        // cout << "row: " << i << " ";
        // for(auto s: sizeOnRow) cout << s << ",";
        // cout << endl;

        for(auto sz : sizeOnRow){

            TileDescriptor temp = *tileIt;
            tileIt++;

            //modify original tile
            Dims p = temp.getPosition();
            Dims s = temp.getSize();
            p.x = pos;
            s.x = sz * objectSize.x;
            temp.setPosition(p);
            temp.setSize(s);
            // push modified tile
            newTiles->push_back(temp);
            pos += sz * objectSize.x; 
        }
        pos = 0; //begin new line
        it += cols;
    }
    // cout << "==========" << endl;
 

    return newTiles;

}




void LoadBalancer::clearArrays(void)
{
    if(import_global_ids != NULL){
        delete[] import_global_ids;
        import_global_ids = NULL;
    }


    if(import_local_ids != NULL){
        delete[] import_local_ids;
        import_local_ids = NULL;
    }

    if(export_global_ids != NULL){
        delete[] export_global_ids;
        export_global_ids = NULL;
    }

    if(export_local_ids != NULL){
        delete[] export_local_ids;
        export_local_ids = NULL;
    }

    if(import_procs != NULL){
        delete[] import_procs;
        import_procs = NULL;
    }

    if(import_to_part != NULL){
        delete[] import_to_part;
        import_to_part = NULL;
    }
    
    if(export_procs != NULL){
        delete[] export_procs;
        export_procs = NULL;
    }

    if(export_to_part != NULL){
        delete[] export_to_part;
        export_to_part = NULL;
    }

}