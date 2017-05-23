/***********************************************
*
*  File Name:       TopologyDescriptor.cpp
*
*  Project:         Dynamic Load Balancing in HPC Applications
*                   DIP (SC@FIT)
*                 
*  Description:     Topology information class implementation
*           
*  Author:          Vojtech Dvoracek
*  Email:           xdvora0y@stud.fit.vutbr.cz
*  Date:            7.3.2017
*
*******************************************/

#include "TopologyDescriptor.h"


using TopologyDescriptor = DLB::TopologyDescriptor;
using TileDescriptor = DLB::TileDescriptor;
using Dims = DLB::Dims;
using TEdge = DLB::TileDescriptor::TEdge;


TopologyDescriptor::TopologyDescriptor(int rank, int worldSize, size_t edgeSize):
rank(rank),
worldSize(worldSize),
edgeSize(edgeSize)
{
    topologyChanged = false;

    // init tiles

    tiles.push_back(TileDescriptor(0, 0,0, edgeSize, edgeSize ));
    for(int i = 1; i < worldSize;i++){
        tiles.push_back(TileDescriptor(i,0,0,0,0));
    }
    myTile = tiles[rank];

     //initialize MPI datatypes
     initDtypes(); 

     myComm = MPI_COMM_NULL;
}

TopologyDescriptor::~TopologyDescriptor(void)
{
    MPI_assert( MPI_Comm_free(&myComm) LOCATION);
}


void TopologyDescriptor::initDtypes()
{
    // TileMsg datatype

    MPI_Datatype types[] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_INT, MPI_INT};

    int lens[] = {1,1,1,1,1,1};
    MPI_Aint displs[6];

    TileMsg tm;


    displs[0] = ( (int *) (&tm.posx) -  (int *) (&tm)) * sizeof(int);
    displs[1] = ( (int *) (&tm.posy) -  (int *) (&tm)) * sizeof(int);
    displs[2] = ( (int *) (&tm.dimx) -  (int *) (&tm)) * sizeof(int);
    displs[3] = ( (int *) (&tm.dimy) -  (int *) (&tm)) * sizeof(int);
    displs[4] = ( (int *) (&tm.rank) -  (int *) (&tm)) * sizeof(int);
    displs[5] = ( (int *) (&tm.hostNumber) -  (int *)(&tm)) * sizeof(int);


    MPI_assert(MPI_Type_create_struct(6, lens, displs, types, &TileMsg_t), "Create struct failed" LOCATION) ;
    MPI_assert( MPI_Type_commit(&TileMsg_t), "Type commit failed" LOCATION);

}


void TopologyDescriptor::setTiles(TileMsg * tmsgs )
{
    if(tmsgs == NULL){
        stringstream ss;
        ss << "TopologyDescriptor::setTiles, rank " << rank <<  " NULL pointer passed";
        throw std::runtime_error(ss.str());
    }

    tiles.clear();

    for(int  i = 0; i < worldSize; i++){

        tiles.push_back(TileDescriptor(tmsgs[i]));
        if( tmsgs[i].rank == rank){
            myTile = tiles.back();
        }
    }

}

void TopologyDescriptor::setTiles(const vector<TileDescriptor> & tds )
{
    tiles = tds;
    auto my = find(tds.begin(), tds.end(), rank);
    myTile = *my;
}


vector<TileDescriptor> & TopologyDescriptor::getTiles(void)
{
    return tiles;
}

Dims TopologyDescriptor::getPosition(void) const
{
    return myTile.getPosition();
}

Dims TopologyDescriptor::getSize(void) const
{
    return myTile.getSize();
}

unsigned TopologyDescriptor::getDomainArea(void) const
{
    return edgeSize * edgeSize;
}

int TopologyDescriptor::getRank(const Dims & position) const
{
    for(auto tile : tiles)
    {
        if(tile.getPosition().x == position.x && tile.getPosition().y == position.y){
            return tile.getRank();
        }
    }
    return TopologyDescriptor::NO_RANK;
}



void TopologyDescriptor::updateTopology(void)
{


    if(myComm != MPI_COMM_NULL)
        MPI_Comm_free(&myComm);

    nData.clear();

    neighbors.clear();
    displs.clear(); 
    counts.clear();

    // order must be preserved
    countNeighborRanks();
    countDispls();
    initComms(); //dependent on neighbors
    midUpdate();




}

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param  [description]
 */
void TopologyDescriptor::countNeighborRanks(void)
{
    neighbors.clear();

    for(auto t : tiles){

        if( myTile != t && myTile.isNeighbor(t) )        
            neighbors.push_back(t.getRank());

    }

    sort(neighbors.begin(), neighbors.end());
}

void TopologyDescriptor::countDispls(void)
{
	unsigned cnt;
    Neighbor nbor;

	for(auto n: neighbors){
		// find neighbor interator by rank
	    auto res = find(tiles.begin(), tiles.end(), n);
	    // get neighbor index
	    unsigned idx = res - tiles.begin();
        unsigned dsp = myTile.getOverlapOffset(tiles[idx], cnt);

        nbor.displ = 2*dsp;
        nbor.count = 2*cnt;
        nbor.wRank = n;

        displs.push_back(2*dsp);
        counts.push_back(2*cnt);

        nData.insert(pair<int, Neighbor>(n, nbor ));
        // nbor.clear();

	}

    // insert counts and dipls for scatter for this tile
    // must be put to appropriate place relative to neighbors
    bool once = true;
    for(auto it = neighbors.begin(); it < neighbors.end();it++){

        if(once && (*neighbors.begin()) > rank){

            counts.insert(counts.begin(), 0);
            displs.insert(displs.begin(), 0);
            once = false;
            
        }else if(once &&  *(neighbors.end() - 1) < rank){

            counts.insert(counts.end(), 0);
            displs.insert(displs.end(), 0);
            once = false;

        }else if(once &&  (*it < rank) && ((it+1) != neighbors.end()) && (*(it+ 1) > rank) )
        {   //somewhere inside
            unsigned idx = it - neighbors.begin();
            counts.insert(counts.begin() + idx +1, 0);
            displs.insert(displs.begin() + idx +1, 0);
            once = false;
        }
    }
}


void TopologyDescriptor::initComms(void)
{
    
    MPI_Comm newcomm;
    MPI_Group world, tmpGroup;
    stringstream ss, css;

    int * scatCnts, *scatDisp;

    // create comm group for rank mapping
    MPI_assert( MPI_Comm_group(MPI_COMM_WORLD, &world), "Comm_group_world" LOCATION);

    //iterate over ranks

    for(int i = 0; i < worldSize;i++){

        auto res = find(neighbors.begin(), neighbors.end(), i);

        if(res != neighbors.end() || i == this->rank){ //neighbor or me

            // cout << "split " << this->rank << " color " << i << endl << flush;

            MPI_assert(MPI_Comm_split(MPI_COMM_WORLD, i, this->rank, &newcomm), "Comm split error" LOCATION);

            ss << "Comm " << i;
            MPI_assert( MPI_Comm_set_name(newcomm, ss.str().c_str()), "Set comm name err" LOCATION);
            ss.str("");

            if(newcomm == MPI_COMM_NULL){
                stringstream ss;
                ss << "initComms: newcomm is NULL, possible fail" << endl;
                throw runtime_error(ss.str());
            }

            int size;
            MPI_assert( MPI_Comm_size(newcomm, &size) LOCATION);
            scatCnts = new int[size];
            scatDisp = new int[size];

            if(i == this->rank){

                myComm = newcomm;
                MPI_assert( MPI_Comm_rank(myComm, &myCommRank), "Rank in my comm failed" LOCATION);
                // MPI_assert( MPI_Comm_group(newcomm, (&tmpGroup)), "tmpGroup failed"  LOCATION);

                for(unsigned i = 0 ; i < counts.size();i++){
                    scatCnts[i] = counts[i];
                    scatDisp[i] = displs[i];
                }

                MPI_assert( MPI_Bcast(scatCnts, size, MPI_INT, myCommRank, myComm) LOCATION);
                MPI_assert( MPI_Bcast(scatDisp, size, MPI_INT, myCommRank, myComm) LOCATION);

                delete[] scatCnts;
                delete[] scatDisp;

               

            }else{ // i is neighbor

                int tmp = 0, root = 0;
                MPI_assert( MPI_Comm_rank(newcomm, &tmp), "Rank in newcomm failed" LOCATION);

                nData.at(i).comm = newcomm;
                nData.at(i).wRank = i;
                nData.at(i).myRank = tmp;

                MPI_assert( MPI_Comm_group(newcomm, (&tmpGroup)), "tmpGroup failed"  LOCATION);
                MPI_assert( MPI_Group_translate_ranks(world, 1, &i, tmpGroup, &root ), "group rank translate " LOCATION);

                if(root == MPI_UNDEFINED)
                    throw runtime_error("rank translation failed");
                // since I don't know it's rank, I have to use ANY_SOURCE
                nData.at(i).root = root;


                MPI_assert( MPI_Bcast(scatCnts, size, MPI_INT, root, newcomm) LOCATION);
                MPI_assert( MPI_Bcast(scatDisp, size, MPI_INT, root, newcomm) LOCATION);

                nData.at(i).scatterCnts = scatCnts;
                nData.at(i).scatterDispls = scatDisp;

            }

        }else{

            // do not participate in this communicator
            // I'm not actually on turn and i is no my neighbour
            // will receive MPI_COMM_NULL
            MPI_assert( MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, this->rank, &newcomm), "Comm split error" LOCATION);

            if(newcomm != MPI_COMM_NULL){
                throw runtime_error("initComms: newcomm not NULL, possible fail");
            }
        } //else




        MPI_Barrier(MPI_COMM_WORLD);

    } //for
}


string TopologyDescriptor::toString(void)
{
    stringstream ss;
    ss << COUTLOC << endl << "DBD for rank " << rank << endl;
    
    for(auto t : tiles){
        ss << t.toString() << endl;
        ss << "=================" << endl;
    }

    ss << "Rank " << rank << " scatter:";
    ss << "Displacements: ";
    
    for(auto o : displs){
        ss << o << " ";
    }

    ss << endl ;    
    ss << "Neighbors: ";

    for(auto n : neighbors){
        ss << n << " ";
    }

    ss << endl ;    
    ss << "Counts: ";

    for(auto c : counts){
        ss << c << " ";
    }

    ss << endl ;    

    return ss.str();

}

string TopologyDescriptor::commsToString(void)
{
    stringstream ss;


    ss << COUTLOC <<  "Rank " << rank << endl;

    char  * name = new char[MPI_MAX_OBJECT_NAME];
    int len;
    int size;

    MPI_assert( MPI_Comm_get_name(myComm, name, &len), "Get name failed" LOCATION);
    MPI_assert( MPI_Comm_size(myComm, &size) LOCATION);

    ss << "myComm: name = " << name <<" rank = "  << myCommRank << " size " << size <<  endl;
    ss << "=======" << endl << "neighbors: " ;

    for(auto n: neighbors) ss << n << " ";
    ss << endl;

    for(auto n : nData){

        ss << "Neighbor: " << n.second.wRank << endl;;

        MPI_assert( MPI_Comm_get_name(n.second.comm, name, &len), "Get name failed" LOCATION );
        ss << name << " ";

        ss << "myRank in " << name << " " << n.second.myRank << endl;
        ss << "root: "    << n.second.root << endl;
        ss << "count: "  << n.second.count << endl;
        ss << "displ: "  << n.second.displ << endl;

    }

    ss << "================================" << endl;

    return ss.str();

}


void TopologyDescriptor::setBorders(void)
{
	Dims p = myTile.getPosition();
    Dims s = myTile.getSize();

    for(int i = 0; i < 4;i++){
    	borders[i] = false;
    }

    if(p.x == 0){
    	borders[TEdge::LEFT] = true;
        // bd.left = 4;
    }
    if(p.y == 0){
    	borders[TEdge::TOP] = true;
        // bd.top = 4;
    }
    if(p.x + s.x == edgeSize ){
    	borders[TEdge::RIGHT] = true;
    }
    if(p.y + s.y == edgeSize){
        // bottom
        borders[TEdge::BOTTOM] = true;
        // bd.bottom = s.y;
    }

}



void TopologyDescriptor::getBorders(BlockData & bd)
{
    Dims s = myTile.getSize();

	bd.top = 2;
    bd.right = s.x + 2;
    bd.bottom = s.y + 2;
    bd.left = 2;

    // set initial
    
  	if(borders[TEdge::TOP]){
  		bd.top = 4;
   	}
   	if(borders[TEdge::RIGHT]){
        bd.right = s.x;
   	}
   	if(borders[TEdge::BOTTOM]){
        bd.bottom = s.y;
   	}
   	if(borders[TEdge::LEFT]){
        bd.left = 4;
   	}

    // flag for computing halos
    bd.topF = borders[TEdge::TOP];
    bd.rightF = borders[TEdge::RIGHT];
    bd.bottomF = borders[TEdge::BOTTOM];
    bd.leftF = borders[TEdge::LEFT];


}

DLB::BlockData TopologyDescriptor::getBlockData(void)
{
    BlockData bd;

    setBorders();
    getBorders(bd);


    // neighbor data
    bd.nData = &nData;

    bd.displs = &displs;
    bd.counts = &counts;
    bd.neighbors = &neighbors;

    bd.myComm = myComm;
    bd.myCommRank = myCommRank;
    

    //middle communicator
    bd.middle = middle;
    bd.COMM_MIDDLE = COMM_MIDDLE;
    bd.midRank = midRank;
    bd.midSize = midSize;

    return bd;
}

void TopologyDescriptor::setTile(const TileDescriptor & t)
{
    myTile = t;
}


void TopologyDescriptor::midUpdate(void)
{

    // is actual tile in the middle ?

    middle = myTile.isMiddle( edgeSize / 2);

    int color = middle ? 1 : MPI_UNDEFINED;

    MPI_assert(MPI_Comm_split(MPI_COMM_WORLD, color, rank, &COMM_MIDDLE ), 
        "Comm_split: error creating middle col communicator" LOCATION );

    if(middle){

        MPI_assert( MPI_Comm_size(COMM_MIDDLE, &midSize) LOCATION );
        MPI_assert( MPI_Comm_rank(COMM_MIDDLE, &midRank) LOCATION ) ;

    }else {

        midSize = -1;
        midRank = -1;
    }

}