# /**
# * @File        Makefile
# * @Author      Jiri Jaros, Filip Vaverka
# * @Affiliation FIT, Brno University of Technology
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Salomon
# * 
# * @Tool        ARC project 2015
# * @Created     10 April  2015, 10:49 AM
# * @LastModif   05 April  2016, 02:37 PM
#
# */


# include /usr/local/include/Makefile.export.zoltan
# include /usr/local/include/Makefile.export.zoltan.macros


#Load following modules first for icpc
#   module load intel/2016.01
#   module load HDF5/1.8.13-intel-2016.01


#SET PATHS
#HDF5_DIR=/usr/local/hdf5-serial

# HDF5_DIR=/apps/libs/hdf5/intel/1.8.13_parallel
# HDF5_MIC_DIR=../hdf5-mic/install
HDF5_DIR=/usr/lib64
ZOLTAN_LIB=~/lib/libzoltan_intel.a  

SRC=../Sources
SRCDLB=$(SRC)/DLB
#HDF5_MIC_DIR=../hdf5-mic/install


#SET COMPILATOR, FLAGS and LIBS
CXX = icpc
MPICXX = mpicxx
SCOREP_CXX = scorep --nouser mpicxx

CXXFLAGS=-W -Wall -Wextra -pedantic \
         -O0 -std=c++11 -ggdb3 \
	     -I../Sources/ \
	     -I. \
	     -I ../Sources/DLB \
   	     -I ../Sources/DLB/Logger/ \
	     -I$(HOME)/lib/trilinos/include \
	     -I$(HOME)/lib/include 


CXXFLAGS_NOMIC=-march=native -I$(HDF5_DIR)/include
#CXXFLAGS_NOMIC=-I$(HDF5_DIR)/include
CXXFLAGS_MIC=-mmic -I$(HDF5_MIC_DIR)/include 

	     
LDFLAGS= -std=c++11 -L$(SRCDLB)/Logger -L$(SRCDLB) -L. -L$(HOME)/lib/trilinos/lib/ -lzoltan -lhdf5 -L$(HDF5_DIR)/lib/ -Wl,-rpath,$(HDF5_DIR)/lib/
LDBOOST= -lboost_unit_test_framework -lboost_system

# LDFLAGS_NOMIC=-march=native -L$(HDF5_DIR)/lib/ -Wl,-rpath,$(HDF5_DIR)/lib/
#LDFLAGS_NOMIC=-L$(HDF5_DIR)/lib/ -Wl,-rpath,$(HDF5_DIR)/lib/
LDFLAGS_MIC=-mmic -L$(HDF5_MIC_DIR)/lib/ -Wl,-rpath,$(HDF5_MIC_DIR)/lib/

LIBS=-lhdf5

DEPS= $(SRC)/MaterialProperties.o $(SRC)/BasicRoutines.o  $(SRCDLB)/Logger/Logger.o \
	  $(SRCDLB)/DynamicBlockDescriptor.o $(SRCDLB)/LoadBalancer.o $(SRCDLB)/PerfMeasure.o $(SRCDLB)/TileDescriptor.o $(SRCDLB)/TopologyDescriptor.o $(SRCDLB)/Dims.o \
	  $(SRCDLB)/TileMsg.h $(SRCDLB)/BlockData.h $(SRCDLB)/Asserts.h

TARGET=PerfMeasureTestbench
TARGET_MIC=arc_proj02_mic


#zoltan macros
OBJS := $(patsubst %.cpp,%.o,$(wildcard ../Sources/DLB/*.cpp))
# OBJS := $(wildcard ../Sources/DLB/*.cpp))

.PHONY: build link 

# all:	$(TARGET)


#push:
#	scp proj02.cpp BlockDescriptor.cpp BlockDescriptor.h Makefile xdvora0y@salomon.it4i.cz:/home/xdvora0y/ARC2/Sources/


# all: TileDescriptorUnit TileDescriptorTestbench DynBlockDescriptorTestbench zoltan ListTestbench
all: TileDescriptorTestbench

TileDescriptorTestbench: TileDescriptorTestbench.cpp 
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o TileDescriptorTestbench TileDescriptorTestbench.cpp $(DEPS)

TileDescriptorUnit: Unittests/TileDescriptorUnit.cpp TileDescriptor.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDBOOST) -o Unittests/TileDescriptorUnit $(SRC)DLB/TileDescriptor.o Unittests/TileDescriptorUnit.cpp

TileDescriptor.o: $(SRC)DLB/TileDescriptor.cpp
	$(CXX) $(CXXFLAGS) -c -o $(SRC)DLB/TileDescriptor.o $(SRC)DLB/TileDescriptor.cpp

ListTestbench: ListTestbench.cpp
	$(CXX) $(CXXFLAGS) -o ListTestbench	ListTestbench.cpp

zoltan: ZoltanTestbench.cpp
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o ZoltanTestbench ZoltanTestbench.cpp 

DynBlockDescriptorTestbench: DynBlockDescriptorTestbench.cpp
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(OBJS) -o DynBlockDescriptorTestbench DynBlockDescriptorTestbench.cpp

LoadBalancerTest: LoadBalancerTest.cpp $(SRCDLB)/LoadBalancer.o
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o LoadBalancerTestLoadBalancerTest.cpp $(DEPS)



zoltan_heap: zoltan_heap.cpp
	$(MPICXX) $(CXXFLAGS) $(LDFLAG S) $(LIBS) -o ZoltanHeap zoltan_heap.cpp $(ZOLTAN_LIB)


clean:
	rm -f TileDescriptorTestbench
	rm -f TileDescriptorUnit
