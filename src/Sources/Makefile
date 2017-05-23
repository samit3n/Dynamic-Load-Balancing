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


#Load following modules first for icpc
#   module load intel/2016.01
#   module load HDF5/1.8.13-intel-2016.01


#SET PATHS
# HDF5_DIR=/usr/local/hdf5-serial

#HDF5_DIR=/apps/all/HDF5/1.8.16-intel-2017.00

HDF5_DIR=/apps/all/HDF5/1.10.0-patch1-intel-2017a
# HDF5_MIC_DIR=../hdf5-mic/install

# HDF5_DIR=/usr/lib64

#zoltan static linkage
# ZOLTAN_LIB=~/lib/libzoltan_intel.a  
# ZOLTAN_LIB=~/lib/libzoltan_ompi_gcc.a
# ZOLTAN_LIB=~/lib/trilinos/lib/libzoltan.a

#HDF5_MIC_DIR=../hdf5-mic/install

#SET COMPILATOR, FLAGS and LIBS
CXX = mpicxx

SCOREP_CXX = scorep --nouser mpicxx

CXXFLAGS=-W -Wall -Wextra -pedantic  \
         -O3 \
         -std=c++11 -g \
	     -DPARALLEL_IO \
	     -I. \
	     -IDLB \
	     -IDLB/Logger \
	     -I$(HDF5_DIR)/include \
	     -I$(HOME)/lib/trilinos/include \
	     -xhost

		# -I/apps/tools/trilinos/11.2.3/include
	     # -march=native \

		 # -fdiagnostics-color=always 

LDFLAGS=-std=c++11 -O3 \
		 -LLogger -LDLB -L.  \
		 -L$(HOME)/lib/trilinos/lib/ \
		 -lzoltan -lhdf5 -L$(HDF5_DIR)/lib/ -Wl,-rpath,$(HDF5_DIR)/lib/ \
	     -xhost


		 #-L/apps/tools/trilinos/11.2.3/lib
	

# original

CXXFLAGS_NOMIC=-march=native -I$(HDF5_DIR)/include
#CXXFLAGS_NOMIC=-I$(HDF5_DIR)/include
CXXFLAGS_MIC=-mmic -I$(HDF5_MIC_DIR)/include

LDFLAGS_NOMIC=-march=native -L$(HDF5_DIR)/lib/ -Wl,-rpath,$(HDF5_DIR)/lib/
#LDFLAGS_NOMIC=-L$(HDF5_DIR)/lib/ -Wl,-rpath,$(HDF5_DIR)/lib/
LDFLAGS_MIC=-mmic -L$(HDF5_MIC_DIR)/lib/ -Wl,-rpath,$(HDF5_MIC_DIR)/lib/

DEPS= dlb_heat.o MaterialProperties.o BasicRoutines.o  DLB/Logger/Logger.o \
	  DLB/DynamicBlockDescriptor.o DLB/LoadBalancer.o DLB/PerfMeasure.h DLB/TileDescriptor.o DLB/TopologyDescriptor.o DLB/Dims.o \
	  DLB/TileMsg.h DLB/BlockData.h DLB/Asserts.h DLB/Neighbor.h DLB/HaloBuffers.h

TARGET=arc_proj02
TARGET_MIC=arc_proj02_mic

NP=16
WI=10
ITER=100
MATER=8192
T=50
M=0.2
OBJ=8

.PHONY: build link

all:	$(TARGET)


#push:
#	scp proj02.cpp BlockDescriptor.cpp BlockDescriptor.h Makefile xdvora0y@salomon.it4i.cz:/home/xdvora0y/ARC2/Sources/

push:
	scp proj02.cpp \
		BlockDescriptor.cpp \
		BlockDescriptor.h \
		Makefile \
		PerfMeasure.h \
		PerfMeasure.cpp \
		xdvora0y@salomon.it4i.cz:/home/xdvora0y/ARC2/Sources/

push2:
	scp proj02.cpp BlockDescriptor.cpp BlockDescriptor.h Makefile xdvora0y@salomon.it4i.cz:/scratch/work/user/xdvora0y/ARC2/Sources/


$(TARGET): $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXXFLAGS_NOMIC) -c -o $@ $< 


#DLB/DynamicBlockDescriptor.o: DLB/DynamicBlockDescriptor.cpp
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXXFLAGS_NOMIC) -c -o $@ $<

DLB/%.o: DLB/%.cpp DLB/%.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXXFLAGS_NOMIC) -c -o $@ $<

DLB/Logger/%.o: DLB/Logger/%.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $<

$(TARGET_MIC): proj02_mic.o MaterialProperties_mic.o BasicRoutines_mic.o BlockDescriptor.o
	$(CXX) $(LDFLAGS) $(LDFLAGS_MIC) $^ $(LIBS) -o $@

%_mic.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_MIC) -c -o $@ $<

test:
	./arc_proj02 -n 1000 -m 0 -w 5 -i ../DataGenerator/material256_100_20.h5 -o result256.h5

testmpi: all $(DEPS)
	mpirun -np $(NP) ./arc_proj02 -n $(ITER)  -v -p -m 1 -w $(WI) -i $(INDATA)/input_data_$(MATER).h5 -o result$(MATER).h5

testdlb: all $(DEPS)
	mpirun -np $(NP) ./arc_proj02 -b -n $(ITER) -X -t $(T) -p -m 1 -w $(WI) -M ${M} -s ${OBJ} -i $(INDATA)/input_data_$(MATER).h5 -o result$(MATER)_X.h5

testmust: all
	mustrun -np 4  --must:nodesize 2 ./arc_proj02 -n 100 -m 1 -w 10 -i ../DataGenerator/material32.h5 -o result.h5

#profile:
#	scorep mpirun -np $(NP) ./arc_proj02 -b -n $(ITER) -X -t $(T) -p -m 1 -w $(WI) -s 8 -i ../Scripts/input_data_$(MATER).h5 -o result$(MATER)_X.h5






testmic:
	./arc_proj02_mic -n 100 -m 0 -w 10 -i ../DataGenerator/material.h5 -o result.h5

clean:
	rm -f *.o
	rm -f DLB/Logger/*.o
	rm -f DLB/*.o
	rm -f *~
	rm -f $(TARGET) $(TARGET_MIC)
