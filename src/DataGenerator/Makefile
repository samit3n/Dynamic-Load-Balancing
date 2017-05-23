# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# * @Affiliation FIT, Brno University of Technology
# * @Email       jarosjir@fit.vutbr.cz
# * @Comments    Linux makefile for Salomon
# * 
# * @Tool        ARC project 2015
# * @Created     15 February  2015, 15:25 AM
# * @LastModif   05 April     2016, 02:37 PM
#
# */


#Load following modules first for GCC
#   module load gcc/4.9.0
#   module load hdf5/1.8.13



#Load following modules first for icpc
#   module load intel/2016.01
#   module load HDF5/1.8.13-intel-2016.01

HDF5_DIR=/apps/all/HDF5/1.8.16-intel-2017.00

# HDF5_DIR=/apps/libs/hdf5/intel/1.8.13_serial
#DF5_DIR=/usr/lib64/openmpi/bin
# HDF5_DIR=/usr/lib64/

CXX		= icpc
#CXX		= g++

CXXFLAGS        = -std=c++11 -O3 -openmp  -I$(HDF5_DIR)/include -I.

TARGET		= arc_generator
LDFLAGS		= -std=c++11 -O3 -openmp  -L$(HDF5_DIR)/lib \
        -Wl,-rpath,$(HDF5_DIR)/lib


all:		$(TARGET)	


$(TARGET):	main.o 					

#---------------------------------------------------------

	$(CXX) $(LDFLAGS) main.o 			\
		-lhdf5_hl -lhdf5 \
		-o $@


$(TARGET).o : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp

test:
	./arc_generator -o material.h5 -N 128 -H 100 -C 20


clean:
	rm -f *.o  $(TARGET)

