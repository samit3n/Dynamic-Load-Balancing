Dynamic Load Balancing
======================

This project implements dynamic load balancing mechanism into
the parallel heat distribution model of a CPU cooler.

This model uses MPI message passing library and HDF5 library for data I/O.

Original model was designed as one of student projects as part of Parallel System Architecture and Programming
course at Faculty of Information Technology, VUT Brno.

This project implements ability for dynamic load balancing, which as runtime optimization
feature, redistributing work among processes based on current load.


Repository strucutre:

 * doc 

 	*	doxygen documentation
 	* 	related pdfs

 * src

 	* 	Sources - project sources
 	*	Scripts - experiments
 	* 	DataGenerator - input data
 	* 	Tests - unittest, testbenches

 * thesis

 	* 	thesis pdf
 	*	thesis sources



Build & Run
-----------

Project may be compiled using Makefile in sources directory.


Dependencies
------------

 * C++11 compatible compiler with mpi support
 * MPI version 2
 * HDF5 I/O library		(https://support.hdfgroup.org/HDF5/)
 * Zoltan library 		(https://trilinos.org/packages/zoltan/)
















