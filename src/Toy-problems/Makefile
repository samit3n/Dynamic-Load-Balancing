CC=mpic++
CFLAGS=-std=c++11

.PHONY: hellow circ_block  circ_noblock circ_noblock vec_prod wall_time

all: hellow circ_block circ_noblock vec_prod

hellow: mpi_hworld.c
	$(CC) $(CFLAGS) mpi_hworld.c -o hworld

circ_block: circle_block.cpp
	$(CC) $(CFLAGS) circle_block.cpp	-o circle_block

circ_noblock: circle_noblock.cpp
	$(CC) $(CFLAGS) circle_noblock.cpp -o circle_noblock

vec_prod: vector_product.cpp
	$(CC) $(CFLAGS) vector_product.cpp -o vector_product

wall_time: wall_time.cpp
	$(CC) $(CFLAGS) wall_time.cpp -o wall_time




