all: Logger-testbench.o Logger.o
	g++ -std=c++11 -I../ Logger-testbench.cpp Logger.o -o testbench

%.o: %.cpp
	g++ -std=c++11 -I../ -c -o $@ $<
