#!/usr/bin/env python3

# returns maximal size of object eg. width of single block
# in regular domain decomposition
# used by experiment scripts

#arg1: number of processes
#arg2: size of domain


from math import *
import sys

procs = int(sys.argv[1])
domsize = int(sys.argv[2])

log = int(log(procs, 2))

objsize = 0

if log % 2 == 0:
	objsize = int(domsize / sqrt(procs))
else:
	objsize = int(domsize / sqrt(procs*2))

print(objsize)
	







