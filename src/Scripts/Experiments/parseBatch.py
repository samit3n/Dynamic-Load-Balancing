#!/usr/bin/env python3

# parsing output files of balancing experiments
# creates plots using matplotlib

import sys
from operator import itemgetter
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from collections import OrderedDict
import re
from itertools import product


matplotlib.rc('font', family='DejaVu Sans')

multipliers=  [0.5, 0.75, 1, 1.5]
cores = [16,64,128,256]

name = sys.argv[1][:sys.argv[1].rfind('.')]
outfname = name + ".txt"
outfplot = name + ".png"

inf = open(sys.argv[1])
outf = open(outfname, 'w')


parsed = list() 	# list of parsed entries
ldata = dict() 		# dict represents single entry

first = True
for line in inf:
	# print(line)

	#if first: #skip first line
	#	first = False
	#	continue

	# measurement delimiter reached
	if "--" in line:
		parsed.append(ldata) #store ldata to list
		ldata = dict() # clear ldata
	else: # store measured data
		# print(line)
		line.replace('\n','')
		# print(line)
		tmp = line.split(':')

		# print(tmp)
		ldata[tmp[0]] = tmp[1][:-1]

# print(parsed)

#declare result dict
times = dict()
for key in product(cores, ["par", "parBal"]):
	times[key] = list()

# pcent = map(lambda m: m*25, multipliers)

for x in parsed:
	# parse multiplier value
	m = re.search(r"-M\s+(\d+\.\d+|\d+)", x["cmd"])
	multiplier = float( m.group(1) )
	# percent = (float(x["SleepTotal[ms]"]) / 1000) / float(x["TotalTime"]) * 100
	# if int(x["MPIProcs"]) != 16:
	times[ ( int(x["MPIProcs"]), x["Mode"] ) ].append( (float(x["TotalTime"]), multiplier * 25) )

	# print(float(x["BalanceTotal"]) / float(x["TotalTime"]) * 100)


for k,v in times.items():
	v.sort(key=itemgetter(1))
	print(str(k) + ":" + str(v) )
	outf.write(str(k) + ":" + ";".join([str(i[0]) for i in v]) + "\n")


print("Percent")
for p in cores:
	par = [ t[0] for t in times[(p, "par")]]
	parBal = [ t[0] for t in times[(p, "parBal")]]
	res =[ (p - pb) for p,pb in zip(par,parBal)]

	print(p)
	print(par)
	print(parBal)
	print(res)
	for x,r in zip(par, res):
		print(r / x * 100)






# print(times)
#for k,v in times.items():
#	tmpmult = [m[1] for m in v]
#	for mul in multipliers:
#		if mul not in tmpmult:
#			times[k].append((0.0, mul))






# plot graph
plt.xlabel(u"Podíl zpoždění [%]")
plt.ylabel(u"Délka výpočtu [s]")
plt.title(u"Délka výpočtu s/bez vyvážení")

lines = list()

colors = ['b', 'g','k','c','r', 'y', 'm']
cnt = 0
col = 0

for k in sorted(times.keys(), key=lambda x: x[0]):
	# print(k)
	print(times[k])
	tms = [x[0] for x in times[k]]
	x = np.array([p[1] for p in times[k]])
	# print(tms)

	lab = " s vyv." if k[1] == "parBal" else " bez vyv."
	if lab == " s vyv.":
		ln, = plt.plot(x,tms, label=str(k[0])+ lab, linewidth=2, ls="solid", color=colors[col])
	else:
		ln, = plt.plot(x,tms, label=str(k[0])+ lab, linewidth=2, ls="dashed",  color=colors[col])

	cnt += 1
	col = (col + 1) if cnt % 2 == 0 else col

	lines.append(ln)

plt.legend(handles=lines, loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()

plt.savefig(outfplot, bbox_inches='tight', dpi=300)
# plt.show()




outf.close()
inf.close()


