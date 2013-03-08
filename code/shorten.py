#!/usr/bin/python
import sys
n = 25

d = {}
execfile("../rawData/2F5short.py",d,d)
sequenceMaps = d["sequenceMaps"]

ks = sequenceMaps.keys()

smOut = {}
for k in ks:
	smOut[k] = {}
	for i in range(25):
		smOut[k][str(i)] = sequenceMaps[k][str(i)]

outfile = open("../rawData/2F5vshort.py",'w')
outfile.write('sequenceMaps = {}\n')
for k in smOut.keys():
	outfile.write("sequenceMaps ['{}'] = ".format(k),)
	outfile.write(str(smOut[k])+"\n")
	
