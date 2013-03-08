#!/usr/bin/python

import sys

#makeResp.py - parse antibody data into binary response vector
#usage: ./makeResp.py titerFile map/fastaFile outputFile
#titerFile - should be a tab-delmited file with columns named "ID" for sample ID and "IC50" for titer data
#map/fastaFile - should be the HyPhy map or fasta alignment file to be used for prediction
#outputFile - a two-column file with ID and "-1" for sensitive (titer <= 25) and "1" for resistant (titer > 25)

arg = sys.argv
if len(arg) != 4:
	print "\nError: Wrong number of arguments, expected 3, found {}\n".format(len(arg)-1)
	print "Usage: ./makeResp.py titerFile map/fastaFile outputFile \ntiterFile - should be a tab-delmited file with columns named \"ID\" for sample ID and \"IC50\" for titer data \nmap/fastaFile - should be the HyPhy map or fasta alignment file to be used for prediction \noutputFile - a single-column file with \"1\" for sensitive (titer <= 25) and \"-1\" for resistant (titer > 25) \n"
	sys.exit(1)

#Set up files
titerFileN , mapFileN , outFileN = arg[1:4]
try:
	titerFile = open(titerFileN ,'r')
except IOError:
	print "Error: Unable to	open titer file\n"
	sys.exit(1)
	
try:
	mapFile = open(mapFileN ,'r')
except IOError:
	print "Error: Unable to	open map/fasta file\n"
	sys.exit(1)
mapFile.close()

try:
	outFile = open(outFileN, 'w')
except IOError:
	print "Error: Unable to	open map/fasta file\n"
	sys.exit(1)
	
#Read in titer data
titerLines = titerFile.readlines()
titerFile.close()

headerLine = titerLines.pop(0).strip("\n").split("\t")
try:
	idCol = headerLine.index("ID")
	titerCol = headerLine.index("IC50")
except ValueError:
	print "Error: unable to parse titer file - must have colunms named \"ID\" and \"IC50\"."
	sys.exit(1)

#Extract and categorize titer data
titers = {}
numSens = 0
for line in titerLines:
	line = line.strip("\n").split("\t")
	titerID = line[idCol].strip()
	titerVal = line[titerCol].strip()
	try:
		if titerVal == ">50" or titerVal == ">25":
			titers[titerID] = 1
		elif titerVal == ">12.5":
			print "Titer Value of >12.5 found for {}, scoring as resistant".format(titerID)
			titers[titerID] = 1
		elif float(titerVal) > 25:
			titers[titerID] = 1
		else:
			titers[titerID] = -1
			numSens += 1
	except ValueError:
		print "Error: non-numeric, non \">50\" or \">25\" value found in titer values, found {}".format(titerVal)
		sys.exit(1)
		
print "found {} titer values, of which {} are sensitive".format(len(titers.keys()), numSens)

#Load .py map file or .fasta/.fas sequence file.  

if mapFileN.find(".py") >= 0:
	d = {}
	execfile(mapFileN,d,d)
	sequenceMaps = d["sequenceMaps"]
	seqIDs = sequenceMaps.keys()
elif mapFileN.find(".fas") >= 0:
	mapFile = open(mapFileN, 'r')
	mapLines = mapFile.readlines()
	mapFile.close()
	seqIDs = {}
	for line in mapLines:
		line = line.strip("\n")
		if line.find(">") >= 0:
			line = line.strip(">")
			seqIDs[line] = 1
	seqIDs = seqIDs.keys()
	
else:
	print "Error, could not identify map/fasta file type.  Must have file extension .py or .fas(ta)\n"
	sys.exit(1)
		
#Output data
outFile.write("ID,sens\n")
for myID in seqIDs:
	try:
		myVal = titers[myID]
		outFile.write("{},{}\n".format(myID,myVal))
	except KeyError:
		print "Error, no titer in titer file for squence ID {}.".format(myID)
		
outFile.close()

