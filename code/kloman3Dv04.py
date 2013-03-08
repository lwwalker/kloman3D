#!/usr/bin/python
##kloman3D -- kernel learning of monoclonal antibody neutralization in 3D

import os
import sys 
import getopt
sys.path.append(".")
from kern3Dv04 import *
import time
import copy

"""
kloman3Dv04.py
Depends on kern3Dv04 or higher
Changes from version 03:
Simplified selection scheme, mostly a simple version for debugging
"""

def usage():
	print "\nUsage: kloman3D.py [optional flags] fastafile gp120file responsefile outfile"
	print "Manditory parameters:"
	print "fastafile -> filename of fasta file, expected that first sequence will be consensus which matches the pdb protein structure data(consensus only sequence read)" 
	print "gp120file -> filename for .pdb file describing structure of gp120"
	print "responsefile -> filename of file containing training classes, should be 2 column comma-delimited file with column 0 being sample names matching fasta file, and column 1 being 1 or -1 indicating classes" 
	print "outfile -> file to write output to\n"
	print "Optional parameters:"
	print "-m --map -> filename of seperate alignment map of sequences to consensus"
	print "-M --matrix -> filename of pre-computed array of kernel matrices"
	print "-g --gp41 -> optional pdb filename for gp41"
	print "-N --norm -> normalize the kernel matrixes"
	print "-r --radius -> The \"bubble\" radius for the kernel (default 10 Angstroms)"
	print "-s --stdev -> Standard deviation for the word kernel function (default 10 Angstroms)"
	print "-x --nocross -> use non-cross-validation SVMs for feature selection"
	print "-a --aminoacid -> filename for an amino acid partial-match matrix"
	print "-c --cSVM -> perform c-svm with the passed regularization parameter"
	print "-n --nuSVM -> perform nu-svm with the passed regularization parameter (default with nu = 0.2)"
	print "-v --saveKM -> Save the list of position-wise kernel matrixes to supplied filename for future use"
	print "-h --help -> Display this message\n"

def main(argv):
	timeStart = time.time()
	argString = "kloman3D.py "+" ".join(argv)+" \nCalled at " + time.strftime("%a, %d %b %Y %H:%M:%S") + "\n"
	
	#Initialize defaults
	myRad = 10
	mySigma = 10
	svmType = "nu"
	myParam = 0.2
	myIter = 0
	mapFileName = None
	AAmap = None
	normlze = False
	loadKM = None
	gp41filename = None
	pdbfilename = None
	nocross = False
	saveKMfilename = None

	#Read in command line arguments
	try:                                
		opts, args = getopt.getopt(argv, "hNxr:s:a:m:M:g:c:n:v:", ["help","norm","nocross","radius","stdev","aminoacid","map","matrix","gp120","gp41","cSVM","nuSVM","saveKM"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-x","--nocross"):
			nocross = True
		if opt in ("-m","--map"):
			try:
				mapFileName = arg
				mapFile = open(mapFileName,'r')
				mapFile.close()
			except IOError:
				print "\nCould not find map file "+arg
				sys.exit(2)
		if opt in ("-M","--matrix"):
			try:
				loadKM = loadKmats(arg)
			except ValueError:
				print "Invalid format in Kernel Matrix Load"
				sys.exit(2)	
		if opt in ("-g","--gp41"):
			try:
				open(arg,'r')
				gp41filename = arg
			except IOError:
				print "\nInvalid gp41 pdb file"
				sys.exit(2)
			if not mapFileName is None:
				print "\nPassed both a structure map and structure files.  One or the other, please."				
		if opt in ("-N","--norm"):
			if not loadKM is None:
				print "Passed pre-computed kernel AND kernel calculation parameters: {0}".format("\"Normalize\"")
				sys.exit(2)
			normlze = True
		if opt in ("-a","--aminoacid"):
			try:
				AAmap = parseAAmap(arg)
			except IOError:
				print"\nCould not find amino-acid partial match file " +arg
				sys.exit(2) 
		if opt in ("-s","--stdev"):
			if not loadKM is None:
				print "Passed pre-computed kernel AND kernel calculation parameters: {0}".format("\"Stdev\"")
				sys.exit(2)		
			try: 
				mySigma = float(arg)
			except ValueError:
				print "\nInvalid kernel sigma input."
				sys.exit(2)
		if opt in ("-r","--radius"):
			if not loadKM is None:
				print "Passed pre-computed kernel AND kernel calculation parameters: {0}".format("\"Radius\"")
				sys.exit(2)		
			try:
				myRad = float(arg)
			except ValueError:
				print "\nInvalid kernel radius input."
				sys.exit(2)
		if opt in ("-c","--cSVM"):				
			try:
				svmType = "C"
				myParam = float(arg)
			except ValueError:
				print "\nInvalid c-svm parameter."	
				sys.exit(2)
		if opt in ("-n", "--nuSVM"):
			try:
				svmType = "nu"
				myParam = float(arg)
			except ValueError:
				print "\nInvalid nu-svm parameter."
				sys.exit(2)
		if opt in ("-v", "--saveKM"):
			try:
				saveKMfilename = arg
			except ValueError:
				print "\nInvalid kernel matrix save filename"
				sys.exit(2)
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
	if len(args) < 4:
		print "\nMissing filename argument, expected: fastafile gp120file responsefile outputfile"
		usage()
		sys.exit(2)

	fastafilename = args.pop(0)
	pdbfilename = args.pop(0)	
	responsefilename = args.pop(0)
	outputfilename = args.pop(0)

	print("Initalizing . . .")
	allResp = parseY(responsefilename)	#Load in titers
	outfile = open(outputfilename, 'w')
	sys.stdout = outfile
	print argString
 	outfile.flush()
 	os.fsync(outfile.fileno())

	#Load kernel matrixes if supplied
	if not loadKM is None:
		print("Loadking Kernel Matrixes")
		kernMats = loadKM
	else:
		print("Calculate kernel matrixes, ")
		print("Radius:"+str(myRad)+", SD:"+str(mySigma)+"\n")
	 	outfile.flush()
	 	os.fsync(outfile.fileno())
	 	
		#Calculate kernel matrixes
		kernMats = reCalcPOS(myRad, mySigma, pdbfilename, fastafilename, mapFileName, AAmap, normlze, gp41filename)

	
	n = len(kernMats)
	kMatDim = kernMats[0].dim

	#Associate titer data with samples in kernel matrixes
	y = []
	for name in kernMats[0].sampleNames:
		try:
			y.append(int(allResp[name]))
		except KeyError:
			print "Unable to find response value for sequence \"{}\"".format(name)
			sys.exit(2)
	
	#Save kernel matrixes if so instructed
	if not saveKMfilename is None:
		saveKmats(kernMats, saveKMfilename) 

	t1 = time.time()
	dirModel = simpleBidirectionalSelection(kernMats, myParam, svmType, y)
	t2 = time.time()

	print "Tuning regularization parameters"
	tuneModel = regularizationTune(dirModel)

	tuneModel.saveModel(outputfilename.split('.')[0]+".model")
	print "Model calculated in {} seconds".format(round(t2-t1,2))	
	print tuneModel.prn()

if __name__ == "__main__":
	main(sys.argv[1:])
