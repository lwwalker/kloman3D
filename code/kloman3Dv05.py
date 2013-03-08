#!/usr/bin/python
##kloman3D -- kernel learning of monoclonal antibody neutralization in 3D

import os
import sys 
import getopt
sys.path.append(".")
from kern3Dv05 import *
import time
import copy

"""
kloman3Dv05.py
Depends on kern3Dv05 or higher
Changes from version 0.4:
-Updated regularization tuning
-Now tunes kernel function radii and standard deviations
-Implements multiple rounds of model selection after tuning
"""

def usage():
	print "\nUsage: kloman3D.py [optional flags] responsefile outfile"
	print "Manditory parameters:"
	print "responsefile -> filename of file containing training classes, should be 2 column comma-delimited file with column 0 being sample names matching fasta file, and column 1 being 1 or -1 indicating classes" 
	print "outfile -> file to write output to\n"
	print "Optional parameters:"
	print "-M --matrix -> filename of pre-computed array of kernel matrices"
	print "-m --map -> filename of seperate alignment map of sequences to consensus"
	print "-f --fasta -> filename for fasta file with reference sequence (if passing map) or all sequences if no map - mandatory of inot loading pre-computed matrixes"
	print "-G --gp120 -> pdb filename for gp120 - mandatory if not loading pre-computed matrixes"
	print "-g --gp41 -> optional pdb filename for gp41"
	print "-N --norm -> normalize the kernel matrixes"
	print "-r --radius -> The \"bubble\" radius for the kernel (default 10 Angstroms), may pass a single value, or comma-delimited list"
	print "-s --stdev -> Standard deviation for the word kernel function (default 10 Angstroms)"
	print "-x --nocross -> use non-cross-validation SVMs for feature selection"
	print "-a --aminoacid -> filename for an amino acid partial-match matrix"
	print "-c --cSVM -> perform c-svm with the passed regularization parameter"
	print "-n --nuSVM -> perform nu-svm with the passed regularization parameter (default with nu = 0.2)"
	print "-v --saveKM -> Save the list of position-wise kernel matrixes to supplied filename for future use"
	print "-h --help -> Display this message\n"

def main(argv):
	timeStart = time.time()
	argString = "kloman3Dv05.py "+" ".join(argv)+" \nCalled at " + time.strftime("%a, %d %b %Y %H:%M:%S") + "\n"
	
	#Initialize variables
	myRads = None
	mySigma = None
	svmType = None
	myParam = None
	mapFileName = None
	AAmap = None
	normlze = False
	loadKMfilename = None
	gp41filename = None
	pdbfilename = None
	nocross = False
	saveKMfilename = None
	pdbfilename = None
	fastafilename = None
	
	#Read in command line arguments
	try:                                
		opts, args = getopt.getopt(argv, "hNxr:s:a:m:M:f:G:g:c:n:v:", ["help","norm","nocross","radius","stdev","aminoacid","map","matrix","fasta","gp120","gp41","cSVM","nuSVM","saveKM"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

	if len(args) < 2:
		print "\nMissing filename argument, expected: responsefile outputfile"
		usage()
		sys.exit(2)

	responsefilename = args.pop(0)
	outputfilename = args.pop(0)
	
	outfile = open(outputfilename, 'w')
	sys.stdout = outfile
	print argString
	print("Initalizing . . .")
	sys.stdout.flush()

	#parse arguments
	for opt, arg in opts:
		#NOCROSS 
		if opt in ("-x","--nocross"):
			nocross = True
			
		#PRECOMPUTED KERNEL MATRIXES
		if opt in ("-M","--matrix"): #using pre-loaded kernel matrixes
			try:
				loadKMfilename = arg
				testFile = open(loadKMfilename, 'r')
				testFile.close()
			except IOError:
				print "\nCould not find kernel matrix file"
				sys.exit(2)				
				
		#PRECOMPUTED COORDINATE MAP
		if opt in ("-m","--map"): #Using a pre-computed coordinate map
			if not(loadKMfilename is None):
				print "\nWarning: pased pre-computed kernel matrixes, map file will be ignored"
			else:
				try:
					mapFileName = arg
					testFile = open(mapFileName,'r')
					testFile.close()
				except IOError:
					print "\nCould not find map file "+arg
					sys.exit(2)
					
		#GP120 PDB FILE
		if opt in ("-G","--GP120"):
			if not(loadKMfilename is None):
				print "\nWarning: pased pre-computed kernel matrixes, GP120 file will be ignored"
			else:
				pdbfilename = arg
				
		#FASTA FILE
		if opt in ("-f", "--fasta"):
			if not(loadKMfilename is None):
				print "\nWarning: pased pre-computed kernel matrixes, fasta file will be ignored"
			else:
				fastafilename = arg
		
		#GP41 PDB FILE
		if opt in ("-g","--gp41"):
			if not(loadKMfilename is None):
				print "\nWarning: pased pre-computed kernel matrixes, gp41 file will be ignored"
			else:
				try:
					gp41filename = arg
					testFile = open(gp41filename, 'r')
					testFile.close()
				except IOError:
					print "\nCould not find gp41 pdb file"
					sys.exit(2)
					
		#NORMALIZE KMs
		if opt in ("-N","--norm"):
			if not(loadKMfilename is None):
				print "\nWarning: passed pre-computed kernel matrixes, KM normalization option will be ignored"
			else:
				normlze = True
				
		#AMINO ACID MAP
		if opt in ("-a","--aminoacid"):
			if not(loadKMfilename is None):
				print "\nWarning: passed pre-computed kernel matrixes, amino acid map will be ignored"
			else:
				try:
					AAmap = parseAAmap(arg)
				except IOError:
					print"\nCould not find amino-acid partial match file " +arg
					sys.exit(2) 
					
		#GAUSSIAN KERNEL SD
		if opt in ("-s","--stdev"):
			if not(loadKMfilename is None):
				print "\nWarning: passed pre-computed kernel matrixes, standard deviation will be ignored"
			else:	
				try: 
					mySigma = float(arg)
				except ValueError:
					print "\nInvalid kernel sigma input."
					sys.exit(2)
		
		#MOLECULAR NEIGHBORHOOD RADIUS
		if opt in ("-r","--radius"):
			if not(loadKMfilename is None):
				print "\nWarning: passed pre-computed kernel matrixes, radius will be ignored"
			else:
				try:
					myRads = map(lambda a: float(a), arg.split(','))
				except ValueError:
					print "\nInvalid kernel radius input."
					sys.exit(2)

		#C-TYPE SVM
		if opt in ("-c","--cSVM"):				
			try:
				svmType = "C"
				myParam = float(arg)
			except ValueError:
				print "\nInvalid c-svm parameter."	
				sys.exit(2)
		
		#NU-TYPE SVM
		if opt in ("-n", "--nuSVM"):
			try:
				svmType = "nu"
				myParam = float(arg)
			except ValueError:
				print "\nInvalid nu-svm parameter."
				sys.exit(2)
				
		#SAVE KM FILENAME
		if opt in ("-v", "--saveKM"):
			if not(loadKMfilename is None):
				print "\nWarning: passed pre-computed kernel matrixes, kernel matrix save filename will be ignored"	
			else:
					saveKMfilename = arg

		#HELP
		if opt in ("-h", "--help"):
			usage()
			sys.exit()

	#Ensure we were passed sufficient information
	if loadKMfilename is None:
		if pdbfilename is None or fastafilename is None:
			print "\nError: Unless pre-computed kernel matrixes are supplied, a GP120 structural file and a fasta file must be given."
			sys.exit(2)
		#Set KM defaults as needed		
		if myRads is None:
			print "Using default radius of 10 angstroms"
			myRads = [10.0]
		if mySigma is None:
			print "Using default sd of 10 angstroms"
			mySigma = 10.0
			
	#Set SVM defaults as needed
	if svmType is None:
		print "Using default nu-SVM with parameter 0.2"
		svmType = "nu"
		myParam = 0.2

	
	#Load response data
	print "Loading response file"
	allResp = parseY(responsefilename)	#Load in titers
	
	#Load kernel matrixes if supplied
	if not loadKMfilename is None:
		print("Loading Kernel Matrixes")
		try:
			kernMats = loadKmats(loadKMfilename)
		except ValueError:
			print "Invalid pre-computed Kernel Matrix file format"
			sys.exit(2)
		referenceSeq=None
		sampleSeqs=None
	else:
		print("Load Sequences . . . ")
		referenceSeq, sampleSeqs = parseSeqData(pdbfilename, fastafilename, mapFileName, allResp.keys(), gp41filename)
	
		print("Calculate kernel matrixes, ")
	 	sys.stdout.flush()
	 	
		#Calculate kernel matrixes
		t0 = time.time()
		kernMats = []
		for aRad in myRads:
			print "Calculating for radus {}".format(aRad)
			kernMats.extend(reCalcPOS(referenceSeq, sampleSeqs, aRad, mySigma, AAmap, normlze))
		t1 = time.time()
		print "Kernel Matrixes calculated in {} seconds".format(round(t1-t0,2))	
	
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

	#Initialize empty model
	workingModel = model(kernMats,[], myParam, svmType, y)
	lastModel = workingModel.copy()


	goFlag = True
	i = 1
	while goFlag:
		print "\nStarting model Selection, cycle {}".format(i)
		t0 = time.time()
		workingModel = simpleBidirectionalSelection(workingModel)
		t1 = time.time()
		print "Model selection performed in {} seconds".format(round(t1-t0,2))	

		print "\nWorking Model:"
		print workingModel.prn()
		sys.stdout.flush()

		print "Tuning regularization parameters, cycle {}".format(i)
		t0 = time.time()
		workingModel = regularizationTune(workingModel)
		t1 = time.time()
		print "Regularization tuning performed in {} seconds".format(round(t1-t0,2))
	
		print "\nWorking Model:"
		print workingModel.prn()
		sys.stdout.flush()

		print "Tuning radius and standard deviation, cycle {}".format(i)
		t0 = time.time()
		workingModel = radiiTune(workingModel, referenceSeq, sampleSeqs, normlze, AAmap)
		t1 = time.time()
		print "Radius/sd tuning performed in {} seconds".format(round(t1-t0,2))
		
		if workingModel.accuracy() <= lastModel.accuracy():
			goFlag = False
			workingModel = lastModel
			print "Cycle {} result ({}) not better than cycle {} ({}), keeping previous result.".format(i,workingModel.acc, i-1, lastModel.acc)
			print workingModel.prn()
		else:
			lastModel = workingModel.copy()
			lastModel.acc = workingModel.acc
			print "Cycle {} result ({}) improved from previous cycle ({}), continue iterating.".format(i,workingModel.acc, lastModel.acc)
			print workingModel.prn()	
		i += 1
		
	workingModel.saveModel(outputfilename.split('.')[0]+".model")

	print "\nFinal Model"
	print workingModel.prn()
	
	sys.exit()

if __name__ == "__main__":
	main(sys.argv[1:])
