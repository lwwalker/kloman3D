#!/usr/bin/python

from kern3Dv04 import *
import getopt
import sys
import random
import time

def usage():
	print "\nUsage: permutationTest.py [optional parameters] kmfilename responsefilename outfilename"
	print "-p --perm -> Number of permutations to test (default 100)"
	print "-c --cSVM -> perform c-svm with the passed regularization parameter"
	print "-n --nuSVM -> perform nu-svm with the passed regularization parameter (default with nu = 0.2)"
	print "-x --nocross -> use non-cross-validation SVMs for feature selection"

def main(argv):
	argString = "permutationTest.py "+" ".join(argv)+" \nCalled at " + time.strftime("%a, %d %b %Y %H:%M:%S") + "\n"

	#Initialize defaults
	svmType = "nu"
	myParam = 0.2
	nocross = False
	saveKMfilename = None
	nPerm = 100
	
	#Read in command line arguments
	try:                                
		opts, args = getopt.getopt(argv, "xp:c:n:", ["nocross","perm","cSVM","nuSVM","saveKM"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
		
	for opt, arg in opts:
		if opt in ("-x","--nocross"):
			nocross = True
		if opt in ("-p","--perm"):
			try:
				nPerm = int(arg)
			except ValueError:
				print "Invalid permutation number passed"
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
	
	if len(args) < 3:
		usage()
		sys.exit(2)
	
	kmfilename = args.pop(0)	
	responsefilename = args.pop(0)
	outputfilename = args.pop(0)	
	
	outfile = open(outputfilename, 'w')
	sys.stdout = outfile
	print argString
	sys.stdout.flush()	
	
	print("Initalizing . . .")
	allResp = parseY(responsefilename)	#Load in titers

	print "Loading kernel matrixes . . ."
	myKM = loadKmats(kmfilename)
	
	#Associate titer data with samples in kernel matrixes
	y = []
	for name in myKM[0].sampleNames:
		try:
			y.append(int(allResp[name]))
		except KeyError:
			print "Unable to find response value for sequence \"{}\"".format(name)
			sys.exit(2)	
			
	includedSites = {}
	finalAcc = []
	for i in range(nPerm):
		print "Running permutation {}".format(i+1)
		random.shuffle(y)
		thisModel = simpleBidirectionalSelection(myKM, myParam, svmType, y)
		tunedModel = regularizationTune(thisModel)
		finalAcc.append(tunedModel.accuracy())
		for site in tunedModel.mySlice:
			try:
				includedSites[site] += 1
			except KeyError:
				includedSites[site] = 1
				
	print "Accuracy:"
	print finalAcc
	print "Max: {}, Mean: {}".format(max(finalAcc),mean(finalAcc))
	print "Included Sites:"
	print includedSites
	
				
def mean(l):
	val = float(sum(l))/len(l) if len(l) > 0 else float('nan')
	return val
	
if __name__ == "__main__":
	main(sys.argv[1:])
