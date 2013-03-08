#!/usr/bin/python

import re
import math
import sys
import random
import os
from multiprocessing import  Queue, JoinableQueue, cpu_count, Process

sys.path.append("/home/lorney/lib/libsvm-3.0/python")
sys.path.append("/home/lwalker/libsvm-3.0/python")
sys.path.append("/home/lwalker/Python-2.7.1/build/lib.linux-x86_64-2.7")


sys.path.append("../code")
from svmutil import *

from itertools import permutations
from math import sqrt, ceil
import copy

"""#################################################################
kern3D-0.4.py -- Tools for kernel analysis of 3D protein structure
Author: Lorne Walker
Changes from previous version: 
1. reviewed, merged and cleaned up code from v0.2, v0.3
see 01codeNotes.txt for details
##################################################################"""

#Initialize processor information
global NUMPROC
NUMPROC = cpu_count()-1
#NUMPROC = 4 #Force to 4 cores for debugging purposes

def queueManager(numProc, myList, function, *args):
	'''queueManager(numProc, myList, function, *args):
	generic function used to start worker processes via the multiprocessing Queue object
	numProc - number of processors to use
	myList - a list of objects to be iterated over
	function - target function
	*args - additional arguments to pass to function

	Return - an unordered(?) list of the results from myList
	'''
	qIn = Queue()
	qOut = JoinableQueue()
	if args:
		arguments = (qIn, qOut,) + args
	else:
		arguments = (qIn, qOut,)
	results = []
	
	# reduce processer count if proc count > files
	
	i = 0
	for l in myList:
		qIn.put((i,l))
		i += 1

	for _ in range(numProc):
		p = Process(target = function, args = arguments).start()
	sys.stdout.write("Progress: {:>3}%".format(0)
)
	curProgress = 0
	lastProgress = 0
	while qOut.qsize() < len(myList):
		#sys.stdout.write("\b\b\b\b{:>3}%".format(int(ceil(100*qOut.qsize()/len(myList)))))
		curProgress = int(ceil(100*qOut.qsize()/len(myList)))
		if curProgress - lastProgress > 10:
			lastProgress += 10
			sys.stdout.write("\nProgress: {:>3}%".format(lastProgress))
			sys.stdout.flush()
	sys.stdout.write("\nProgress: {:>3}%".format(100))
	#sys.stdout.write("\b\b\b\b{:>3}%".format(100))
	sys.stdout.write("\n")
	for _ in range(len(myList)):
		# indicate done results processing
		results.append(qOut.get())
		qOut.task_done()
	#tell child processes to stop
	for _ in range(numProc):
		qIn.put('STOP')

	orderedRes = [None]*len(results)
	for i, res in results:
		orderedRes[i] = res

	qOut.join()

	qIn.close()
	qOut.close()
	return orderedRes

def saveKmats(kmls, fn):
	"""saveKmats(kmls, fn):
	Takes an array of kernel matricies and saves them to file
	kmls -- an array of kernel matrix objects
	fn -- the filename to save to
	
	No return value
	"""
	outFile = open(fn, 'w')
	for km in kmls:
		outFile.write(km.pyDump())
	outFile.close()

def loadKmats(fn):	
	"""loadKmats(fn)
	Takes a file of the format created by saveKmats and loads it as an array of 
	kernel matrices
	
	fn - filename
	
	Return:array of kernel matrices	
	"""
	kRes = []
	inFile = open(fn,'r')
	for line in inFile:
		kRes.append(kMatrix.pyLoad(line))
	inFile.close()
	return kRes


def pValCalc(allCross, onVal):
	"""pValCalc(allCross, onVal)
	Calculates permutation p-value
	
	Manditory Parameters:
	allCross -- The values from the null distribution
	onVal -- the achived value for comparison
	
	Return:
	The estimated p-value
	"""
	
	allCross.sort()
	myInd = 0
	for cr in allCross:
		if float(cr) < float(onVal):
			myInd += 1
	rval = float(len(allCross) + 1 - myInd)/len(allCross)
	return min([rval, 1-rval])
	
def uniqueify(arr):
	""" uniqueify(arr)
	Return an array containing only the unique elements of original array (not necessiarially in the same order as passes
	
	Please know that this will eliminate equivalent primitives or duplicates of the same actual object.  It will not "uniqueify" two distinct objects that evaluate to equality.
	"""
	myDict = {}
	for a in arr:
		myDict[a] = 1
	return(myDict.keys())

def meanSDNoNone(arr):
	"""meanSDNoNone(arr):
	Evaluate the mean and standard deviation of a list which may contain "None" values
	
	Manditory parameters:
	arr -- the array of values to be evaluated
	
	Returns:	
	A tuple of the form (mean, SD)
	Note: raises a ValueError if all passed values are "None"
	"""
	dnm = 0
	mySum = 0
	for a in arr:
		if not(a is None):
			mySum += a
			dnm += 1
	if dnm == 0:
		raise ValueError
	myMean = float(mySum)/dnm
	
	mySqDev = 0
	for a in arr:
		if not(a is None):
			mySqDev += pow(a - myMean, 2)
	mySD = sqrt(mySqDev/(dnm-1))

	return (myMean, mySD)
	
def close_enough(x1,x2): 
	"""close_enough(x1,x2)
	an equality-like comparitor for floats

	Manditory parameters:
	x1,x2 -- two floating point values
	
	Return:
	A boolean value indicating whether the values are nearly equal
	"""
	return abs(x1-x2) < 0.0001
	
def prepareKM(myModel):
	"""prepareKM(myModel)
	Several possible uses: 1) combine kernel matricies by (possibly weighted) summation 2)
	format KM(s) for passing to libSVM
	
	Manditory Arguments:
	myModel -- a model object containing the needed parameters
	
	Returns:
	A list of dictionarys suitable for passing to libsvm as a custom kernel matrix
	
	"""
	if not(myModel.calcKM is None):
		return myModel.calcKM
	
	kmLs = []
	for i in myModel.mySlice:
		kmLs.append(myModel.km[i])
		
	weights = myModel.weights
	
	#Error checking, make sure all the parameters are consistent
	l = len(kmLs)
	if not(weights is None) and len(weights) != l:
		print "Error: prepareKM.  kmLS of length {0}, weights of length {1}".format(l,len(weights))
		raise ValueError
	if not(weights is None) and min(weights) < 0:
		print "Error: prepareKM. negative weights not allowed"
		raise ValueError
	
	#Initialize return value
	ret = []
	d = kmLs[0].dim
	
	#Handle weights
	if weights is None:
		weights = [1] * l
	
	#sum up kernel matrix
	for i in range(0,d): #For each row in the matrix
		myLine = {}
		myLine[0] = float(i+1) #The first value is a line counter 1-indexed
		for j in range(0,d): #For each column
			myTot = 0 
			for k in range(0,l): #Sum over all the matricies
				myTot += kmLs[k].getVal(i,j) * weights[k]
			myLine[j+1] = myTot
		ret.append(myLine)
	myModel.calcKM = ret
	return ret

def regularizationTune(mod):
	"""regularizationTune(mod)
	Iteratively tries a number of regularization parameters via independent cross-validation runs
	
	Paremeters:
	mod: the model you wish to tune
	
	Return:
	A model object with fine-tuned regularization parameter
	"""
	rng = [0.1,0.33,0.67,1.0,3.3,6.7,10.0]

	newMod = mod
	oldMod = model(mod.km,[], mod.pr, mod.svmType, mod.y)
	i = 1
	
	print "Starting with {}".format(newMod.accuracy())
	while newMod.accuracy() > oldMod.accuracy():
		print "Regularization round {}".format(i)
		oldMod = newMod
		regOpt = [i * oldMod.pr for i in rng]
		modelList = []
		for aReg in regOpt:
			modelList.append(model(mod.km, copy.copy(mod.mySlice), aReg, mod.svmType,mod.y))
			
		newModel = candidateListEvaluate(modelList)[0]
		i += 1

	return oldMod

def simpleBidirectionalSelection(kernMats, myParam, svmType, y, maxIter = None, start="none"):
	""" simpleBiirectionalSelection(kernMats, myParam, svmType, y, maxIter=25, start="none")
	Performs a simple feature selection search, starting with either all or none of the 
	avaliable feature.  In each step consider all possible models arrived at by
	adding or removing a single feature, pick the one with best cross-validation accuracy
	stop when current model is the best amongst possibilities, or maximum number of 
	iterations reached.

	Parameters:
	kernMats: an array of kMatrix objects, representing the avaliable features
	myParam: Regularization parameter
	svmType: "nu" or "C" designating SVM type
	y: response vector {1, -1}
	maxIter: maximum number of iterations to perform (None = no maximum)
	start: either "none" or "all" designating whether to start with all or no features

	Return:
	A model object representing the product of the selection	
	"""

	start=start.lower()

	if start == "none":
		curMod = model(kernMats,[], myParam, svmType, y)
	elif start == "all":
		curMod = model(kernMats,range(len(kernMats)), myParam, svmType, y)
	else:
		print "start value of {} passed to simpleBidirectionalSelection, should be \"none\" or \"all\"".format(start)
		raise ValueError

	runSVM(curMod)

	continueFlag = True
	i = 1

	while continueFlag and (maxIter is None or i <= maxIter):
		print "Model selection round {}".format(i)
		nextModUp = simpleDirectionalSelection(curMod, selectType = "up")
		nextModDown = simpleDirectionalSelection(curMod, selectType = "down")
		
		if nextModDown.accuracy() >= curMod.accuracy():
			if nextModUp.accuracy() > nextModDown.accuracy():
				curMod = nextModUp
			else:
				curMod = nextModDown
		elif nextModUp.accuracy() > curMod.accuracy():
			curMod = nextModUp
		else:
			continueFlag = False				
		i += 1
	return curMod

def simpleDirectionalSelection(baseModel, selectType = "up", nfold=10):
	"""simpleDirectionalSelection(baseModel, selectType="up", nfold=10,)
	Feature selection function.  Pass a current working model. Depending on the select type,
	either perform one step forward or reverse selection, returning a single best model with
	one more or one fewer locus included.
	
	Manditory Parameters:
	baseModel: our current working model
	
	Optional Parameters:
	nfold: the fold-validation to use while selecting features
	selectType: "up" for adding a single site, "down" for subtracting a single site.
	
	retuns a tuple containing (array containing indexes included in best model, model (cross-validation) accuracy)
	"""
	selectType = selectType.lower()
	n = len(baseModel.km)
	candList = []
	
	if selectType == "up":
		for i in range(0,n): #make a list of all the augmented models we wish to evaluate
			if baseModel.mySlice.count(i) == 0:
				myModel = model(baseModel.km, copy.copy(baseModel.mySlice), baseModel.pr, baseModel.svmType, baseModel.y)
				myModel.addKM(i)
				candList.append(myModel)
				
	elif selectType == "down": #make a lost of all the simplified models we wish to evaluate
		for i in baseModel.mySlice:
			myModel = model(baseModel.km, copy.copy(baseModel.mySlice), baseModel.pr, baseModel.svmType, baseModel.y)
			myModel.rmKM(i)
			candList.append(myModel)
			
	else:
		print "Un-recognized selectType parameter: {0}".format(selectType)
		raise ValueError
	
	try:
		return candidateListEvaluate(candList)[0]
	except IndexError:
		return model(None, [], None, None, None, None)


def candidateListEvaluate(candList, numReturn=1):
	"""candidateListEvaluate(candList, numReturn=1)

	Manditory Parameters:
	candList: a list of candidate model objects to evaluate
	
	Optional Parameters:
	numReturn: Then number of best models to return
	
	Return Value:
	An array of best fitting model objects, up to length of numReturn
	"""
	global NUMPROC #Number nodes to split work out into

	if len(candList) <= numReturn: #If we ask for more results than we have options
		return candList

	allResults = queueManager(NUMPROC, candList, runSVMqueue)

	for i in range(len(allResults)):
		candList[i].acc = allResults[i]

	return [x for (y,x) in sorted(zip(allResults,candList), reverse = True)][:numReturn]

	
def runSVMqueue(qIn, qOut):
	"""runSVMqueue(qIn, qOut)
	Wrapper function for multiprocessing SVMs

	"""
	for i, m in iter(qIn.get, 'STOP'):
		qOut.put((i, runSVM(m)))


def runSVM(myModel):
	"""runSVM(myModel, qIn= None, qOut=None):
	Evaluate the 10-fold crossvalidation accuracy of a ksvm based on the subset of 
	kernel matrices indicated by the myModel variable (or contained in the qIn variable)

	Manditory Parameters:
	myModel - a model object that describes the candidate we're evaluating

	Returns:
	SVM cross-validation accuracy (float)
	
	"""
		
	#The rest is for evaluation of a single model
	if len(myModel.mySlice) == 0:
		return 0

	x = prepareKM(myModel)
	
	prob = svm_problem(myModel.y, x)
	paramFlag = "-n"
	svmFlag = 1
	if myModel.svmType == "C":
		paramFlag= "-c"
		svmFlag = 0
	elif myModel.svmType != "nu":
		print "Defaulting to nu-SVM"
	crossFlag = "-v {}".format(myModel.nfold)
	if myModel.nfold == 1:
		crossFlag = ""
	param = svm_parameter("-q -s {} -t 4 {} {} {} ".format(svmFlag,crossFlag,paramFlag, myModel.pr))
	try:
		if myModel.nfold == 1:
			mod = svm_train(prob, param)	
			cross = svm_predict(myModel.y, x, mod)[1][0]
		else:
			cross = svm_train(prob,param)
	except ValueError:
		cross = None

	return cross

def getNextSeq(fileLines):
	"""getNextSeq(fileLines)
	pull (and parse) the next sequence out of the list of lines from a fasta file
	
	Manditory parameters:
	fileLines -- the lines from a fasta file
	
	Return:
	a tuple: (sequence name, AA sequence)
	"""
	while(fileLines[0].find(">") < 0 and len(fileLines) > 0):
		fileLines.pop(0)
	if(len(fileLines) == 0):
		seqName = None
	else:
		seqName = fileLines.pop(0).strip().strip(">")
	inputSeq = ""
	while(len(fileLines) > 0 and fileLines[0].find(">") < 0):
		inputSeq += fileLines.pop(0).strip()
	return (seqName,inputSeq)
			
def parseY(respFileName):
	"""parseY(respFileName)
	Prepare the response data
	
	Manditory parameters:
	respFileName -- file name of the file containing response data
	
	Returns: 
	dictionary with sample names as keys and 1 or -1 as values
	"""
	resFile = open(respFileName, 'r')
	resFileLines = resFile.readlines()
	resFileLines.pop(0) #Assume the first line is a label

	retVal = {}
	for line in resFileLines:
 		line = line.strip("\n").strip().split(',')
		if line != "": 
			retVal[line[0].strip()] = line[1].strip()

	return(retVal)
			
def readInConsensus(fLines):
	"""readInConsensus(fLines)
	Read in and parse the consensus sequence from the fasta file lines
	
	Mandatory Parameters:
	fLines -- the list of lines from the fasta file
	
	Return:
	a refSeq object
	"""
	refName, referenceAASeq = getNextSeq(fLines)
	print "Identified reference sequence: "+refName
	return refSeq(referenceAASeq)
	
def parseAAmap(aaFileName):
	"""parseAAmap(aaFileName)
	Read in an amino acid map file
	
	Manditory parameters:
	aaFileName -- filename for AA map
	
	Return:
	a dictionary of dictionaries of similarity scores between AA's
	"""
	aaFile = open(aaFileName, 'r')
	aaLines = aaFile.readlines()
	aaFile.close()

	aaMap = {}

	headerLine = aaLines.pop(0)
	headerLine = headerLine.strip().split(',')

	i = 0
	for line in aaLines:
		myName = headerLine[i]
		line = line.strip().split(',')
		aaMap[myName] = {}
		j =0
		for item in line:
			try:
				myVal = float(item)
			except ValueError:
				"Non-numeric entry in amino acid match file: "+str(item)
				sys.exit(2)
			aaMap[myName][headerLine[j]] = myVal
			j += 1
		i += 1
		
	return aaMap

def reCalcPOS(myRad, sigma, pdbFN, fastFN, mapFN=None, AAmap=None, nrm=False, gp41filename=None):
	"""reCalcPOS(myRad, sigma, pdbFN, fastFN, mapFN, AAmap, nrm=False, gp41filename=None)
	Calculate the positional kernel matrixes based on the input parameters

	Manditory parameters:
	myRad -- the radius of the positional bubble in Angstroms
	sigma -- the standard deviation of the penalty function for AA mismatch
	pdbFN -- the filename of the gp120 PDB filename
	fastFN -- the filename of the fastafile.  First row is consensus (only row read if a map included)
	
	Optional parameters:
	mapFN -- a filename for map of sample sequences onto consensus
	AAmap -- a dictionary of dictionaries, expressing the pairwise AA similarity on a [0,1] scale (default: None)
	nrm -- a boolean indicating whether to normlalize the kernel matrix or not (default: False)
	gp41filename -- a filename for a gp41 pdb file, not manditory (default: None)
	
	Return:
	a kMatResult object
	"""
	global NUMPROC
	
	#Input FASTA file
	try:
		fastaFile = open(fastFN,'r')
	except IOError:
		print "\nCouldn't open fasta file.\n"
		sys.exit(2)
	fastaLines = fastaFile.readlines()
	fastaFile.close()

	#Read in consensus sequence
	referenceSeq = readInConsensus(fastaLines)

	#Map consensus sequence onto pdb coordinates
	if not (pdbFN is None):
		referenceSeq.associateCoordinates(pdbFN,1)
		
	#Add gp41 coordinates, if avaliable
	if not (gp41filename is None):
		referenceSeq.associateCoordinates(gp41filename,2)

	#Set up sample sequences
	sampleSeqs = []
	
	if not(mapFN is None): #If there's a pre-calculated map
		print "Importing map from "+mapFN+". . . "
		d = {}
		execfile(mapFN,d,d)
		sequenceMaps = d["sequenceMaps"]
		for k in sequenceMaps.keys():
			sampleSeqs.append(seq(k,sequenceMaps[k]))
	else: #Estimate directly from the remaining sequences in the fasta file
		while (len(fastaLines) > 0):
			myName, mySeq = getNextSeq(fastaLines)
			if((not myName is None) and (mySeq != "")):
				sampleSeqs.append(seq(referenceSeq, myName, mySeq))
	print "Imported "+str(len(sampleSeqs))+" sample sequences"

	print "Calculate kernel matrixes:"
	#Calculate the KM for each position in the reference sequence
	kernResult = queueManager(NUMPROC, range(referenceSeq.seqLen), kernMatCalcQueue, referenceSeq, sampleSeqs, myRad, sigma, AAmap)
	
	if nrm: #Normalize if requested
		for km in kernResult:
			km.normalize()
	return(kernResult)
	
def kernMatCalcQueue(qIn, qOut, rSeq, sSeqs, rad, sig, AAmap):
	"""kernMatCalcQueue(qIn, qOut, rSeq, sSeqs, rad, sig, AAmap)
	Wrapper method for using queueManager to run kernel matrix calculations
	"""
	for i, p in iter(qIn.get, 'STOP'):
		qOut.put((i,kernMatCalc(p, rSeq, sSeqs, rad, sig, AAmap)))
	
def kernMatCalc(pos,rSeq,sSeqs,rad,sig,AAmap):
	"""kernMatCalc(iterObj):
	A parallelizable method for calculating a single kernel matrix
	
	Manditory Parameters:
	iterObj -- a tuple or list that contains the parameters in the following order:
	(pos,rSeq,sSeqs,rad,sig,AAmap)
	pos is the zero-indexed position in the reference sequence at which we're calculating the kernel matrix
	
	See description for calculatePositional KM for description of other parameters
	
	Return:
	A tuple containing (kMatrix object, averageWordLength)
	"""
	myKmat = kMatrix(len(sSeqs))
	wordDict = {}
	wordList = []
	
	seqNames = []

	for i in range(0,len(sSeqs)): # A dictionary of what words are in which samples
		seqNames.append(sSeqs[i].name)
		myWord = sSeqs[i].word(pos, rSeq, rad)
		try:
			matchWord = wordList[wordList.index(myWord)]
			wordDict[matchWord].append(i)
		except ValueError:
			wordDict[myWord] = [i]
			wordList.append(myWord)

	for word1 in wordDict.keys(): 
		for word2 in wordDict.keys():
			my_val = wordKernelFunction(word1,word2,sig,AAmap) #This calculates the kernel function
			for ind1 in wordDict[word1]:
				for ind2 in wordDict[word2]:
					myKmat.putVal(ind1,ind2,my_val)

	myKmat.sampleNames = seqNames
	myKmat.rad = rad
	myKmat.sd = sig
	myKmat.pos = pos

	return myKmat	
	
def dist(x1,x2,ind1,ind2,chain1,chain2):
	"""dist(x1,x2,ind1,ind2,chain1,chain2)
	Calculates the distance between two points
	
	Mandatory Parameters:
	x1,x2 -- the coordiantes of the two points
	ind1, ind2 -- the indeces of the two points
	chain1, chain2 -- indicators of which chain the coordinates come from 
	
	Return:
	the distance in Angstroms
	"""
	if x1 == (-1,-1,-1) or x2 == (-1,-1,-1) or chain1 != chain2 or chain1 == -1:
		return abs(ind1-ind2) * 3.8
	else:
		return math.sqrt(pow(x1[0] - x2[0],2) + pow(x1[1] - x2[1],2) + pow(x1[2] - x2[2],2))

def wordKernelFunction(w1,w2,sig,AAmap=None):
	"""wordKernelFunction(w1,w2,sig,AAmap=None)
	Given two word objects, calculate the kernel function

	Manditory parameters: 
	w1,w2 -- the two words to compare by kernel function
	sig -- the standard deviation of the distance penalty function

	Optional parameters:
	AAmap -- a dictionary of dictionaries holding the seq objects (default: None)

	Return:
	a float value for the kernel function
	"""
	total_val = 0
	if (len(w1) == 0) and (len(w2) == 0):
		return 1
	elif len(w1) == 0:
		return 0
	elif len(w2) == 0:
		return 0
	else:
		for i in range(0,len(w1)):
			for j in range(0,len(w2)):
				if w1.seq[i] == w2.seq[j]:
					myMult = 1
				else: 
					myMult = 0
				if not AAmap is None:
					try:
						myMult = AAmap[w1.seq[i]][w2.seq[j]]
					except (KeyError, TypeError):
						pass
				myDist = dist(w1.coord[i],w2.coord[j],w1.pseudoIndex[i],w2.pseudoIndex[j], w1.chain[i], w2.chain[j])
				total_val += myMult*math.exp(-math.pow(myDist,2)/(2*sig))
		return float(total_val)/(len(w1)*len(w2))

class refSeq(object):
	"""class refSeq
	A class for a reference sequence
	
	Constructor:
	refSeq() for an empty object
	refSeq(seq): seq = character string for the reference sequence
	
	Methods:
	associateCoordinates(self, pdL,i)
	pseudoCoord(self, pseudoInd)
	pseudoChain(self,pseudoInd):
	"""

	def __init__(self, *params):
		"""	Constructor for class refSeq:
	
		Optional Parameters:
		seq -- a character string of the reference sequence
	
		Returns:
		refSeq object
		"""
		if len(params) < 1:
			inSeq = ""
		else:
			inSeq = params[0]
		
		#Instance variables
		#seq: a string containing the AA sequence with gaps removed
		#seqAlignIndex: an array with the indexes of all non-gaps in the reference sequence
		#seqLen: length of the AA sequence
		self.alignLen = len(inSeq)		
		self.seqAlignIndex = [match.start() for match in re.finditer(r"[A-Z]",inSeq)]
		self.seq = re.sub(r"[-\*\.]","",inSeq)
		self.seqLen = len(self.seq)
		self.coord = self.seqLen*[(-1,-1,-1)]
		self.chain = self.seqLen*[-1]
		#Check to ensure we have an index for each AA
		if self.seqLen != len(self.seqAlignIndex):
			raise ValueError("Invalid reference sequence.  Should only contain A-Z and .*-")
			
	def associateCoordinates(self, pdbFilename,chainNum):
		"""refSeqInstance.associateCoordinates(pdbFilename,chainNum)
		Map the reference sequence onto the pdb lines
	
		Mandatory parameters:
		pdbFilename: filename of the pdb file
		chain: chain number to assign coordinates to
		
		Return: no value returned	
		"""
		try:
			pdbFile = open(pdbFilename,'r')
		except IOError:
			print "\nCouldn't open pdb file.\n"
			sys.exit(2)
		pdbL = pdbFile.readlines()
		pdbFile.close()
		
		findAtom = re.compile("ATOM")
		for line in pdbL:
			line = line.strip()
			if not(findAtom.match(line) is None):
				atomType = line[12:15].strip()
				if (atomType == "CA") and (not(re.search("[ABC]",line[21]) is None)):
					atomIndex = int(line[23:26].strip())
					atomX = float(line[31:38].strip())
					atomY = float(line[39:46].strip())
					atomZ = float(line[47:54].strip())		
					self.coord[atomIndex] = (atomX,atomY,atomZ)
					self.chain[atomIndex] = chainNum
	
	def pseudoCoord(self, pseudoInd):
		"""refSeqInstance.pseudoCoord(pseudoInd)
		Calculate the approximated coordinates of a residue based on it's estimated position
		
		Manditory parameters:
		pseudoInd: the pseudo-index of the position to be referenced (need not be an integer)
		
		Return: 
		a tuple (X,Y,Z)
		"""
		#print(str(pseudoInd))
		if close_enough(pseudoInd, int(pseudoInd)): #if it's a whole number index
			return self.coord[int(pseudoInd)]
		elif int(pseudoInd) + 1 > self.seqLen - 1:
			return (-1,-1,-1)
		elif self.coord[int(pseudoInd)] == (-1,-1,-1) or self.coord[int(pseudoInd)+1] == (-1,-1,-1):
			return (-1,-1,-1)
		elif self.chain[int(pseudoInd)] != self.chain[int(pseudoInd)+1]:
			return (-1,-1,-1)
		else:
			fract = pseudoInd - int(pseudoInd)
			coord1 = self.coord[int(pseudoInd)]
			coord2 = self.coord[int(pseudoInd)+1]
			return (coord1[0] + fract*(coord2[0] - coord1[0]), coord1[1] + fract*(coord2[1] - coord1[1]), coord1[2] + fract*(coord2[2] - coord1[2]))
			
	def pseudoChain(self,pseudoInd):
		"""refSeqInstance.pseudoChain(pseudoInd)
		Figure the source chain for a position in the reference sequence
		
		Manditory parameters:
		pseudoInd: the pseudo-index of the position to be referenced (need not be an integer)
			
		Return:
		Chain number, if present, -1 if no chain at that position	
		"""
		if close_enough(pseudoInd, int(pseudoInd)): #if it's a whole number index
			return self.chain[int(pseudoInd)]
		elif int(pseudoInd) + 1 > self.seqLen - 1:
			return -1
		elif self.chain[int(pseudoInd)] == self.chain[int(pseudoInd)+1]:
			return self.chain[int(pseudoInd)]
		else:
			return -1

class seq(object):
	"""class seq:
	A class for sample sequences
	
	Constructor:
	seq(name, seqMap) -- for sequences with pre-calculated position map onto reference
	seq(refSeq, name,  mySeq) -- for for calculation directly from the sequence
	name -- sequence name
	AND
	seqMap -- a dictionary of position and pseudo-index
	OR
	refSeq -- the reference sequence to map to
	mySeq -- amino acid sequence in question
	
	Instance Methods:
	calcPseudoCoords(self, rSeq)
	mapPseudoIndex(self, rSeq)
	word(self, pos, rSeq, rad)
	"""
	def __init__(self,*params): #If 3 parameters, called directly from alignment, if 2, from HyPhy Map.  Kind of hackish way to do it.  Should consider switching to named parameters.
		"""Constructor for seq object
		seq(name, seqMap)
		OR
		seq(refSeq, name, mySeq)
		
		Mandatory parameters:
		name -- Sequence name
		AND
		seqMap -- a dictionary of position and pseudo-index
		OR
		refSeq -- the reference sequence to map to
		mySeq -- amino acid sequence in question
		
		Return:
		seq object
		"""
		if len(params) == 2: #map from HyPhy
			self.name, myDict = params
			self.seqAlignIndex = None #not from an alignment
			self.seq = ""
			self.pseudoIndex = []
			moreFlag = True
			i = 0
			while moreFlag:
				try:
					nextResidue = myDict[str(i)]
					self.seq += nextResidue['AA']
					self.pseudoIndex.append(nextResidue['C'])
					i += 1
				except KeyError:
					moreFlag = False
			self.seqLen = len(self.seq) 
			
		elif len(params) == 3: #Directly from alignment
			rSeq, self.name, mySeq = params
			self.alignLen = len(mySeq)
			self.seqAlignIndex = [match.start() for match in re.finditer(r"[A-Z]",mySeq)]
			self.seq = re.sub(r"[-\*\.]","",mySeq)
			self.seqLen = len(self.seq)
			self.pseudoIndex = self.seqLen*[-1]
	
			self.mapPseudoIndex(rSeq)
		else:
			raise ValueError("seq constructor requires 2 or 3 parameters, gave " + int(len(params)))
			
		self.pseudoCoords = []
		self.pseudoChain = []
		
	def calcPseudoCoords(self, rSeq):
		"""seqInstance.calcPseudoCoords
		Updates the seq object's pseudo-coordinates by mapping it to a reference sequence
	
		Mandatory parameters:
		rSeq -- reference sequence object
	
		Return: no returned value
		"""

		for i in range(0,self.seqLen):
			self.pseudoCoords.append(rSeq.pseudoCoord(self.pseudoIndex[i]))
			self.pseudoChain.append(rSeq.pseudoChain(self.pseudoIndex[i]))

	def mapPseudoIndex(self, rSeq):
		"""seqInstance.mapPseudoIndex
		Mapping the sequence pseduo indexes (fractional indexes to account for indels) according to reference sequence
		
		Mandatory parameters:
		rSeq -- reference sequence object
		
		Return: no returned value
		"""
		i = 0
		j = 0
		#iterate over each position in the sequences
		while(i < self.seqLen and j < rSeq.seqLen):
			#if both indicies point to the same alignment position
			if self.seqAlignIndex[i] == rSeq.seqAlignIndex[j]: 
				self.pseudoIndex[i] = j
				i += 1
				j += 1
			#introducing an gap into the sample sequence
			elif self.seqAlignIndex[i] > rSeq.seqAlignIndex[j]: 
				start_j = j
				while(self.seqAlignIndex[i] > rSeq.seqAlignIndex[j] and j < rSeq.seqLen):
					j += 1
				#This means that the locations of residues (i-1):i should be be spread out evenly amongst the locations (j_start-1):j
				self.pseudoIndex[i-1] = start_j - 2 + (j-start_j+3.0)/3
				self.pseudoIndex[i] = start_j- 2 + 2*(j-start_j+3.0)/3
				i += 1
				j += 1
			#introducing an insertion into the sample sequence
			elif self.seqAlignIndex[i] < rSeq.seqAlignIndex[j]:
				start_i = i
				while(self.seqAlignIndex[i] < rSeq.seqAlignIndex[j] and i < self.seqLen):
					i += 1
				#This means that residues start_i:(i-1) should be fitted in between (j-1) and j
				insertLength = i-start_i	
				for i_ins in range(0,insertLength):
					self.pseudoIndex[start_i+i_ins] = j + (i_ins + 1.0)/(insertLength + 1.0)
					
	def word(self, pos, rSeq, rad):
		"""seqInstance.word(pos, rSeq, rad):
		Determine the Amino-Acid word centered at a particular point in the reference sequence
		
		Manditory parameteres:
		pos -- pseudo-index to center the sphere at
		rSeq -- relevant reference sequence
		rad -- radius of sphere in Angstroms
		
		Return:
		a word object
		"""
		if len(self.pseudoCoords) == 0:
			self.calcPseudoCoords(rSeq)
		myWord = word()
		for i in range(0,self.seqLen):
			myPseudoCoord = rSeq.pseudoCoord(self.pseudoIndex[i])
			try:
				myDist = dist(rSeq.coord[pos],self.pseudoCoords[i],pos,self.pseudoIndex[i],rSeq.chain[pos],self.pseudoChain[i])
			except IndexError:
				print "Index error: pos:{}, i:{}".format(pos,i)
				sys.exit(2)

			if myDist <= rad:	
				myWord.seq.append(self.seq[i])
				myWord.coord.append(myPseudoCoord)
				myWord.pseudoIndex.append(self.pseudoIndex[i])
				myWord.chain.append(self.pseudoChain[i])
		return myWord
	
class kMatrix(object):
	"""class kMatrix
	A class for kernel matrices, stores matrix compactly due to constraint of symmetry
	
	Constructor:
	kMatrix(dimension) 
	dim: integer: number of rows columns in the matrix
	
	Instance Methods:	
	getVal(self,i,j)
	putVal(self,i,j,val)
	normalize(self)
	rString(self)
	pyDump(self)
	
	Class Methods:
	pyLoad(cls, kmStr)
	"""
	def __init__(self, *params):
		"""kMatrix(dimension)
		Constructor for class kMatrix
		
		Manditory Parameters:
		dimension - number of rows/columns to include
		"""
		if len(params) < 1:
			raise ValueError("Need to supply kMatrix constructor with dimension.")
		else:
			myDim = int(params[0])
			try:
				self.dim = myDim
				self.values = myDim*(myDim+1)/2*[0]
			except:
				raise ValueError("Invalid dimension passed to kMatrix")
		self.sampleNames = ["NA"]*self.dim
		self.rad = None
		self.sd = None
		self.pos = None
			
	def descStr(self):
		"""kMatrixInstance.descStr():
		Formats a string describing the kernel matrix

		Parameters: None

		Return: String describing object
		"""
		return "pos:{};rad:{};sd:{}".format(self.pos,self.rad,self.sd)
			
	def getVal(self,i,j):
		"""kMatrixInstance.getVal(i,j)
		Method for reading from the matrix
		
		Manditory parameters:
		i,j: coordinates of the value to retrieve
		
		Return:
		float value at position [i,j] in the matrix
		"""
		return self.get(self.values,i,j)
		
	def putVal(self,i,j,val):
		"""kMatrixInstance.putVal(i,j,val)
		Method for writing to the matrix
		
		Manditory parameters:
		i,j: coordinates of the value to set
		val: float value to place at that position		
		
		Returns: no return value
		"""
		self.put(self.values,i,j,float(val))
		
	def get(self,vals, i,j):
		"""kMatrixInstance.get(values,i,j)
		Internal method for getting kernel matrix values
		"""
		if j > i:
			return self.get(vals,j,i)
		else:
			return vals[i*(i+1)/2+j]
	
	def put(self,vals,i,j,val):
		"""kMatrixInstance.put(values,i,j,val)
		Internal method for setting kernel matrix values
		"""	
		if j > i:
			self.put(vals,j,i,val)
		else:
			vals[i*(i+1)/2+j] = val
			
	def normalize(self):
		"""kMatrixInstance.normalize()
		Normalizes the kernel matrix in place
		
		No parameters or return value		
		"""
		newValues = self.dim*(self.dim+1)/2*[0]
		for i in range(0,self.dim):
			for j in range(0,i+1):
				myVal = self.getVal(i,j)/sqrt(self.getVal(i,i)*self.getVal(j,j))
				self.put(newValues,i,j,myVal)
				
		self.values = newValues

	def rString(self):
		"""kMatrixInstance.rString()
		Convert kernel matrix into a string for importing into R
		
		No parameters
				
		Return: an R-formatted string
		"""
		myArr = []
		for i in range(0,self.dim):
			for j in range(0,self.dim):
				myArr.append(str(self.getVal(i,j)))
		return "structure(c("+", ".join(myArr)+"), .Dim=c("+str(self.dim)+"L, "+str(self.dim)+"L))"
		
	def pyDump(self):
		"""kMatrixInstance.pyDump()
		Return a string representation of kernel matrix that can be written to file, loadable through pyLoad
		
		No parameters
		
		Return: a python-formatted string
		"""
		myStr = ":".join([str(self.dim), str(self.pos), str(self.rad), str(self.sd)])+":"
		myStr = myStr + ','.join(self.sampleNames) + ":"
		myStr = myStr + ','.join(map(lambda s: str(s), self.values)) + "\n"
		return myStr
		
	def pyLoad(cls, kmStr):
		"""kMatrix.pyLoad(kmStr)
		
		Manditory Parameters:
		kmStr -- a string of the format created by pyDump
		
		Return: 
		a kMatrix object		
		"""
		try:
			dim, pos, rad, sd, names, val = kmStr.strip("\n").split(":")
			myKM = kMatrix(int(dim))
			myKM.rad = float(rad)
			myKM.sd = float(sd)
			myKM.pos = int(pos)
			myKM.sampleNames = names.split(",")
			myKM.values = map(lambda s: float(s), val.split(','))
			return myKM
		except ValueError:
			print "Invalid kernel matrix file format"
			sys.exit(2)
			
	pyLoad = classmethod(pyLoad)

class word(object):
	"""class word
	A class for amino-acid words

	Constructtor: word()
	
	Instance Methods:
	len(wordInstance)
	prrnt(self)	
	"""
	def __init__(self, *params):
		"""word()
		Initializer for the class word
		
		No parameters
		"""
		self.seq = []
		self.pseudoIndex = []
		self.coord = []
		self.chain = []
	def __len__(self):
		"""len(wordInstance)
		Overrides the len() operator
		
		Return:
		Length in amino acids of described word
		"""
		return len(self.seq)
	def prnt(self):
		"""wordInstance.prnt()
		Provides a string of information about this particular word
		
		Return:
		string of the form "CDE: (1,2,3)"
		"""
		outString = "".join(self.seq)+": ("
		lst = []
		for pI in self.pseudoIndex:
			lst.append(str(pI))
		outString += ",".join(lst)
		outString += ")"
		return outString
		

class model(object):
	""" class model
	An object representing a given model, designates the kernel matrix (or composite of matricies), outcome variable and model parameters.  Once the model has been tested, can also remember its cross-validation result.
	"""
	def __init__(self, km, mySlice, pr, svmType, y, nfold = 10, weights = None):
		"""model constructor
		Mandatory parameters
		km: list of kernel matrices to be considered
		mySlice: list of indicies of matrices to be included in the model
		pr: model regularization parameter
		svmType: either "nu" or "C" type svm
		y: response vector
		
		Optional parameters:
		nfold: number of fold validation to use to assess the model
		weights: if we want to make our predictor kernel matrix a weighted sum of the contributors, pass an array of weights here
		"""
		self.km = km
		self.mySlice = mySlice
		self.weights = weights
		self.sortMe()
		self.pr = pr
		self.svmType = svmType
		self.y = y
		self.nfold = nfold
		self.acc = None
		self.calcKM = None
		
		#Quick checks for parameter consistency
		if not(weights is None) and len(weights) != len(mySlice):
			print "Inconsistent length of weights {} and km index {} vectors".format(len(weights),len(mySlice))
			raise ValueError
		if len(mySlice) > 0 and max(mySlice) >= len(km):
			print "Out-of-bounds kernel matrix index: {}, length of options: {}".format(max(mySlice),len(km))
	
	def saveModel(self, fileName):
		""" model.save(fileName)
		Saves model in a file format readable by the loadModel method

		Parameters: fileName -> name of a file to save to

		Return: None
		"""
		outf = open(fileName, 'w')
		outf.write("y;{}\n".format(','.join(map(lambda i: str(i), self.y))))
		if not(self.weights is None):
			outf.write("weights;{}\n".format(','.join(map(lambda i: self.weights[i], self.mySlice))))
		if not(self.acc is None):
			outf.write("acc;{}\n".format(self.acc))
		outf.write("nfold;{}\n".format(self.nfold))
		outf.write("svmType;{}\n".format(self.svmType))
		outf.write("pr;{}\n".format(self.pr))
		for i in self.mySlice:
			outf.write(self.km[i].pyDump())
		outf.close()
		
		
	def loadModel(cls, fileName):
		"""model.loadModel(fileName)
		Loads model from file format as created by saveFile method

		Parameters: fileName -> name of file to load

		Return: A model object
		"""
		inf = open(fileName, 'r')
		lns = inf.readlines()
		
		kms = []
		pr = None
		svmType = None
		y = None
		nfold = None
		weights = None
		acc = None
		
		for l in lns:
			l = l.strip().split(";")
			if l == "":
				pass
			elif len(l) == 2:
				if l[0] == "y":
					y = map(lambda i: int(i),l[1].split(','))
				elif l[0] == "acc":
					acc = float(l[1])
				elif l[0] == "weights":
					weights = map(lambda i: float(i),l[1].split(','))
				elif l[0] == "nfold":
					nfold = int(l[1])
				elif l[0] == "svmType":
					svmType = l[1]
				elif l[0] == "pr":
					pr = float(l[1])
			else:
				kms.append(kMatrix.pyLoad(l[0]))
		
		myM = model(kms, range(len(kms)), pr, svmType, y, nfold, weights)
		myM.acc = acc
		
		inf.close()
		
		myM.sortMe()

		return myM 
	
	loadModel = classmethod(loadModel)
	
	def __eq__(self,other):
		"""
		model1 == model2
	
		Returns True if the two models reference the same features and have the same elements

		"""
		try:
			return self.asStr() == other.asStr()
		except AttributeError:
			return False
	
	
	def addKM(self, i, wt=None):
		"""
		modelInstance.addKM(i)
		Adds the kernel matrix at index i to the model, 
	
		Parameters:
		i -> index of the parameter to remove
		wt -> an Optional weight if we're looking at a weighted model	

		Return: no value		
		"""
		self.mySlice.append(i)
		if not self.weights is None:
			self.weights.append(wt)
		self.sortMe()
	
	def rmKM(self, res):
		"""
		modelInstance.rmKM(i)
		Removes the feature at index i (in the original list of all possible features) from the model
		
		Parameters: res -> The index of the feature to remove

		Return: No value
		"""
		try:
			i = self.mySlice.index(res)
			del self.mySlice[i]
			if not self.weights is None:
				del self.weights[i]
		except ValueError:
			print "tried to remove KM {} from model, KM not found".format(res)
	
	def accuracy(self):
		"""
		modelInstance.accuracy()
		Returns model accuracy.  If it has not been set yet, run the SVM so we know.

		Parameters: no parameters
	
		Returns: Float value <= 1 designating cross-validation accuracy
		"""
		if len(self.mySlice) == 0:
			self.acc = 0.0

		if self.acc is None:
			self.acc = runSVM(self)
		return self.acc
	
	def asStr(self):
		""" modelInstance.asStr()
		Returns a unique string that describes (uniquely) the model features and weights, to 
		see if we've already considered this model
	
		"""
		self.sortMe()
		if self.weights is None:
			return(','.join([str(x) for x in self.mySlice]))
		else:
			return(','.join([str(x) for x in self.mySlice])+":"+','.join([str(x) for x in round(self.weights,2)]))
		
	def sortMe(self):
		"""modelInstance.sortMe()
		Sorts (in place) the indexes of the contained model features, as well as the
		 associated weights.  This allows or easy determination if two models are identical.

		Parameters: no parameters

		Return: no value

		"""
		if not self.weights is None:
			self.weights = [wt for (ind,wt) in sorted(zip(self.mySlice,self.weights))]
		self.mySlice.sort()

	def prn(self):
		"""modelInstance.prn()
		Creates a string that describes the model in question in a user-friendly format
		
		Parameters: no parameters

		Return: string description of model
		"""
		myStr =  "Model includes kernel matrices #{}, with cross-validation accuracy:{}\n".format(self.asStr(), self.accuracy())
		j = 1
		for i in self.mySlice:
			if self.weights is None:
				myStr += "KM{}:{}\n".format(j, self.km[i].descStr())
			else:
				myStr += "KM{}:{}, weight:{}\n".format(j, self.km[i], self.weights[i])
			j += 1
			
		return myStr
