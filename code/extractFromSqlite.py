#!/usr/bin/python

"""extractFromSqlite.py
Script for pulling a fasta file and response data from Lance's sqlite database.  
Syntax: ./extractFromSqlite.py sqliteFileName referenceFastaFileName abName
sqliteFileName         - filename for the sqlite database file
referenceFastaFileName - filename for a fasta file where the first sequence is the reference sequence to be used for structural data
abName                 - name of the antibody to extract
"""

import sys
import sqlite3 as lite

def main(args):
	if len(args) < 3:
		print "Three parameters required, {} supplied".format(len(args))
		sys.exit(2)
		
	sqliteFileName = args.pop(0)
	referenceFastaFileName = args.pop(0)
	abName = args.pop(0)
	
	fastaOut = open("{}.fasta".format(abName), 'w')
	respOut = open("{}.csv".format(abName), 'w')
	
	
	print referenceFastaFileName
	#Get reference sequence, put it at the top of the fasta file
	refFile = open(referenceFastaFileName, 'r')
	refLines = refFile.readlines()
	refName = refLines.pop(0).strip("\n")
	refSeq = ""
	while len(refLines) > 0 and refLines[0].find(">") < 0:
		aLine = refLines.pop(0).strip("\n").strip()
		refSeq = refSeq+aLine
	fastaOut.write(refName+"\n"+refSeq+"\n")
	
	try:
		con = lite.connect(sqliteFileName)
		
		cur = con.cursor()
		
		cur.execute("SELECT * FROM neut WHERE ANTIBODY='{}';".format(abName))
		neuts = cur.fetchall()
		
		cur.execute("SELECT * FROM SEQUENCE;".format(abName))
		seqs = cur.fetchall()
	
	except lite.Error, e:
		print "Error %s:" % e.args[0]
		sys.exit(2)

	neutMap = {}
	for n in neuts:
		neutMap[str(n[1])] = str(n[3])
		
	seqMap = {}
	for s in seqs:
		seqMap[str(s[1])] = str(s[2])
		
	respOut.write("ID\tseq\tANTIBODY\tIC50\n")
	
	i = 0
	for aID in neutMap.keys():
		try:
			mySeq = seqMap[aID]
			respOut.write("\t".join([aID,mySeq,abName,neutMap[aID]])+"\n")
			fastaOut.write(">{}\n{}\n".format(aID,mySeq))
			i += 1
		except KeyError:
			print "No sequence for {}".format(aID)

	print "Found {} sequence/titer pairs\n".format(i)

	fastaOut.close()
	respOut.close()
	
if __name__ == "__main__":
	main(sys.argv[1:])	
