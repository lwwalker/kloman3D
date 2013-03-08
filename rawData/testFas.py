#!/usr/bin/python
import sys

lines = open(sys.argv[1], 'r').readlines()

dup = 0
mat = 0

n = {}
name = ""
for l in lines:
	l = l.strip("\n")
	if l != "":
		if l.find(">") >= 0:
			name = l.strip(">")
			print name

		else:
			if name in n.keys():
				dup += 1
				newS = l
				oldS = n[name]
				if newS == oldS:
					print "Found duplicate for {}, seqs match.".format(name)
					mat += 1
				else:
					print "Found duplicate for {}, seqs do not match.".format(name)
			else:
				n[name] = l
			
print "Found {} duplicates, of which {} matched".format(dup, mat)

