"""

getNames.py

Take in a fasta file and output the species names to a text file in the format "Genus_species", "Genus_species", ...

"""

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO


def main(inFasta):
	# Read the input fastas in order
	seqContainer, records, recordLen = intake(inFasta)

	# For every sequence in fasta, save id
	names = []
	for sequence in records:
		names.append(sequence.id)
	#print(names)

	# Write out in desired format
	outgroup = ""
	for name in names:
		outgroup = addToString(outgroup, name)
	outgroup = outgroup[:-2]
	#print(outgroup)
	g = open(no_ext(inFasta)+"_seqNames.txt", "w")
	g.write(outgroup)
	g.close()


def addToString(inStr, addPart):
	return inStr+"\""+addPart+"\", "
	

def intake(inFasta):
	# Read the input fasta in order
	seqContainer = []
	records = list(SeqIO.parse(inFasta, "fasta"))
	recordLen = len(records)
	return seqContainer, records, recordLen


def no_ext(inStr):
	"""
	Takes an input filename and returns a string with the file extension removed.

	"""
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find(".", prevPos+1)
	return inStr[0:prevPos]


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(str(sys.argv[1]))
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()