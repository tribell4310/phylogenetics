"""

concatSeqs.py

Take two fasta files ("1" and "2").  Concatentate the sequence in 1 onto the sequence in 2.  Obviously, fileheaders must match
exactly.  If there are mismatches, only output the mutually exclusive set.

"""


import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO


def main(inTarget, inQuery):
	# Read the input fastas in order
	seqContainerTarget, recordsTarget, recordLenTarget = intake(inTarget)
	seqContainerQuery, recordsQuery, recordLenQuery = intake(inQuery)

	# For every sequence in target, check for sequence in query
	recordsOut = []
	for item in recordsTarget:
		found_yet = False
		for possibleMatch in recordsQuery:
			if item.id == possibleMatch.id:
				# Found a match - concatente sequences and add 
				seq_list = [item.seq, possibleMatch.seq]
				new_seq = "".join([str(seq_rec) for seq_rec in seq_list])
				possibleMatch.seq = Seq(new_seq)
				recordsOut.append(possibleMatch)
				found_yet = True
				break
		if found_yet == False:
			print("Sequence not found: "+item.id)

	# Output the matching records to a new fast file
	with open(no_ext(inTarget)+"_"+no_ext(inQuery)+".fasta", "w") as output_handle:
		for i in range(0, len(recordsOut)):
			SeqIO.write(recordsOut[i], output_handle, "fasta")


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
	if len(sys.argv) == 3:
		main(str(sys.argv[1]), sys.argv[2])
	else:
		print("Check argument inputs: foo.py Fasta1 Fasta2")
		exit()