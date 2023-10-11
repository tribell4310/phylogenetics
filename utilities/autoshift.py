"""

autoshift.py

Take two sequences and apply the gap pattern in sequence 1 to sequence 2.

"""

import sys
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2


def main(inFasta):
	# Read the two items
	headers = []
	seqs = []
	with open(inFasta, "r") as f:
		records = SeqIO.parse(f, "fasta")

		# Do they still have the gaps?
		for record in records:
			headers.append(str(record.id))
			seqs.append(str(record.seq))

		alignment = pairwise2.align.globalxx(seqs[0], seqs[1])
		aln1 = str(alignment[0].seqA)
		aln2 = str(alignment[0].seqB)
		print(" ")

	with open(no_ext(inFasta)+"_aln.fasta", "w") as g:
		g.write(">"+headers[0]+"\n"+aln1+"\n>"+headers[1]+"\n"+aln2+"\n")	


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


if __name__ == "__main__":
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check usage: foo.py inFasta")