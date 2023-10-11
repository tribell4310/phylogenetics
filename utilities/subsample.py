"""

subsample.py
py3

Randomly subsamples entries from a fasta file.

"""


from Bio import SeqIO
import sys
import random


def main(inFasta, inSize):
	# Get an initial count of the sequences
	numSeqs = fastaCount(inFasta)
	print("\nInput fasta has "+str(numSeqs)+" sequences.")
	if inSize <= numSeqs:
		print("Randomly subsampling "+str(inSize)+" sequences from this pool.")
	else:
		print("Error: user requested more sequences than are present in the input fasta.  Exiting...")
		exit()

	# Load fasta sequences into a container list using Biopy seq objects
	records = list(SeqIO.parse(inFasta, "fasta"))

	# Make a randomized list of 0/1 keys.
	keyList = []
	for i in range(0, inSize):
		keyList.append(1)
	for i in range(0, numSeqs-inSize):
		keyList.append(0)
	random.shuffle(keyList)

	# Append seqs to newList
	recordContainer = []
	for i in range(0, numSeqs):
		if keyList[i] == 1:
			recordContainer.append(records[i])

	# Write out to a new fasta structure
	with open(inFasta[:-6]+"_subsampled.fasta", "w") as output_handle:
		SeqIO.write(recordContainer, output_handle, "fasta")
	print("... Done.")


def fastaCount(inFasta):
	f = open(inFasta, "r")
	lines = f.readlines()
	myCounter = 0
	for line in lines:
		if line[0] == ">":
			myCounter += 1
	f.close()
	return myCounter


if __name__ == '__main__':	
	if len(sys.argv) == 3:
		main(sys.argv[1], int(sys.argv[2]))
	else:
		print("Check argument inputs: foo.py inFasta numberSeqs")
		exit()