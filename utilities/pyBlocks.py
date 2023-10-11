"""

pyBlocks.py

Removes alignment slivers and sections automatically using a user-specified completeness threshold.

"""


import sys


def main(inFasta, userLimit):
	# Read the fasta
	headers, seqs = fileIntake(inFasta)
	clean_headers = []
	clean_seqs = []
	for i in range(0, len(headers)):
		clean_headers.append(headers[i].replace("\n", ""))
		clean_seqs.append(seqs[i].replace("\n", ""))

	# Verify that all seqs are same length (i.e. this is an alignment)
	for i in range(1, len(seqs)):
		if len(seqs[i]) != len(seqs[i-1]):
			print("Input must be formatted as a fasta alignment.  Exiting...")
			exit()

	# Work through column by column and determine whether it should be preserved in the final file.
	keepList = [] # This list keeps track of which alignement positions pass
	numSeqs = len(clean_seqs)
	for i in range(0, len(clean_seqs[0])):
		counter = 0
		for j in range(0, numSeqs):
			if clean_seqs[j][i] != "-":
				counter += 1
		if float(counter) / float(numSeqs) >= userLimit:
			keepList.append(1) # Keep this position
		else:
			keepList.append(0) # Exclude this position

	# Reconstruct the sequences using only the positions that passed threshholding
	new_seqs = []
	for i in range(0, len(clean_seqs)):
		new_seqs.append("")

	for i in range(0, len(keepList)):
		if keepList[i] == 1:
			for j in range(0, len(new_seqs)):
				new_seqs[j] += clean_seqs[j][i]

	# Write out
	fastaWriter(headers, new_seqs, inFasta, "pyBlocks")

	# Report out
	print("\nInput alignment contains "+str(len(clean_seqs[0]))+" columns.")
	print("Output alignment contains "+str(len(new_seqs[0]))+" columns.")



def fileIntake(inFasta):
	f = open(inFasta, "r")
	headers, seqs = customParse(f.readlines())
	f.close()

	return headers, seqs


def fastaWriter(headers, seqs, inFasta, suffix):
	g = open(no_ext(inFasta)+"_"+suffix+".fasta", "w")
	for i in range(0, len(headers)):
		g.write(">"+headers[i]+"\n")
		g.write(seqs[i]+"\n")
	g.close()


def customParse(inLines):
	"""
	Split fasta rawlines into headers and sequences.

	"""
	# Iterate through the list and find headers by their leading carat
	seqFlag = False
	accumSeq = ""
	headers = []
	seqs = []
	for i in range(0, len(inLines)):
		if inLines[i][0] == ">":
			if seqFlag == True:
				seqs.append(accumSeq.replace("\n", ""))
				accumSeq = ""
				seqFlag = False
			headers.append(inLines[i][1:].replace("\n", ""))
		else:
			seqFlag = True
			accumSeq += inLines[i]
			# Catch the final sequence in the fasta file and ensure it's written
			if i == len(inLines) - 1:
				#print(accumSeq)
				seqs.append(accumSeq.replace("\n", ""))
				accumSeq = ""
	
	# Don't allow to progress if parsing fails
	if len(headers) != len(seqs):
		raise Exception("Invalid fasta input - not the same number of sequences as sequence headers.")

	return headers, seqs


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


def fastaCountFromFile(inFile):
	# Read the file
	f = open(inFile, "r")
	lines = f.readlines()
	f.close()

	# Count the lines
	myCounter = 0
	for line in lines:
		if line[0] == ">":
			myCounter += 1

	return myCounter


if __name__ == "__main__":
	if len(sys.argv) == 3:
		main(sys.argv[1], float(sys.argv[2]))
	else:
		print("Check argument inputs: foo.py inFasta threshold")
		exit()