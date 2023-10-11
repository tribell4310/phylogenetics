"""

slashRemove_os.py
py3

Find and remove slash nomenclature from a fasta file or Newick tree.  The newick function could be 
easily adapted to work with any arbitrary string.  WARNING: This script will break for trees that 
have sequences with duplicate names.

Gets rid of problems in the sequence too.  Random backslashes

"""

import sys
import string


def main(inFile):
	# Detect whether inFile is a fasta or a Newick tree
	f = open(inFile, "r")
	lines = f.readlines()

	if lines[0][0] == ">":
		# It's a fasta file
		newLines = changeFasta(lines)
		g = open(inFile[:-6]+"_noSlash.fasta", "w")
	elif lines[0][0] == "(" and lines[0][-1] == ";":
		# It's not a fasta file, presumably a newick tree
		newLines = changeTree(lines)
		g = open(inFile[:-4]+"_noSlash.txt", "w")
	else:
		# If we've made it here, this is a badly-formatted input.  Alert the user and exit.
		print("This doesn't appear to be a fasta file or newick tree.  Exiting...")
		exit()

	# Write out a new file
	for line in newLines:
		g.write(line)

	# Clean up
	g.close()
	f.close()
	

def changeFasta(lines):
	# Re-print the fasta titles excluding the slash and everything after
	collateLines = []
	for i in range(0, len(lines)):
		if lines[i][0] == ">":
			slashLoc = lines[i].find("/")
			collateLines.append(lines[i][0:slashLoc].replace("\\", "")+"\n")
		else:
			collateLines.append(lines[i].replace("\\", ""))
	return collateLines


def changeTree(lines):
	# Scan the newick string for words, then replace the string
	allLetters = string.ascii_letters

	# Populate a dictionary that references every instance of a sequence label to its version without a slash
	matchingDict = {}

	stopFlag = False
	for i in range(0, len(lines[0])):
		if stopFlag == False:
			if lines[0][i] in allLetters:
				colonLoc = lines[0][i:].find(":")
				matchingDict[lines[0][i:i+colonLoc]] = naiveRemoveSlash(lines[0][i:i+colonLoc])
				stopFlag = True
		if lines[0][i] == ":":
			stopFlag = False

	# Use the dictionary to perform replacement.
	workingStr = lines[0]
	for item in matchingDict:
		workingStr = workingStr.replace(item, matchingDict[item])
			
	return [workingStr]


def naiveRemoveSlash(inStr):
	slashLoc = inStr.find("/")
	return inStr[0:slashLoc]


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(str(sys.argv[1]))
	else:
		print("Check argument inputs: foo.py inFile.  inFile must be a fasta file or Newick tree.")
		exit()
