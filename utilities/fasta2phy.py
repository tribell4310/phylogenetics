"""

fasta2phy.py
py3

Naively convert fasta to phylip format.

"""

import sys
import os
import string
from Bio import AlignIO
from Bio import SeqIO


def main(inFasta):
	fa2phy(inFasta)


def fa2phy(inFasta):
	# Do the naive Biopy conversion step.  This produces a file that doesn't meet phyML parsing standards
	myAlign = AlignIO.read(inFasta, "fasta")
	with open(inFasta[:-6]+"_TEMP_BAD.phy", "w") as handle:
		count = SeqIO.write(myAlign, handle, "phylip-sequential")

	# Re-load this file without Biopy to reformat
	f = open(inFasta[:-6]+"_TEMP_BAD.phy", "r")
	lines = f.readlines()
	lowers = string.ascii_lowercase
	collateLines = []
	for i in range(0, len(lines)):
		if lines[i][0] in lowers:
			collateLines.append(lines[i][0:10]+"\n")
			collateLines.append(lines[i][10:])
		else:
			collateLines.append(lines[i])
	f.close()

	# Write this file as the correct, re-formated file
	g = open(inFasta[:-6]+".phy", "w")
	for i in range(0, len(collateLines)):
		g.write(collateLines[i])
	g.close()

	# Remove the temporary file
	os.remove(inFasta[:-6]+"_TEMP_BAD.phy")


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()