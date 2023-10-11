"""

nj_tree.py

Python method for making a quick neighbor-joining tree from an aligned fasta file.

"""

import sys
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo
from Bio import AlignIO


def main(inAln):
	# Read the alignment file - hardcoded for fasta cause it's what i use...
	my_aln = AlignIO.read(open(inAln), "fasta")

	# Build a NJ tree
	print("Calculating distance")
	constructor = DistanceTreeConstructor()
	calculator = DistanceCalculator('identity')
	dm = calculator.get_distance(my_aln)
	print("Constructing tree")
	njtree = constructor.nj(dm)
	print("Done.")

	# Write out
	with open(no_ext(inAln)+".newick", "w") as handle:
		count = Phylo.write(njtree, handle, "newick")



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


if __name__=="__main__":
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check inputs: foo.py inAlign.fasta")