"""

pseudoBL.py
py3

Takes a tree with no branch lengths (such as ROTL synthetic guidetree) and adds a user-specified psuedo-length to every branch.

"""

import sys


def main(inTree, pseudoLength):
	# Read the tree
	f = open(inTree, "r")
	lines = f.readlines()
	myTree = lines[0]

	# Insert branch lengths.
	myTree = myTree.replace(",", ":"+pseudoLength+",")
	myTree = myTree.replace(")", ":"+pseudoLength+")")

	# Write out
	g = open(inTree[:-4]+"_pseudoBLs.txt", "w")
	g.write(myTree)
	g.close()
	f.close()


if __name__ == '__main__':	
	if len(sys.argv) == 3:
		main(str(sys.argv[1]), str(float(sys.argv[2])))
	else:
		print("Check argument inputs: foo.py inTree pseudoLength")
		exit()