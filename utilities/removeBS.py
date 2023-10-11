"""

removeBS.py
py3

Remove branch support values from a Newick tree so it can pass FastML parsing.

"""

import sys


def main(inTree):
	noBS(inTree)


def noBS(inTree):
	#Open the treefile
	f = open(inTree, "r")
	treevals = f.read().strip()

	# Delete all floating point numbers that appear between an end paren and a colon.
	# Find those indeces and their values, and store them in a dictionary.  Then do a trival str find/replace cycle.
	corrDict = {}
	for i in range(0, len(treevals)):
		if treevals[i] == ")":
			colonLoc = treevals.find(":", i)
			corrDict[i] = treevals[i:colonLoc+1]

	# Sanitize the dictionary to remove empty listings
	cleanDict = {}
	for item in corrDict:
		if corrDict[item].strip() != "":
			cleanDict[item] = corrDict[item]

	# Srting replacement cycle
	for item in cleanDict:
		newStr = treevals.replace(cleanDict[item], "):", 1)
		treevals = newStr

	# Treevals no longer contains branch supports.  Write to file.
	g = open(inTree[:-4]+"_noBS.txt", "w")
	g.write(treevals)
	g.close()
	f.close()


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check argument inputs: foo.py inTree")
		exit()