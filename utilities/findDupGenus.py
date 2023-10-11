"""

findDupGenus.py

Find sequences with duplicate genus.  Print to terminal.

"""


import sys


def main(inFasta):
	f = open(inFasta, "r")
	lines = f.readlines()
	
	titleList = []

	# Find the list indeces of the duplicate sequences.
	for i in range(0, len(lines)):
		if lines[i][0] == ">":
			# Find the _
			usLoc = lines[i].find("_")
			titleList.append(lines[i][1:usLoc])

	titleList = sorted(titleList)

	for i in range(0, len(titleList)):
		if titleList[i] in titleList[i+1:]:
			print(titleList[i])


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(str(sys.argv[1]))
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()