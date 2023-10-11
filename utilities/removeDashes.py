"""
removeDashes.py
py3

"""

import sys


def main(inFasta):
	f = open(inFasta, "r")
	lines = f.readlines()
	
	collateLines = []
	for i in range(0, len(lines)):
		collateLines.append(remove_dash(lines[i]))

	g = open(inFasta[:-6]+"_dashesRemoved.fasta", "w")
	for line in collateLines:
		if line != "\n":
			g.write(line)
	g.write("\n")
	
	g.close()
	f.close()

def remove_dash(inLine):
	newLine = inLine.replace("-", "")
	return newLine

if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()