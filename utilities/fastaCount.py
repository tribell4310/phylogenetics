"""

fastaCount.py

Counts entries in a fasta.

"""


import sys


def main(inFasta):
	f = open(inFasta, "r")
	lines = f.readlines()
	myCounter = 0
	for line in lines:
		if line[0] == ">":
			myCounter += 1

	print(str(myCounter)+" entries in "+str(inFasta))


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print ("Check argument inputs: foo.py inFasta")
		exit()