"""

iqtree_state2fasta.py

Reads an iqtree asr state file and converts it into a fasta with node labels.

"""

import sys
import csv


def main(inState):
	with open(inState) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter='\t')
		dataDict = {}
		for row in csv_reader:
			if len(row) > 4:
				if row[0] not in dataDict:
					dataDict[row[0]] = ""
				dataDict[row[0]] += row[2]

	write_to_fasta(dataDict, inState)


def write_to_fasta(seqDict, inName):
	collect = []
	for label in seqDict:
		collect.append(">"+label.replace(" ", "_"))
		collect.append(seqDict[label])

	g = open(inName+"_MLseqs.fasta", "w")
	for item in collect:
		g.write(item+"\n")
	g.close()


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(str(sys.argv[1]))
	else:
		print("Check argument inputs: foo.py stateFile")
		exit()