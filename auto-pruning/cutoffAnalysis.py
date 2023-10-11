"""
cutoffAnalysis.py
py3

Exhaustively compares all sequences in an alignment and makes a global profile of pairwise similarity.
Derived from the pruneSimilarGlobal.py script.

"""


import sys
import os
from Bio import SeqIO
from Bio import pairwise2
import matplotlib
import matplotlib.pyplot as plt


def main(inFasta):
	# Read the input fasta in order
	seqContainer = []
	records = list(SeqIO.parse(inFasta, "fasta"))
	recordLen = len(records)

	# Compare each item to its previous entry, load keep (True) / delete (False) info into a flagList.
	# By default, always keep the first entry.
	flagList = []
	cutoffList = []
	for i in range(1, 10000):
		flagList.append(0)
		cutoffList.append(i/10000.)

	for i in range(0, len(records)):#for every record
		print("Scanning sequence "+str(i+1)+" of "+str(recordLen)+"...")
		similarity = 0.0
		for j in range(i+1, len(records)):#compare to every subsequent record
			# Find the highest similarity with any other sequence
			local_similarity = alignPairs(records[i].seq, records[j].seq)
			if local_similarity > similarity:
				similarity = local_similarity

		# Record whether this sequence would have passed each cutoffvalue
		for k in range(0, len(cutoffList)):# for each cutoffvalue
			if similarity < cutoffList[k]:
				flagList[k] += 1

	# Print results
	g = open(inFasta[:-6]+"_cutoffAnalysis.csv", "w")
	for i in range(0, len(cutoffList)):
		g.write(str(cutoffList[i])+","+str(flagList[i])+"\n")
	g.close()

	# Plot
	matplotlib.use("Agg")
	plt.plot(cutoffList, flagList)
	plt.xlabel("Pairwise identity threshhold")
	plt.ylabel("Surviving sequences")
	plt.title("Auto-prune cutoff analysis\n"+inFasta[:-6])
	plt.savefig(inFasta[:-6]+".png")


def alignPairs(seq1, seq2):
	alignment = pairwise2.align.globalxx(seq1, seq2)
	return sorted(alignment)[0][2] / max(len(seq1), len(seq2))


def checkWhitelist(whitelist, iden):
	# Take the id and convert it to Genus nomenclature
	scoreLoc = iden.find("_")
	if scoreLoc == -1:
		genus = iden
	else:
		genus = iden[:scoreLoc]

	# Check the whitelist
	is_in_list = False
	for item in whitelist:
		if item.capitalize().strip() == genus.capitalize().strip():
			is_in_list = True
			print("Found in whitelist: "+genus)
			break

	return is_in_list


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(str(sys.argv[1]))
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()