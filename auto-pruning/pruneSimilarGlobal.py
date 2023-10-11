"""

pruneSimilarGlobal.py

 - This script analyzes sequences in a fasta file without assuming they are in tree order (i.e. exhaustive O(mn) comparison).
 - If any two adjacent sequences are >NN percent identical, then the second one gets discarded.
 - Outputs a new fasta suffix _idPruned.fasta
 - Optimized for Python 3

 - Now finally incorporates optional whitelisting!  Checks for a file called "passlist.txt" in the cwd.
    - Includes legacy support for "whitelist.txt"

"""


import sys
import os
from Bio import SeqIO
from Bio import pairwise2


def main(inFasta, idLimit):
	# Read the input fasta in order
	seqContainer = []
	records = list(SeqIO.parse(inFasta, "fasta"))
	#for protein in records:
	#	seqContainer.append(protein)
	recordLen = len(records)

	# Check for the existance of a file called "passlist.txt".  If it exists, load the contents into a file.
	# Includes legacy support for the deprecated term "whitelist"
	white_flag = False
	if os.path.isfile("passlist.txt") == True:
		print("\nPasslist detected with the following sequences:")
	elif os.path.isfile("whitelist.txt") == True:
		print("\nPasslist detected as whitelist.txt.\nWhitelist is a deprecated term, please consider converting to passlist.\nThe following sequences were detected:")
		white_flag = True
	else:
		g = open("passlist.txt", "w")
		g.close()
		
	if white_flag == True:
		f = open("whitelist.txt", "r")
	else:
		f = open("passlist.txt", "r")
	passlist = f.readlines()
	f.close()
	clean_passlist = []
	for item in passlist:
		clean_passlist.append(item.strip())
	for item in clean_passlist:
		print(item)
	print(" ")

	# Compare each item to its previous entry, load keep (True) / delete (False) info into a flagList.
	# By default, always keep the first entry.
	flagList = []

	for i in range(0, len(records)):
		print("Scanning sequence "+str(i+1)+" of "+str(recordLen)+"...")
		if checkPasslist(clean_passlist, records[i].id) == True:
			flagList.append(True)
		else:
			continueFlag = True
			for j in range(i+1, len(records)):
				if continueFlag == True:
					if alignPairs(records[i].seq, records[j].seq) > idLimit:
						flagList.append(False)
						print("Discard")
						continueFlag = False
			if continueFlag == True:
				flagList.append(True)

	# Go through the list and copy the "keeper" sequences to a new record
	recordContainer = []
	for i in range(0, len(flagList)):
		if flagList[i] == True:
			recordContainer.append(records[i])
	newRecordLen = len(recordContainer)

	# Write out
	with open(inFasta[:-6]+"_idPrunedGlobal_"+str(100*idLimit)+".fasta", "w") as output_handle:
		for i in range(0, len(recordContainer)):
			SeqIO.write(recordContainer[i], output_handle, "fasta")

	# Report back
	print("Done.  Input threshhold was "+str(idLimit)+" percent identity.")
	print("Output contains "+str(newRecordLen)+" of "+str(recordLen)+" original entries.")


def alignPairs(seq1, seq2):
	alignment = pairwise2.align.globalxx(seq1, seq2)
	return sorted(alignment)[0][2] / max(len(seq1), len(seq2))

def checkPasslist(passlist, iden):
	# Take the id and convert it to Genus nomenclature
	scoreLoc = iden.find("_")
	if scoreLoc == -1:
		genus = iden
	else:
		genus = iden[:scoreLoc]

	# Check the passlist
	is_in_list = False
	for item in passlist:
		if item.capitalize().strip() == genus.capitalize().strip():
			is_in_list = True
			print("Found in passlist: "+genus)
			break

	return is_in_list

if __name__ == '__main__':	
	if len(sys.argv) == 3:
		main(str(sys.argv[1]), float(sys.argv[2]))
	else:
		print("Check argument inputs: foo.py inFasta identityLimit\nThis program supports passlisting.\nTo exclude sequences from parsing, include the genus in Genus nomenclature on separate lines in the file passlist.txt in this directory.")
		exit()