"""

automodel_setup.py

Converts a fasta file into a set of .ali files for modeller, with names based on fasta headers.
Output lines need to be all-caps and 75 residues long max.

"""

import sys
import os
from Bio import SeqIO


def main(inFasta):
	# Check for the required output folders, create them if they do not exist.
	neededDirs = ["ali_files"]
	for neededDir in neededDirs:
		if os.path.exists(neededDir) == False:
			os.mkdir(neededDir)

	# Read the input fasta in order
	seqContainer = []
	records = list(SeqIO.parse(inFasta, "fasta"))
	recordLen = len(records)

	for record in records:
		f = open("ali_files/"+str(record.id.replace("|", "_"))+".ali", "w")
		f.write(">P1;"+str(record.id)+"\nsequence:"+str(record.id)+":: :: ::: 0.00: 0.00\n")
		seq_lines = chunker(record.seq)
		for i in range(0, len(seq_lines)):
			f.write(seq_lines[i]+"\n")
		f.close()


def chunker(inSeq):
	working_seq = inSeq.upper()+"*"
	out_lines = []
	temp_seq = ""
	for i in range(1, len(working_seq)+1):
		temp_seq += working_seq[i-1]
		if i % 75 == 0:
			out_lines.append(temp_seq)
			temp_seq = ""
	if len(temp_seq) > 0:
		out_lines.append(temp_seq)

	return out_lines


if __name__ == '__main__':	
	if len(sys.argv) == 2:
		main(sys.argv[1])
	else:
		print("Check argument inputs: foo.py inFasta")
		exit()