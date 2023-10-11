"""

remove_seq_fuzzy.py

Takes a list of sequence labels and a fasta and uses fuzzy string matching to find the sequence header in the title.

Outpus a fasta file with the query sequences removed from the input fasta.

Prints to terminal when fuzzy matching is used, or any sequences that can't be matched.

"""


import sys
from Bio import SeqIO
from fuzzywuzzy import fuzz


def main(inFasta, inQuery):
	# Read the input fasta in order
	seqContainer = []
	records = list(SeqIO.parse(inFasta, "fasta"))
	recordLen = len(records)

	# For pretty formatting later, get the longest item in the list
	longest_id = 0
	for record in records:
		if len(no_slash(record.id)) > longest_id:
			longest_id = len(no_slash(record.id))

	# Read the inQuery into an object
	queries = []
	with open(inQuery, "r") as f:
		lines = f.readlines()
		for i in range(0, len(lines)):
			cleanStr = lines[i].replace(" ", "_").replace("\n", "")
			if len(cleanStr) > 0:
				queries.append(cleanStr)

	# First round: definite name matching.
	satisfiedQueries = []
	recordsToRemove = []

	for query in queries:
		for record in records:
			if query == no_slash(record.id):
				satisfiedQueries.append(no_slash(record.id))
				recordsToRemove.append(record)
				break

	# Update the pool of records that remain in the pool
	remaining_records = find_remainder(records, recordsToRemove)

	# Second round: fuzzy name matching
	if len(satisfiedQueries) == len(queries):
		print("\nAll sequences found without fuzzy string matching.")
		with open(no_ext(inFasta)+"_queriesRemoved.fasta", "w") as output_handle:
			for i in range(0, len(remaining_records)):
				SeqIO.write(remaining_records[i], output_handle, "fasta")
			print("\n...done.")
		exit()

	else:
		print("\n"+str(len(satisfiedQueries))+" of "+str(len(queries))+" were removed by definite string matching.\nBeginning fuzzy matching...\n\n"+rightpad("QUERY", longest_id+2)+rightpad("MATCHED RECORD", longest_id+2)+"SCORE")
		fuzz_records_to_remove = []
		for query in queries:
			if query not in satisfiedQueries:
				for record in remaining_records:
					score = fuzz.partial_ratio(query, record.id)
					if score > 80:
						fuzz_records_to_remove.append(record)
						satisfiedQueries.append(query)
						alert_out(query, record, score, longest_id+2)
						break

	# Update the pool of records that remain in the pool
	remaining_records_2 = find_remainder(remaining_records, fuzz_records_to_remove)

	# Third round: surrender
	if len(satisfiedQueries) == len(queries):
		print("\nAll sequences matched.")

	else:
		print("\nNot all sequences could be matched - "+str(len(queries)-len(satisfiedQueries))+" remain.")
		for query in queries:
			if query not in satisfiedQueries:
				print(query)

	# Write out
	with open(no_ext(inFasta)+"_queriesRemoved.fasta", "w") as output_handle:
		for i in range(0, len(remaining_records_2)):
			SeqIO.write(remaining_records_2[i], output_handle, "fasta")

	print("\n...done.")


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


def alert_out(query, record, score, pad):
	print(rightpad(str(query), pad)+rightpad(str(no_slash(record.id)), pad)+str(score))


def rightpad(inStr, inLen):
	while len(inStr) < inLen:
		inStr = inStr + " "
	return inStr


def find_remainder(records, recordsToRemove):
	# Define the records that remain after this
	remaining_records = []

	for record in records:
		found_flag = False
		for recordToRemove in recordsToRemove:
			if record.id == recordToRemove.id:
				found_flag = True
				break
		if found_flag == False:
			remaining_records.append(record)

	return remaining_records


def no_slash(inStr):
	return inStr[:inStr.find("/")]


if __name__ == '__main__':	
	if len(sys.argv) == 3:
		main(str(sys.argv[1]), str(sys.argv[2]))
	else:
		print("Check argument inputs: foo.py inFasta inQuery")
		exit()