"""

extract_subtree_leaves.py

User provides a tree, a fasta alignment, and a node label.  Extract the set of 
nodes in the tree downstream of the specified node and write them out as a new fasta.

"""


import sys
from Bio import Phylo
from Bio import SeqIO


def main(inFasta, inTree, inQuery):
	# Read the fasta into a sequence object
	records = SeqIO.parse(open(inFasta, "r"), "fasta")

	# Read the tree into one big string
	with open(inTree, "r") as tree_raw:
		tree = tree_raw.read()

	# Sanitize the input query
	query = inQuery.strip()

	# Search the tree for the user query label
	node_loc = tree.find(query)
	if node_loc == -1:
		print(query+" not found.  Exiting...")
		exit()
	else:
		# Extract the substring where the parentheses are not matched
		subtree = get_subtree(tree, node_loc)

	# Identify all the internal labels
	leaves = get_internal_labels(subtree)

	# Search the fasta file for the leaf labels, write out to new fasta
	with open(no_ext(inFasta)+"_subtreeItems.fasta", "w") as g:
		for record in records:
			if record.id.replace("|", "_") in leaves:
				g.write(">"+str(record.id)+"\n"+str(record.seq)+"\n")
	print("\t...done.")


def get_subtree(tree, start):
	# Look both ways, find the way where the paren opens AWAY from the node label
	counter = 1
	direction_to_proceed = 0
	while direction_to_proceed == 0:
		# Get the string characters in tree at equal distances away from start point
		back_assess = tree[start-counter]
		front_assess = tree[start+counter]
		if back_assess == ")":
			direction_to_proceed = -1
			break
		elif front_assess == "(":
			direction_to_proceed = 1
			break
		counter += 1

	# Walk in that direction keeping track of the parentheses
	# The loop has re-closed when the parentheses are all matched up.
	paren_count = 0
	override = True
	if direction_to_proceed == 1:
		while paren_count != 0 or override == True:
			override = False
			#print(paren_count, tree[start+counter])
			if tree[start+counter] == "(":
				paren_count += 1
			elif tree[start+counter] == ")":
				paren_count -= 1
			counter += 1
		return tree[start:start+counter+1]

	elif direction_to_proceed == -1:
		while paren_count != 0 or override == True:
			override = False
			#print(paren_count, tree[start-counter])
			if tree[start-counter] == ")":
				paren_count += 1
			elif tree[start-counter] == "(":
				paren_count -= 1
			counter += 1
		return tree[start-counter+1:start]


def get_internal_labels(tree):
	leaves = []
	counter = 0
	# Pattern to search for (x: or ,x: where x=label
	for i in range(0, len(tree)):
		if tree[i] in ["(", ","]:
			# Find the next colon
			colon_loc = tree.find(":", i+1)
			if colon_loc == -1:
				break
			else:
				leaves.append(tree[i+1:colon_loc])
				counter += 1
	# Resolve duplicate string + paren problems
	clean_leaves = []
	for leaf in leaves:
		new_str = ""
		for i in range(0, len(leaf)):
			if leaf[i] not in ["(", ")"]:
				new_str += leaf[i]
		clean_leaves.append(new_str)
	clean_unique_leaves = set(clean_leaves)
	print("Number of leaves found:\t"+str(len(clean_unique_leaves)))
	
	return leaves


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


if __name__=="__main__":
	if len(sys.argv) == 4:
		main(sys.argv[1], sys.argv[2], sys.argv[3])
	else:
		print("Check usage: foo.py inFasta inTree inQuery")
