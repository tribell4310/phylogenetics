"""

iqtree_nodeLogos_subset.py

Convert IQTREE ASR output into human-readable outputs.
Use a comparative alignment to define only certain positions for analysis
Produce histograms and weblogos.

"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("Agg")
mpl.rc("figure", max_open_warning = 0)
import logomaker as lm
import statistics
import csv
from Bio import SeqIO


def main(inState, inCompareFasta, customNodes):
	# Check for the required output folders, create them if they do not exist.
	neededDirs = ['WebLogos', 'Histograms', 'FakeBFactors', 'RawOutputs']
	for neededDir in neededDirs:
		if os.path.exists(neededDir) == False:
			os.mkdir(neededDir)

	# Read in the fasta alignment, take all the names and make a list.  Discard these from analysis.
	#tipNodes = fastaIndeces(inFasta)

	# Load marginal probabilities file into a container object (indexed list of probabilities for every internal node)
	print("\nLoading data from IQ-TREE output statefile...")
	#probContainer = loadInternalNodes(inMargProb, tipNodes)
	with open(inState) as csv_file:
		probContainer = loadInternalNodes(csv_file)
	print("... "+str(len(probContainer))+" sequence positions imported.\n")

	# Identify the internal node labels, store in a list
	if customNodes != False:
		internalNodes = sorted(customNodes)
	else:
		internalNodes = sorted(probContainer[0].keys())

	# Add logic here to adjust internalNodes so it only contains the positions of interest
	# In other words, create a list defining the position of interest
	pos_of_interest = compare_parse(inCompareFasta)

	# For every internal node, pull out the appropriate pandas series at every position and
	# store as 2D DataFrames within a dictionary container.
	# Adjusted to take the pos_of_interest and only load subset positions into nodeDF_dict
	nodeDF_dict = {}
	for node in internalNodes:
		nodeDF_dict[node] = nodeToDataFrame(probContainer, pos_of_interest, node)

	# Create a MaxML version of each dataframe.
	maxlDF_dict, maxlTuple_dict, maxlStr_dict = extractionContainer(internalNodes, nodeDF_dict)

	# Write out the max likelihood values into a list
	bfactorContainer(internalNodes, maxlTuple_dict)

	# Output #1 - Weblogos for each of the nodes.
	weblogoContainer(internalNodes, nodeDF_dict, maxlDF_dict)

	# Output #2 - Posterior probability histograms for the internal nodes
	histogramContainer(internalNodes, maxlTuple_dict)

	# Output #3 - ML reconstructed sequences and overall probability
	summaryContainer(internalNodes, maxlStr_dict, maxlTuple_dict)

	# Output #4 - Raw outputs with max likelihood for each node
	rawOutputs(internalNodes, maxlStr_dict, maxlTuple_dict)

def compare_parse(inFasta):
	"""
	Returns a list of bools.  0 = no analyze; 1 = yes analyze
	"""

	with open(inFasta, "r") as f:
		records = SeqIO.parse(f, "fasta")
		counter = 0
		# Get the two sequences as strings.
		for record in records:
			counter += 1
			if counter == 1:
				seq1 = str(record.seq)
			else:
				seq2 = str(record.seq)

	# Check for sanitary input
	if len(seq1) != len(seq2):
		print("Error: comparative fasta sequences not the same length.  Exiting...")
		exit()

	# Compare the two seqs pos by pos to determine whether or not there is a "-" character
	compare_pattern = []
	for i in range(0, len(seq1)):
		if (seq1[i] != "-") and (seq2[i] != "-"):
			compare_pattern.append(1)
		else:
			compare_pattern.append(0)

	return compare_pattern




def bfactorContainer(internalNodes, tuple_dict):
	for node in internalNodes:
		g = open("FakeBFactors/"+str(node)+".txt", "w")
		relevant_set = list(tuple_dict[node])
		for i in range(0, len(relevant_set)):
			g.write(str(relevant_set[i][1])+"\n")
		g.close()


def maxlDataFrame(inDF, node):
	# Take the dataframe and iterate across every seq position to find the ML sequence
	# Format into:
	# ... a list of tuples containing (ML_residue, ML_prob)
	# ... the ML sequence as a string
	# ... a new data frame where every probaility is empty except for the ML sequence
	ml_tuples = []
	ml_str = ""
	ml_DF = inDF.copy()

	for i in range(0, len(inDF)):
		# Format into list of tuples
		ml_tuples.append((inDF.iloc[i].idxmax(), inDF.iloc[i][inDF.iloc[i].idxmax()]))

		# Format into a string
		ml_str += inDF.iloc[i].idxmax()

		# Format into a new ML dataframe
		maxIndex = inDF.iloc[i].idxmax()
		for residue in ml_DF.iloc[i].index:
			if residue != inDF.iloc[i].idxmax():
				ml_DF.iloc[i][residue] = 0.0

	return ml_tuples, ml_str, ml_DF

def drawWebLogo(inDict, node, colorScheme, outLabel, rangeStart, rangeEnd):
	# Condition plot size on the sequence range being analyzed (keeps the residue letters visible)
	if int(rangeEnd) - int(rangeStart) <= 50:
		myLogo = lm.Logo(inDict[node][int(rangeStart):int(rangeEnd)], color_scheme=colorScheme)
	else:
		myLogo = lm.Logo(inDict[node][int(rangeStart):int(rangeEnd)], color_scheme=colorScheme, figsize=[(len(inDict[node])/5.0), 2])
	myLogo.style_spines(visible=False)
	myLogo.ax.set_ylabel("Posterior Probability", labelpad=-1)
	myLogo.style_xticks(anchor=0, spacing=10, rotation=0)
	plt.savefig("WebLogos/"+node+"_"+outLabel+"_"+str(rangeStart)+"_"+str(rangeEnd)+".png")
	plt.close()

def nodeToDataFrame(probContainer, inPattern, node):
	# Extract one node's data from an indexed list of pandas series
	# Return it as a 2D pandas DataFrame
	#print(len(probContainer))
	#print(len(inPattern))
	total_pos_to_consider = 0
	for i in range(0, len(inPattern)):
		total_pos_to_consider += inPattern[i]

	# Check for length mismatch between statefile and comparative alignment
	if len(probContainer) != len(inPattern):
		print("Length of sites in statefile does not provide user-provided comparative alignment.  Please check your input files.\nExiting...")
		exit()

	nodeContainer = []
	for i in range(0, len(probContainer)):
		if inPattern[i] == 1:
			nodeContainer.append(probContainer[i][node])
	node_df = pd.DataFrame(nodeContainer, index=range(1, total_pos_to_consider+1))
	
	return node_df

def loadInternalNodes(inCsv):
	"""
	This function has been modified from the original version.
	Now reads an IQTREE ASR state output file and parses it into a dictionary.
	"""
	csv_reader = csv.reader(inCsv, delimiter='\t')
	residue_order = []
	nodeDict = {}
	static_array = []#This loads the csv into memory so I can make multiple loop passes after the csv iterator is gone.
	#first_line_read_flag = False
	# First, get the order in which residues are listed
	# And define the starting indeces for each node in a nodeDict object
	# And load the csv object into memory in static_array.
	row_counter = 0
	for row in csv_reader:
		if len(row) > 4:
			#for i in range(0, len(row)):
			static_array.append(row)
			

			if row_counter == 0:
				# Create a list containing the order of residues in the statefile
				for i in range(3, len(row)):
					residue_order.append(row[i][2]) #This is parsing out the one-letter AA code: "p_A" -> "A"
			else:
				if row[1] != "Site":
					if row[0] not in nodeDict:
						nodeDict[row[0]] = row_counter
			row_counter += 1
	
	#print(nodeDict)
	#print(residue_order)
	#print(static_array)

	# Extract a list of start positions from the nodeDict
	start_pos = []
	for item in nodeDict:
		start_pos.append(nodeDict[item])
	

	# Next, need to figure out the length of the alignment.
	# This can be figured by walking forward through rows until a site number duplicates.
	sites = []
	for i in range(1, row_counter):
		#print(static_array[i][1])
		if int(static_array[i][1]) not in sites:
			sites.append(int(static_array[i][1]))
		else:
			break
	#print("C")
	#print(sites)
	alignment_length = max(sites)
	#print(alignment_length)
		#if static_array[nodeDict[item]-1][1] != "329":
		#	print(static_array[nodeDict[item]-1])

	# Walk through the static array populating each alignment position for each node one at a time.
	indexedNodeProbs = []
	for i in range(0, alignment_length):
		nodePosDict = {}
		for pos in start_pos:
			# Go to that row, pull a 20-item list that is now linked ot resudie_order
			probs = static_array[pos+i][3:]
			node_id = static_array[pos+i][0]
			#print(node_id)
			#print(len(probs))

			# Populate a miniDict with it
			miniDict = {}
			for j in range(0, len(probs)):
				miniDict[residue_order[j]] = float(probs[j])
			#print(miniDict)
			#indexedNodeProbs.append()
			nodePosDict[node_id] = pd.Series(miniDict)
			#print(nodePosDict.keys())

		# Add all the data for that alignment position to the container list object
		indexedNodeProbs.append(nodePosDict)
	#for item in indexedNodeProbs:
	#	print(item.keys())
	#print(len(indexedNodeProbs))

	return indexedNodeProbs


	# Load
	#f = open(inProbFile, "r")
	#lines = f.readlines()

	# Find all of the starting positions
	#startPosList = [0]
	#for i in range(0, len(lines)):
	#	if lines[i] == "\n":
	#		if lines[i-1][0:8] == "of node:":
	#			startPosList.append(i+1)
	
	# Create an indexed container object to hold the pandas series for each position
	#indexedNodeProbs = []
	
	# For every position in the sequence, extract the corresponding block of lines, pass
	# it to the iterateOverOnePos to get the pandas series, and append it to the container.
	#for i in range(0, len(startPosList)-1):
	#	relevantLines = lines[startPosList[i]:startPosList[i+1]-1]
	#	indexedNodeProbs.append(iterateOverOnePos(relevantLines, tipNodes))
	
	#return indexedNodeProbs
	# The above code block needs to be spun into a seperate function.  This function
	# will be a container that parses the blocks of text and feeds input for parsing.


def iterateOverOnePos(lineBlock, tipNodes):
	# Initiate a dictionary to hold the probabilities for every node
	nodeProbs = {}

	# Skipping the first line (position number), make a dictionary of p() strings ONLY for
	# the internal nodes (exclude the tip nodes)
	for i in range(1, len(lineBlock)):
		firstColon = lineBlock[i].find(":")
		secondColon = lineBlock[i].find(":", firstColon+1)
		nodeLabel = lineBlock[i][firstColon+2:secondColon]
		if nodeLabel not in tipNodes:
			p_string = lineBlock[i][secondColon+2:]
			p_string.strip()
			nodeProbs[nodeLabel] = pd.Series(popProbDict(p_string))

	return nodeProbs

def popProbDict(p_string):
	# Standard amino acid dictionary
	std_AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

	# Populate a dictionary with float zeroes for every AA
	probDict = {}
	for AA in std_AAs:
		probDict[AA] = 0.0

	# Parse the p() string into floats, load into the dictionary
	p_list = p_string.split(" ")
	p_list = p_list[:-1]
	for value in p_list:
		key = value[2]
		prob = float(value[5:])
		probDict[key] = prob

	return probDict

def fastaIndeces(fastaFile):
	# Take a fasta file, return a list of index names
	# Load
	f = open(fastaFile, "r")
	lines = f.readlines()
	indeces = []

	# Find the indeces, add to the list
	for line in lines:
		if line[0] == ">":
			indeces.append(line[1:].strip())
	f.close()
	
	return indeces

def extractionContainer(nodeList, nodeDF_dict):
	print("Extracting ML sequences for internal nodes...")
	maxlDF_dict = {}
	maxlTuple_dict = {}
	maxlStr_dict = {}
	counter = 1

	for node in nodeList:
		print("... "+str(node)+" - node "+str(counter)+" of "+str(len(nodeList))+" ...")
		maxl_tuples, maxl_str, maxl_DF = maxlDataFrame(nodeDF_dict[node], node)
		maxlDF_dict[node] = maxl_DF
		maxlTuple_dict[node] = maxl_tuples
		maxlStr_dict[node] = maxl_str
		counter += 1
	print("... extraction complete.\n")

	return maxlDF_dict, maxlTuple_dict, maxlStr_dict

def weblogoContainer(nodeList, nodeDF_dict, maxlDF_dict):
	# Output the total probabilites in color and the ML probabilities in b/w
	print("Drawing WebLogos for internal node sequences...")
	counter = 1
	for node in nodeList:#Revert to internalNodes to restart the loop
		print("... "+str(node)+" - node "+str(counter)+" of "+str(len(nodeList))+" ...")
		# Use these lines to make weblogos for the full sequence
		drawWebLogo(nodeDF_dict, node, 'chemistry', 'allProbs', 0, len(nodeDF_dict[node]))
		drawWebLogo(maxlDF_dict, node, 'black', 'ML_seq', 0, len(nodeDF_dict[node]))

		# Use these lines for your own range - change 30, 71 to whatever range you want to analyze.
		#drawWebLogo(nodeDF_dict, node, 'chemistry', 'allProbs', 30, 71)
		#drawWebLogo(maxlDF_dict, node, 'black', 'ML_seq', 30, 71)

		counter += 1
		
	print("... WebLogo drawing complete.\n")


def h_index(inList):
	inList.sort()
	list_len = len(inList)
	proportional_increment = 1.0 / float(list_len)
	for i in range(0, len(inList)):
		target = 1.0 - (i * proportional_increment)
		value = target - inList[i]
		if value <= 0:
			return (value + inList[i])
	
	# Backup for failure (this shouldn't ever execute, but is here for edge cases.)
	if sort_list[0] < 0.01:
		return 0.0
	else:
		return 1.0


def histogramContainer(nodeList, tuple_dict):
	print("Drawing probability histograms for internal node sequences...")
	counter = 1
	for node in nodeList:#Revert to internalNodes to restart the loop
		print("... "+str(node)+" - node "+str(counter)+" of "+str(len(nodeList))+" ...")
		
		# Create a list containing all the probability values
		# THIS INCLUDES A SINGLE 0.0 PSEUDOCOUNT SO ALL THE HISTOGRAMS SCALE ACCURATELY
		allProbVals = [0.0]
		for j in range(0, len(tuple_dict[node])):
			allProbVals.append(tuple_dict[node][j][1])
		mean_pp = statistics.mean(allProbVals)
		#print(allProbVals[1:])
		h_ind = h_index(allProbVals[1:])
		
		# Plot values on a histogram
		n, bins, patches = plt.hist(allProbVals, 20, facecolor='blue', edgecolor='black', alpha=0.8)
		plt.title("ML Reconstruction for node "+str(node)+" - H-index = "+str(np.around(h_ind, decimals=3)))
		plt.xlabel("Posterior Probability")
		plt.ylabel("Counts")
		plt.savefig("Histograms/"+node+"_histogram.png")
		plt.close()

		counter += 1

	print("... histogram drawing complete.\n")

def summaryContainer(nodeList, string_dict, tuple_dict):
	print("Writing a final report for the max likelihood reconstructions...")
	f = open("NODE_SUMMARY.txt", "w")
	g = open("MAXL_NODES.fasta", "w")

	for node in nodeList:#Revert to internalNodes to restart the loop
		# Calculate mean posterior probability
		allProbVals = []
		out_seq = ""
		for j in range(0, len(tuple_dict[node])):
			allProbVals.append(tuple_dict[node][j][1])
			out_seq += tuple_dict[node][j][0]
		mean_pp = statistics.mean(allProbVals)

		# Write out to file
		# Change the print statements to f.write
		f.write("Internal Node: "+node+"\n")
		f.write(string_dict[node]+"\n")
		f.write("Mean Posterior Probability: "+str(mean_pp)+"\n\n")

		# Write fasta
		g.write(">"+node+"\n")
		g.write(out_seq+"\n")

	f.close()
	g.close()
	print("... summary complete.\n")

def rawOutputs(nodeList, string_dict, tuple_dict):
	print("Writing max post probs for all selected nodes...")
	
	# For each node, write out a csv containing the (AA, maxPP) tuple vals
	counter = 1
	for node in nodeList:
		print("... "+str(node)+" - node "+str(counter)+" of "+str(len(nodeList))+" ...")
		with open("RawOutputs/"+node+".csv", "w", newline="") as csvfile:
			csvwriter = csv.writer(csvfile)
			for j in range(0, len(tuple_dict[node])):
				csvwriter.writerow(tuple_dict[node][j])
		counter += 1

if __name__ == '__main__':	
	if len(sys.argv) == 3:
		main(sys.argv[1], sys.argv[2], False)
	elif len(sys.argv) > 3:
		main(sys.argv[1], sys.argv[2], sys.argv[3:])
	else:
		print("Check argument inputs: foo.py inState inCompareFasta [nodes]")
		exit()