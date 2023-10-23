# phylogenetics

Various utilities to accelerate phylogenetic analysis written in python and R.  These scripts are supplied as-is, and you are welcome to use, modify, and distribute them with or without attribution.  I can't guarantee support for these scripts, but if you submit an issue, I'll do my best to have a look.

## utilities
 - autoshift.py
	 - Take two sequences and apply the gap pattern in sequence 1 to sequence 2.
 - concatSeqs.py
	 - Take two fasta files ("1" and "2").  Concatentate the sequence in 1 onto the sequence in 2.  Fileheaders must match exactly.  If there are mismatches, only output the mutually exclusive set.
 - extract_subtree_leaves.py
	 - User provides a tree, a fasta alignment, and a node label.  Extract the set of nodes in the tree downstream of the specified node and write them out as a new fasta.
 - fasta2phy.py
	 - Naively convert fasta to phylip format.
 - fastaCount.py
	 - Counts entries in a fasta.
 - findDupGenus.py
	 - Find sequences with duplicate genus.  Print to terminal.
 - getNames.py
	 - Take in a fasta file and output the species names to a text file in the format "Genus_species", "Genus_species", ...
 - ncbi2gs.py
	 - Rename NCBI nomenclature to Genus_species nomenclature.
 - nj_tree.py
	 - Python method for making a quick neighbor-joining tree from an aligned fasta file.
 - pseudoBL.py
	 - Takes a tree with no branch lengths (such as a synthetic guidetree) and adds a user-specified psuedo-length to every branch.
 - pyBlocks.py
	 - Removes alignment slivers and sections automatically using a user-specified completeness threshold.
 - remove_seq_fuzzy.py
	 - Takes a list of sequence labels and a fasta and uses fuzzy string matching to find the sequence header in the title. Outpus a fasta file with the query sequences removed from the input fasta. Prints to terminal when fuzzy matching is used, or any sequences that can't be matched.
 - removeBS.py
	 - Remove branch support values from a Newick tree.
 - removeDashes.py
	 - Removes all alignment gap dashes from a fasta file.
 - slashRemove.py
	 - Find and remove slash nomenclature from a fasta file or Newick tree.  The newick function could be easily adapted to work with any arbitrary string.  This script will break for trees that have sequences with duplicate names.
 - slashRemove_mac.py
	 - slashRemove.py, but modified to work on OSX systems.
 - subsample.py
	 - Randomly subsamples entries from a fasta file.

## auto-pruning
 - cutoffAnalysis.py
	 - Exhaustively compares all sequences in an alignment and makes a global profile of pairwise similarity.
 - pruneSimilar.py
	 - Analyzes sequences in a fasta file that has been sorted into a tree-order. If any two adjacent sequences are >NN percent identical, then the second one gets discarded. Outputs a new fasta suffix _idPruned.fasta
 - pruneSimilarGlobal.py
	 - Same as pruneSimilar.py, but all sequences are exhaustivley compared - slower, but more thorough.

## iqtree-scripts
 - iqtree_nodeLogos.py
	 - Convert IQTREE ASR output into human-readable outputs.
 - iqtree_state2fasta.py
	 - Reads an iqtree asr state file and converts it into a fasta with node labels.

## trees
 - basic_tanglegram.r
	 - Generate a tanglegram comparison between two newick trees with identical tip labels.
 - make_syn_tree.r
	 - Generate a synthetic tree for tanglegram comparison with an empirical phylogeny using the Open Tree of Life (OTOL) database.
 - root.r
	 - Root a tree against a set of outgroup tip labels.
 - tnrs_resolution.r
	 - Validation script to ensure all tip labels used in a phylogeny can be unambiguously matched to a known species in OTOL.

## modeller
 - See separate readme file in the modeller subdirectory.

> Written with [StackEdit](https://stackedit.io/).
