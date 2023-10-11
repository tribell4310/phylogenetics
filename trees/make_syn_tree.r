# make_syn_tree.r
#
# Makes a synthetic tree with tip labels matching the headers
# in the fasta file (Genus_species nomenclature).

# Imports
library(ape)
library(seqinr)
library(rotl)
library(stringr)
library(dendextend)

# User inputs.  If using OTOL-resolved sequences, inCSV can be left unchanged.
inCSV <- "psiblast_hmmer_v16_noSlash_aln_taxa.csv"
inFasta <- "psiblast_hmmer_v16_noSlash_aln.fasta"
outTree <- "psiblast_hmmer_v16_noSlash_aln_synTree.txt"

# Load the csv file and fasta files. 
raw_taxa <- read.csv(inCSV, sep=",")
seqs <- read.fasta(inFasta, seqtype="AA")
species <- getName(seqs)
lower_species <- tolower(species)

# Select the subset of taxa rows that are included in the fasta file.
rows_to_delete <- c()
for (i in 1:length(raw_taxa$search_string)) {
  if (raw_taxa$search_string[i] %in% lower_species) {} else {
    rows_to_delete <- c(rows_to_delete, i)
  }
}
taxa <- raw_taxa[-rows_to_delete, ]

# Induce a tree from the modified taxa file
tree <- tol_induced_subtree(ott_ids=taxa$ott_id, label_format="name")

# Convert relevant features to mutable lists
taxa_ss <- character()
taxa_ott_id <- character()
for (i in 1:length(taxa[["search_string"]])) {
  taxa_ss[i] <- taxa[["search_string"]][i]
  taxa_ott_id[i] <- taxa[["ott_id"]][i]
}
tree_tiplabel <- tolower(tree$tip.label)
taxa_ott_id <- unname(taxa_ott_id)

# Some of the labels come back as:
#   Query: "Homer_simpson"
#   Response: "Homer_simpson (Species in kingdom cartoon)"
# Delete the "(Species in kingdom cartoon)" part in tree_tiplabel
for (i in 1:length(tree_tiplabel)) {
  myStr <- tree_tiplabel[i]
  for (j in 1:nchar(myStr)) {
    if (substring(myStr, j, j) == "(") {
      newStr <- substring(myStr, 1, j-2)
      tree_tiplabel[i] <- newStr
      break
    }
  }
}

# Match tree_tiplabel (lowercase) to taxa_ss
matches <- character()
notMatches <- character()
for (val in tree_tiplabel) {
  if (val %in% taxa_ss) {
    matches <- c(matches, val)
  } else {
    notMatches <- c(notMatches, val)
  }
}

# Use matches to identify the initial query taxa that correspond to notMatches
notMatchPool <- character()
for (val in taxa_ss) {
  if (val %in% matches) {
    # do nothing, I'm too lazy to negate the %in% function
  } else {
    notMatchPool <- c(notMatchPool, val)
  }
}

# If notMatches is not empty, alert the user
if (length(notMatchPool != 0)) {
  print("Caution!  The following sequences aren't exact matches and will break tangelgrams.")
  for (val in notMatchPool) {
    print(val)
  }
}

# Change back from "genus" to "Genus" nomenclature
tree_tips_repaired_title <- character()
for (i in 1:length(tree_tiplabel)) {
  tree_tips_repaired_title[i] <- str_to_title(tree_tiplabel[i])
}

# Put the repaired names into the tree object
tree[["tip.label"]] <- tree_tips_repaired_title

# Write out a synthetic newick tree with labels that match the taxa object.
write.tree(tree, file = outTree)
