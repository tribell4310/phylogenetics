# 210420_tnrs_resolution.R
# This script is designed to assist with semi-automatic resolution of
# species against OTOL TNRS.  It can be used when initiating a new repository,
# but should not be necessary once the taxa table is fully validated.
# Thus, TNRS resolution is separated into a one-time task that produces an
# object to be used for specific tree generation.

# Dependencies
library(ape)
library(seqinr)
library(rotl)
library(stringr)
library(dendextend)

# Read in seqs, species, taxa YOUR FILENAME HERE!!!
seqs <- read.fasta("psiblast_hmmer_v41_noSlash_testPrank_aln.fasta", seqtype="AA")
outTreeName <- "psiblast_hmmer_v41_noSlash_testPrank_aln_synTree.txt"
species <- getName(seqs)
taxa <- tnrs_match_names(species)
tree <- tol_induced_subtree(ott_ids=taxa$ott_id, label_format="name")
# write.csv(taxa, "psiblast_hmmer_v21_noSlash_aln_taxa.csv")
# taxa <- read.csv("psiblast_hmmer_v21_noSlash_aln_taxa.csv", sep=",")

# Convert relevant lists to 1D vectors
taxa_ss <- character()
taxa_ott_id <- character()
for (i in 1:length(taxa[["search_string"]])) {
  taxa_ss[i] <- taxa[["search_string"]][i]
  taxa_ott_id[i] <- taxa[["ott_id"]][i]
}
tree_tiplabel <- tolower(tree$tip.label)
taxa_ott_id <- unname(taxa_ott_id)

# I now have three de-named 1D vectors with the correlated array indeces

# Pre-processing step - some of the labels come back as
# Query: "Homer"
# Response: "Homer (Genus in simpson)"
# Just delete the "(Genus in simpson)" part in tree_tiplabel
#print(tree_tiplabel)
#for (i in 1:length(tree_tiplabel)) {
#  myStr <- tree_tiplabel[i]
#  for (j in 1:nchar(myStr)) {
#    if (substring(myStr, j, j) == "_") {
#      newStr <- substring(myStr, 1, j-1)
#      tree_tiplabel[i] <- newStr
#      break
#    }
#  }
#}

# For every item in the tree, check to see if it has an exact match in
# taxa_ss.  If not, then flag it for repair.

#toRepair <- character()
#for (val in tree_tiplabel) {
#  foundFlag <- FALSE
#  for (j in 1:length(taxa_ss)) {
#    if (taxa_ss[j] == val) {
#      foundFlag <- TRUE
#      break
#    }
#  }
#  if (foundFlag == FALSE) {
#    toRepair <- c(toRepair, val)
#  }
#}
#print(toRepair)

# Now - search for matches between tree_tiplabel (lowercase) and taxa_ss
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

# Now need to correlate the notMatches and notMatchPool list
# Can do this by re-doing the tol query against items in notMatchPool
# Need ott id's.  Go back to taxa to get this.
# Find the items in notMatchPool in taxa_ss
notMatchNodeIDs <- character()
this_counter = 0
#print(length(notMatchPool))
#print(notMatchPool)
for (val in notMatchPool) {
  #print(val)
  index <- match(val, taxa_ss)
  match_ott_id <- taxa_ott_id[index]
  match_tni <- rotl::tol_node_info(match_ott_id)
  #print(match_tni)
  this_counter = this_counter + 1
  #print(this_counter)
  #Sys.sleep(2)
  
  if (substring(match_tni[["node_id"]], 1, 3) == "ott") {
    notMatchNodeIDs <- c(notMatchNodeIDs, tolower(match_tni[["taxon"]]$unique_name))
  } else {
    notMatchNodeIDs <- c(notMatchNodeIDs, match_tni[["node_id"]])
  }
}

# We now have a correlated list of correct and incorrect names
# Now to repair the tree tip labels!
#I'll name a new list called tree_tips_repaired
tree_tips_repaired <- tree_tiplabel
for (i in 1:length(notMatchNodeIDs)) {
  myLoc <- match(notMatchNodeIDs[i], tree_tiplabel)
  tree_tips_repaired[myLoc] <- notMatchPool[i]
}
#print(tree_tiplabel)
#print(tree_tips_repaired)

# Change back from "genus" to "Genus" nomenclature
tree_tips_repaired_title <- character()
for (i in 1:length(tree_tips_repaired)) {
  tree_tips_repaired_title[i] <- str_to_title(tree_tips_repaired[i])
}
#print(tree_tips_repaired_title)

# Now need to put this repaired form into the tree object
tree[["tip.label"]] <- tree_tips_repaired_title

# Identify any differences between the two trees.
#catch_outputs <- c(catch_outputs, "NEWTREE")
if (length(tree_tips_repaired_title) != length(taxa_ss)) {
  print("One or more taxa were not included in the synthetic tree - grrr.....")
  print("Remove these from the alignment and try again.")
  for (val in taxa_ss) {
    if (val %in% tree_tips_repaired) {} else {
      print(val)
      #catch_outputs <- c(catch_outputs, val)
    }
  }
}


# Tree comparison test function
for (i in 1:length(tree_tips_repaired_title)) {
  if (tree_tips_repaired_title[i] %in% species) {
    # Do nothing
  } else {
    print(tree_tips_repaired_title[i])
  }
}

for (i in 1:length(species)) {
  if (species[i] %in% tree_tips_repaired_title) {
   # Do nothing
  } else {
     print(species[i])
   }
}

# Write out the tree as a new object.  YOUR FILENAME HERE!
write.tree(tree, file = outTreeName)
