# root.r
# Simple tree rooting protocol

# Dependencies
library(ape)
library(seqinr)
library(rotl)

# Input parameters
inName <- "230927_TAB_ASREx2_v19_noSlash_aln.fasta.treefile"
outName2 <- "230927_TAB_ASREx2_v19_noSlash_aln.fasta_root.treefile"
#outgroup <- c("Allomyces_macrogynus", "Batrachochytrium_dendrobatidis", "Gorgonomyces_haynaldii")
outgroup <- c("Klebsormidium_nitens", "Marchantia_polymorpha", "Digitaria_exilis", "Malus_domestica", "Tarenaya_hassleriana", "Punica_granatum")

# Rooting
empTree <- read.tree(inName)
empTree <- unroot(empTree)
empTreeRoot <- root(empTree, outgroup, resolve.root=TRUE)
write.tree(empTreeRoot, file=outName2)


