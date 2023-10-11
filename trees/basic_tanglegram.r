# basic_tanglegram.r
#
# Compares an inferred phylogeny to a synthetic OTOL tree.
# Use make_syn_tree.r in advance to construct the synthetic tree.
# Use pseudo_branch_lengths.py in advance to give the synTree faux branch lengths.

# Dependencies
library(ape)
library(seqinr)
library(rotl)
library(stringr)
library(dendextend)

# User inputs - empirical tree, synthetic tree, and outgroup
empTree <- read.tree("psiblast_hmmer_v39_noSlash_aln_testPrank_aln.fasta.tree_root.txt")
synTree <- read.tree("psiblast_hmmer_v39_noSlash_aln_testMafft.fasta.tree_root.txt")
#outgroup <- c("Nisaea_denitrificans", "Minwuia_thermotolerans", "Phaeospirillum_fulvum", "Fodinicurvata_sediminis", "Halovulum_dunhuangense", "Paracoccus_sp.", "Sagittula_marina", "Kiloniella_litopenaei")
#outgroup <- c("Sneathiella_aquimaris", "Skermanella_sp._TT6", "Nisaea_denitrificans", "Minwuia_thermotolerans", "Polymorphum_sp.", "Phyllobacterium_calauticae", "Rhizobium_sp._Leaf306", "Magnetospirillum_fulvum", "Azospirillum_sp._TSO35-2", "Fodinicurvata_sediminis", "Oceanibium_sediminis", "Halovulum_dunhuangense", "Pararhodobacter_marinus", "Sagittula_marina", "Cognatishimia_sp._F0-27", "Aestuariivita_sp.", "Kiloniella_litopenaei")
#outgroup <- c("Aestuariivita_sp.", "Sagittula_marina", "Cognatishimia_sp._F0-27", "Pararhodobacter_marinus", "Oceanibium_sediminis", "Halovulum_dunhuangense")
#outgroup <- c("Citrus_sinensis", "Adiantum_capillus-veneris")
#outgroup <- c("Diacronema_lutheri", "Thraustotheca_clavata", "Bremia_lactucae", "Achlya_hypogyna", "Pythium_brassicum", "Albugo_candida", "Fistulifera_solaris", "Nitzschia_inconspicua", "Chaetoceros_tenuissimus", "Thalassiosira_oceanica")
outgroup <- c("Nitzschia_inconspicua", "Bremia_lactucae", "Ceratodon_purpureus", "Amborella_trichopoda")
outgroup_syn <- c("Nitzschia_inconspicua", "Bremia_lactucae")

# Root synTree and prepare for tanglegram
synTree <- unroot(synTree)
synTreeRoot <- root(synTree, outgroup_syn, resolve.root=TRUE)
for (i in 1:length(synTreeRoot$edge.length)) {
  if (synTreeRoot$edge.length[i] == 0.0) { synTreeRoot$edge.length[i] = 0.001 }
}
synTreeRootUlt <- chronos(synTreeRoot)
synTreeRootUltBin <- multi2di(synTreeRootUlt) # This arbitrary inference is a problem for tanglegrams...
synDend <- as.dendrogram(synTreeRootUltBin)
synLabels <- synTreeRootUltBin[["tip.label"]]

# Root empTree and prepare for tanglegram
empTree <- unroot(empTree)
empTreeRoot <- root(empTree, outgroup_syn, resolve.root=TRUE)
for (i in 1:length(empTreeRoot$edge.length)) {
  if (empTreeRoot$edge.length[i] == 0.0) { empTreeRoot$edge.length[i] = 0.001 }
}
empTreeRootUlt <- chronos(empTreeRoot)
empDend <- as.dendrogram(empTreeRootUlt)
empLabels <- empTreeRootUlt[["tip.label"]]
empLabelsFilter <- character()
extraEmpLabels <- character()
for (label in empLabels) {
  if (label %in% synLabels) {
    empLabelsFilter <- c(empLabelsFilter, label)
  } else {
    extraEmpLabels <- c(extraEmpLabels, label)
  }
}
synLabelsReorder <- character()
for (i in 1:length(empLabelsFilter)) {
  for (j in 1:length(empLabelsFilter)) {
    if (synLabels[j] == empLabelsFilter[i]) {
      synLabelsReorder <- c(synLabelsReorder, synLabels[j])
      #print(j)
      break
    }
  }
}

# Sanity check - if this isn't broken, then synLabelsReorder should be
# the same length as empLabelsFilter
if (length(synLabelsReorder) != length(empLabelsFilter)) {
  print("DANGER WILL ROBINSON - ERROR IN LIST RE-ORDERING!!!")
}

# Remove dendextend Namespace and reattach it.
unloadNamespace("dendextend")
suppressPackageStartupMessages(attachNamespace("dendextend"))

# Re-order synDend according to synLabelsReorder list
synDendRot <- rotate(synDend, synLabelsReorder)

# Make a preliminary tanglegram
dend12 <- untangle_step_rotate_1side(empDend, synDendRot)
tanglegram(dend12[[2]], dend12[[1]])

# Final tanglegram
FINAL <- untangle_step_rotate_1side(dend12[[2]], dend12[[1]])
tanglegram(FINAL[[1]], FINAL[[2]], lab.cex = 0.5) #Changing lab.cex will change size of labels.

# Optional entanglement analysis
tempEntangle <- entanglement(FINAL[[1]], FINAL[[2]], L=as.numeric(0.5), leaves_matching_method="labels")
round(tempEntangle, digits=5)

