################################################################################
# Script: Source Phylogeny
# March 2022 & revised August 2024
################################################################################




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
library(dplyr)
library(tidyverse)
library(fishtree) # Fish Tree of Life API
library(ape)
library(Rcpp)
library(readxl)
library(hablar)
library(phytools)
library(missForest)
library(rfishbase)
library(data.table)
library(data.tree)
library(ggtree)

rm(list=ls())
myd <- getwd()




#===============================================================================
# Read data:--------------------------------------------------------------------
#===============================================================================
path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"
load(paste0(path_lists, "/l67.RData"))    # main
load(paste0(path_lists, "/l53.RData"))    # excludes assemblages with no data in the contemporary period
load(paste0(path_lists, "/NDT67.RData"))  # for checks (rawdata list)




#===============================================================================
# Source Phylogeny: ------------------------------------------------------------
#===============================================================================
# https://rmacroecology.netlify.app/2018/02/12/mapping-plant-richness/
# The source above might be interesting to explain why you iterate over 100 trees (phylogenetic uncertainty?)
# But also it shows how to deal with the fact that sometimes some species are not in the phylogenetic tree.
# Attribute the average branch length values of the genus to the species not considered in the tree (TBC!!!)
# Add a tip to a tree: http://blog.phytools.org/2012/11/adding-single-tip-to-tree.html


## Example for understanding outputs of fishtree: ------------------------------
# Try example:
#fishI <- c("Astyanax bimaculatus", "Astyanax aeneus", "Poecilia reticulata", "Poecilia butleri")
#treeI <- fishtree_complete_phylogeny(fishI)
#treeI <- treeI[[1]]
#is.ultrametric(treeI) # TRUE
#plot(treeI)
#edgelabels(round(treeI$edge.length,1), col="red", cex=1, adj = c(0.95, 0.75)) 
# NOTE: In an ultrametric tree, species in the same Genus have the same last branch length.


## Sps list: -------------------------------------------------------------------
vecsps <- names(l67)
vecsps <- vecsps[!vecsps %in% c("SiteNameE", "DrainageBasinE", "Period", "Gila sp")] # 94, OK
gila_sp <- c("Gila conspersa", "Gila minacae", "Gila pulchra")
vecsps <- c(vecsps, gila_sp)
sort(vecsps) 


## Fishtree: -------------------------------------------------------------------
vecsps[vecsps=="Tampichthys dichroma"] <- "Tampichthys dichromus"
vecsps[vecsps=="Dajaus monticola"] <- "Agonostomus monticola"
vecsps[vecsps=="Amphilophus istlanus"] <- "Cichlasoma istlanum"
vecsps[vecsps=="Mayaheros beani"] <- "Cichlasoma_beani"
# C. mezquital is considered to be C. jordani (https://fishtreeoflife.org/taxonomy/family/Atherinopsidae/)
vecsps <- vecsps[!vecsps %in% c("Chirostoma mezquital")]


fishtree <- fishtree_complete_phylogeny(vecsps)  # No warning, all sps found with stochastic polytomy resolution.
fishtree_molecular <- fishtree_phylogeny(vecsps) # Warning: Requested 96 but only found 85 species


setdiff(fishtree[[1]]$tip.label, fishtree_molecular$tip.label) # Sps with no molecular data


save(fishtree, file=paste0(myd, "/PhyloObjects/fishtree.RData"))                         # MultiPhylo non-ultrametric object (100 results)
save(fishtree_molecular, file=paste0(myd, "/PhyloObjects/fishtree_molecular.RData"))     # MultiPhylo non-ultrametric object (only sps with molecular data)




#===============================================================================
# Ultrametricity: --------------------------------------------------------------
#===============================================================================
# https://www.r-bloggers.com/2021/07/three-ways-to-check-and-fix-ultrametric-phylogenies/
# A tree is ultrametric if the root-to-tip distance is equal among all tips
packageVersion("ape") 

# For a tree to be ultrametric:
# Default tolerance in ape is the square root of R’s machine epsilon, defined as
# the smallest positive floating-point number x such that 1 + x ≠ 1
sqrt(.Machine$double.eps) # 1.490116e-08

# From 4.0 uses the relative difference of the minimum and maximum root-to-tip 
# distances to check for ultrametricity

# e.g:
ltree_root_to_tip <- lapply(fishtree, function(x) dist.nodes(x)[1:Ntip(x), Ntip(x)+1])
ltree_rel_diff <- lapply(ltree_root_to_tip, function(x) {max(x)-min(x)/max(x)})
ltree_rel_diff > 1.490116e-08 
# Always true, so that is why the trees aren't ultrametric. Close, but not ultrametric.
# The relative difference is larger than the default tolerance (ape > v4)


## Force ultrametric: ---------------------------------------------------------- 
# Seems like the option force ultrametric is valid and OK. For a given tree, this 
# function can find the set of edge lengths with implied distances with minimum 
# sum-of-squared differences to the true distances 
# - in this case the patristic distances on our phylogeny.


fishtree_ultra <- lapply(fishtree, function(x) {force.ultrametric(x)})
fishtree_molecular_ultra <- force.ultrametric(fishtree_molecular) 


# NOTE##########################################################################
# force.ultrametric: should only be used to coerce a phylogeny 
# that fails "is.ultrametric" due to rounding
################################################################################

save(fishtree_ultra, file=paste0(myd,"/PhyloObjects/fishtree_ultra.RData"))                     # MultiPhylo object of ultrametric trees.
save(fishtree_molecular_ultra, file=paste0(myd,"/PhyloObjects/fishtree_molecular_ultra.RData")) # MultiPhylo object of ultrametric trees (molecular tree)




# ==============================================================================
## Average Gila sp: ------------------------------------------------------------
# ==============================================================================
fishtreeAv <- fishtree
fishtree_ultraAv <- fishtree_ultra
# Gila sp to the average branch length of the three Gilas used as proxies. 

g <- fishtreeAv[[1]]$tip.label[grep("Gila", fishtreeAv[[1]]$tip.label)]
gU <- fishtree_ultraAv[[1]]$tip.label[grep("Gila", fishtree_ultraAv[[1]]$tip.label)]
new_tips <- function(x=NULL, y=NULL, z=NULL){
  
  # Obtain the branch lengths for all species in the Genus:
  bl <- list()
  for (i in 1:length(x)){
    tree <- x[[i]]
    value <- which(tree$tip.label %in% y)
    sps <- tree$tip.label[value]
    names(tree$edge.length) <- tree$edge[, 2]
    branch_length <- tree$edge.length[as.character(value)]
    bl[[i]] <- branch_length
  }  
  
  # Create the individual tips
  sp_tips <- list()
  for (i in 1:length(bl)){
    b <- bl[[i]]
    b2 <- rep(b, 1)
    tip <-list(edge=matrix(c(2,1),1,2), tip.label=z, edge.length=b2, Nnode=1)
    sp_tips[[i]] <- tip
  }
  
  sp_tipsII <- lapply(sp_tips, function(x) {class(x) <- "phylo"; x})
  return(sp_tipsII)
}

tipsG_sp <- new_tips(x=fishtreeAv, y=g, z="Gila_sp")
tipsG_ultra_sp <- new_tips(x=fishtree_ultraAv, y=gU, z="Gila_sp")




#### add Gila sp ####

# non-ultrametric: 
ggtree(fishtreeAv[[67]], layout="fan") + geom_text(aes(label=node), size=2, hjust=-.3) + geom_tiplab(size=2)          # check which node (67 or any number!)
fishtree_gila_sp <- mapply(function(x,y) {bind.tip(x, "Gila_sp", edge.length=as.numeric(mean(y$edge.length)), where=102)}, fishtreeAv, tipsG_sp, SIMPLIFY = FALSE)
fishtree_gila_sp <- lapply(fishtree_gila_sp, function(x) {drop.tip(x, c("Gila_conspersa", "Gila_minacae", "Gila_pulchra"))})
ggtree(fishtree_gila_sp[[67]], layout="fan") + geom_text(aes(label=node), size=2, hjust=-.3) + geom_tiplab(size=2)    # check with node


# ultrametric:
fishtree_gila_spU <- mapply(function(x,y) {bind.tip(x, "Gila_sp", edge.length=as.numeric(mean(y$edge.length)), where=102)}, fishtree_ultraAv, tipsG_ultra_sp, SIMPLIFY = FALSE)
fishtree_gila_spU <- lapply(fishtree_gila_spU, function(x) {drop.tip(x, c("Gila_conspersa", "Gila_minacae", "Gila_pulchra"))})
ggtree(fishtree_gila_spU[[67]], layout="fan") + geom_text(aes(label=node), size=2, hjust=-.3) + geom_tiplab(size=2)    # check with node




#### add C mezquital ####

# non ultrametric:
fishtree_mod <- lapply(fishtree_gila_sp, function(x) {bind.tip(x, "Chirostoma_mezquital", edge.length=as.numeric(mean(x$edge.length)), where=40)})
ggtree(fishtree_mod[[67]], layout="fan") + geom_text(aes(label=node), size=2, hjust=-.3) + geom_tiplab(size=2)    # check with node


# ultrametric:
fishtree_modU <- lapply(fishtree_gila_spU, function(x) {bind.tip(x, "Chirostoma_mezquital", edge.length=as.numeric(mean(x$edge.length)), where=40)})
ggtree(fishtree_modU[[67]], layout="fan") + geom_text(aes(label=node), size=2, hjust=-.3) + geom_tiplab(size=2)    # check with node


save(fishtree_mod, file=paste0(myd,"/PhyloObjects/fishtree_mod.RData"))          # MultiPhylo object (additions of G sp & C mezquital)
save(fishtree_modU, file=paste0(myd,"/PhyloObjects/fishtree_modU.RData"))        # MultiPhylo object (additions of G sp & C mezquital)




# End of script ################################################################
