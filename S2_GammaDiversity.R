################################################################################
# Script: Gamma diversity
# March 2022 & revised August 2024
################################################################################




# ==============================================================================
# Libraries:--------------------------------------------------------------------
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ape)
library(ggtree)
library(purrr)
library(stringr)
library(raster)
library(gridExtra)
library(picante)
library(hillR)
library(mFD)
library(forcats)
library(taxize)
library(cluster)
library(factoextra)
library(ggrepel)
library(patchwork)
library(ggpubr)


rm(list=ls())
myd <- getwd()


path_plots <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/Plots"



#===============================================================================
# Read Data: -------------------------------------------------------------------
#===============================================================================
path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits"
path_phylo <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/PhyloObjects"
path_chapter6 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2"

load(paste0(path_lists, "/l67.RData"))        # main
load(paste0(path_lists, "/l53.RData"))        # excludes assemblages with no data in the contemporary period
load(paste0(path_lists, "/NDT67.RData"))      # for checks (rawdata list)

load(paste0(path_traits, "/PCA.RData"))       # for plotting trait space
load(paste0(path_traits, "/tt2.RData"))       # trait data & grouping of species
load(paste0(path_traits, "/dist_mat1.RData")) # trait diversity species x species mat
load(paste0(path_traits, "/divF0.RData"))     # assemblage trait diversity

load(paste0(path_phylo, "/fishtree_modU.RData")) # phylo tree for plotting




#===============================================================================
# Gamma Phylogenetic Diversity (PLOT): -----------------------------------------
#===============================================================================
set.seed(123)
sample(1:100, 1) # 31
tree <- fishtree_modU[[31]] # for plotting

rownames(tt2) <- gsub(" ", "_", rownames(tt2))
setdiff(sort(unique(tree$tip.label)), sort(unique(rownames(tt2))))
setdiff(sort(unique(rownames(tt2))), sort(unique(tree$tip.label)))

tree$tip.label[tree$tip.label=="Agonostomus_monticola"] <- "Dajaus_monticola"
tree$tip.label[tree$tip.label=="Cichlasoma_beani"] <- "Mayaheros_beani"
tree$tip.label[tree$tip.label=="Cichlasoma_istlanum"] <- "Amphilophus_istlanus"
tree$tip.label[tree$tip.label=="Tampichthys_dichromus"] <- "Tampichthys_dichroma"

tips <- data.frame("long"=tree$tip.label, "short"=rep(NA, length(tree$tip.label)))   # to use abbreviations
tips$short <- paste0(substr(tips$long,1,1), "_", str_split_fixed(tips$long, "_", 2)[,2])
tree$tip.label <- tips$short[match(tree$tip.label, tips$long)]


tt2$Short <- paste0(substr(rownames(tt2),1,1), "_", str_split_fixed(rownames(tt2), "_", 2)[,2]) # create common abbrev col
identical(sort(unique(tree$tip.label)), sort(unique(tt2$Short)))  # TRUE



## Introduced split by source: -------------------------------------------------
#cls <- list("Remaining"=tt2$Short[tt2$SourceII=="NR"],
#            "Extirpated"=tt2$Short[tt2$SourceII=="E"],
#            "Introduced A"=tt2$Short[tt2$SourceII=="IA"],
#            "Introduced B"=tt2$Short[tt2$SourceII=="IB"])
#tree <- groupOTU(tree, cls)

#ggtree(tree, layout="fan") + geom_text(aes(label=node), size=2, hjust=-.3) + geom_tiplab(size=2)    # check nodes
#(p1 <- ggtree(tree, layout="circular", aes(color=group, linetype=group)) + 
#    geom_tiplab(size=3.8) +
#   scale_color_manual(values=c("darkgray","#0072B2","#D55E00","#E69F00"), 
#                      breaks=c("Remaining","Extirpated", "Introduced A", "Introduced B"),
#                      labels=c("Remaining","Extirpated", "Introduced A", "Introduced B")) +
#   scale_linetype_manual(values=c("dashed","solid", "solid", "solid"), 
#                          breaks=c("Remaining","Extirpated", "Introduced A", "Introduced B"),
#                          labels=c("Remaining","Extirpated", "Introduced A", "Introduced B"))+
#    labs(title="", color="Status", linetype="Status") +
#    theme(plot.margin = unit(c(14,0,14,8), "mm"),
#          legend.position = "bottom"))                      # [MAIN GAMMA TREE PLOT]
#ggsave(p1, filename=paste0(path_plots, "/Gamma_Phylo_Sources.jpg"), width = 8, height = 8)


# In order to highlight the goodeidae:
#(pGoodeidae <- p1 + geom_hilight(node = 152,
#                        fill = "#CC79A7",
#                        alpha = 0.4,
#                        extend = 0.0017))                   # [GOODEIDAE]
#(pMenidia <- p1 + geom_hilight(node = 131,
#                               fill = "#CC79A7",
#                               alpha = 0.4,
#                               extend = 0.0017))            # [SILVERSIDES]
#(pBoth <- p1 + 
#    geom_hilight(node = 152,
#    fill = "#CC79A7",
#    alpha = 0.4,
#    extend = 0.0017) +
#    geom_hilight(node = 131,
#    fill = "#CC79A7",
#    alpha = 0.4,
#    extend = 0.0017))                                       # [BOTH]



## Introduced (I): -------------------------------------------------------------
clsI <- list("Remaining"=tt2$Short[tt2$SourceIIBasic=="NR"],
             "Extirpated"=tt2$Short[tt2$SourceIIBasic=="E"],
             "Introduced"=tt2$Short[tt2$SourceIIBasic=="I"])
treeI <- groupOTU(tree, clsI)

(p2 <- ggtree(treeI, layout="circular", aes(color=group)) + 
    geom_tiplab(size=3.8) +
    scale_color_manual(values=c("darkgray", "#F55F51", "#0072B2"), 
                       breaks=c("Remaining", "Introduced","Extirpated"),
                       labels=c("Remaining", "Introduced","Extirpated")) +
    #scale_linetype_manual(values=c("solid", "solid", "solid"), 
    #                      breaks=c("Remaining", "Introduced","Extirpated"),
    #                      labels=c("Remaining", "Introduced","Extirpated"))+
    labs(title="", color="Status", linetype="Status") +
    theme(plot.margin = unit(c(14,0,14,8), "mm"),
          legend.position = "none"))





#===============================================================================
# Gamma Trait Diversity (PLOT): ------------------------------------------------
#===============================================================================
tt2$SourceIIBasic <- plyr::revalue(tt2$SourceIIBasic, c("NR"="Native Remaining",
                                                        "E"="Extirpated",
                                                        "I"="Introduced"))
(p3 <- fviz_pca_biplot(PCA,
                    title = "",
                    label = "var", # hide individual labels
                    habillage = tt2$SourceIIBasic, # color by groups
                    palette = c("#0072B2", "#F55F51", "darkgray"),
                    addEllipses = TRUE, # TRUE for concentration ellipses
                    legend.title="Status in 2005",
                    col.var = "black",
                    ellipse.type = "convex"
) + scale_shape_manual(values=c(19,19,19,19)) + 
  theme_classic() +
  theme(legend.position = "none"))
(p4 <- p3 + theme(legend.position = "bottom",
                  legend.text = element_text(size=12),
                  legend.title = element_text(size=14)))




#===============================================================================
# Gamma PLOT: ------------------------------------------------------------------
#===============================================================================
(main_gamma <- ggarrange(p2,
                         p3,
                         common.legend = T, 
                         legend="bottom",
                         labels = c("A)", "B)"),
                         align="hv",
                         font.label = list(size = 16, color = "black", face = "bold", family = NULL, position = "top")))
ggsave(main_gamma, filename=paste0(path_plots, "/Gamma_Phylo_Trait.jpg"), width = 15, height = 8)
ggsave(p2, filename= paste0(path_plots, "/Gamma_Phylo_Legend.jpg"), width = 8, height = 8) 
ggsave(p4, filename= paste0(path_plots, "/Gamma_Trait_Legend.jpg"), width = 8, height = 8) 




#===============================================================================
# Gamma Phylogenetic Diversity (METRICS): --------------------------------------
#===============================================================================
# Use the Faith PD Index (Faith 1992):------------------------------------------
# Sum of total phylogenetic diversity in a tree (total branch length)


## Prep: -----------------------------------------------------------------------
setdiff(sort(unique(fishtree_modU[[1]]$tip.label)), sort(unique(tips$long)))
fishtree_mod <- lapply(fishtree_modU, function(x) {
  x$tip.label[x$tip.label=="Agonostomus_monticola"] <- "Dajaus_monticola";x
  x$tip.label[x$tip.label=="Cichlasoma_beani"] <- "Mayaheros_beani";x
  x$tip.label[x$tip.label=="Cichlasoma_istlanum"] <- "Amphilophus_istlanus";x
  x$tip.label[x$tip.label=="Tampichthys_dichromus"] <- "Tampichthys_dichroma";x
})

fishtree_mod <- lapply(fishtree_mod, function(x) {x$tip.label <- tips$short[match(x$tip.label, tips$long)];x})
length(sort(fishtree_mod[[6]]$tip.label)) # 95



length(clsI$Remaining)       # 49
49/95*100 
length(clsI$Extirpated)      # 28
28/95*100 
length(clsI$Introduced)      # 18
18/95*100  

PD_All <- do.call(rbind, lapply(fishtree_mod, function(x) {sum(x$edge.length)}))
mean(PD_All) # 2848.026


## By Group: -------------------------------------------------------------------
tree_NR <- lapply(fishtree_mod, function(x) {keep.tip(x, c(clsI$Remaining))})
PD_NR <- do.call(rbind, lapply(tree_NR, function(x) {sum(x$edge.length)}))
mean(PD_NR)             # 1704.01
1704.01/2848.026*100    # 60 %

tree_E <- lapply(fishtree_mod, function(x) {keep.tip(x, c(clsI$Extirpated))}) 
PD_E <- do.call(rbind, lapply(tree_E, function(x) {sum(x$edge.length)}))
mean(PD_E)              # 1271.227
1271.227/2848.026*100   # 45 %

tree_I <- lapply(fishtree_mod, function(x) {keep.tip(x, c(clsI$Introduced))}) 
PD_I <- do.call(rbind, lapply(tree_I, function(x) {sum(x$edge.length)}))
mean(PD_I)              # 1009.417
1009.417/2848.026*100   # 35 %



## By Survey: ------------------------------------------------------------------
hnc <- subset(l67, l67$Period=="HNC") 
str(hnc)
hnc <- hnc[, !sapply(hnc, function(x) all(x == 0))]     # 77, OK
hnc <- hnc[,!names(hnc) %in% c("SiteNameE", "DrainageBasinE", "Period")] # 74, OK

vec_hnc <- names(hnc)
vec_hnc <- gsub(" ", "_",vec_hnc)
vec_hnc <- paste0(substr(vec_hnc,1,1), "_", str_split_fixed(vec_hnc, "_", 2)[,2])
setdiff(sort(unique(vec_hnc)), sort(unique(fishtree_mod[[1]]$tip.label))) # 0, as it should
setdiff(sort(unique(fishtree_mod[[1]]$tip.label)), sort(unique(vec_hnc))) # 0, as it should


tree_HNC <- lapply(fishtree_mod, function(x) {keep.tip(x, vec_hnc)}) 
tree_HNC  <- do.call(rbind, lapply(tree_HNC, function(x) {sum(x$edge.length)}))
mean(tree_HNC)   # 2240.087

tree_HNB <- lapply(fishtree_mod, function(x) {keep.tip(x, c(clsI$Extirpated, 
                                                            clsI$Remaining))}) 
PD_HNB  <- do.call(rbind, lapply(tree_HNB, function(x) {sum(x$edge.length)}))
mean(PD_HNB)     # 2261.82

#PD_ContN idem Remaining (NR)

tree_ContAll <- lapply(fishtree_mod, function(x) {keep.tip(x, c(clsI$Introduced,
                                                                clsI$Remaining))}) 
PD_ContAll  <- do.call(rbind, lapply(tree_ContAll, function(x) {sum(x$edge.length)}))
mean(PD_ContAll) # 2298.611




#===============================================================================
# Gamma Trait Diversity (METRICS): ---------------------------------------------
#===============================================================================
nat <- tt2$Short[tt2$SourceIIBasic=="Native Remaining"] # 49
ext <- tt2$Short[tt2$SourceIIBasic=="Extirpated"]       # 28
int <- tt2$Short[tt2$SourceIIBasic=="Introduced"]       # 18
49 + 18 + 28 # 95, OK

coms <- data.frame(matrix(ncol=95, nrow=7))
names(coms) <- c(nat, ext, int)
rownames(coms) <- c("NR", "E", "I", "HNC", "HNB", "ContALL", "ALL")

coms[1,] <- rep(c(1,0), times=c(49,46))      # remaining native (also ContN)
coms[2,] <- rep(c(0,1,0), times=c(49,28,18)) # extirpated
coms[3,] <- rep(c(0,1), times=c(77,18))      # introduced

setdiff(sort(unique(c(nat, ext))), sort(unique(vec_hnc))) # "A_diazi", "A_zacapuensis", "P_turrubarensis", OK
coms[4,] <- rep(c(1,0), times=c(77,18))            # historical conservative 
coms[4,names(coms) %in% "A_diazi"] <- 0            # historical conservative
coms[4,names(coms) %in% "A_zacapuensis"] <- 0      # historical conservative 
coms[4,names(coms) %in% "P_turrubarensis"] <- 0    # historical conservative 

coms[5,] <- rep(c(1,0), times=c(77,18))      # historical broad (nat rem + ext)
coms[6,] <- rep(c(1,0,1), times=c(49,28,18)) # contemporary  (nat rem + int)
coms[7,] <- rep(1, times=95)                 # ALL


coms <- coms[,order(names(coms))]
coords <- PCA$x
rownames(coords) <- gsub("\\.", "_", rownames(coords))
coords <- coords[order(rownames(coords)),]
identical(names(coms), rownames(coords)) # TRUE


gamma_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = coords,
  asb_sp_w         = as.matrix(coms),
  ind_vect         = c("fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = F) # compute functional richness
gamma_fd$functional_diversity_indices$fric*100 # percentages (to not confuse with progress "done" messages)

gamma_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = coords[,c(1:2)],
  asb_sp_w         = as.matrix(coms),
  ind_vect         = c("fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = F) # compute functional richness
gamma_fd$functional_diversity_indices$fric*100 # percentages for a 2D space (Dim1 & Dim2)




# End of script ################################################################
