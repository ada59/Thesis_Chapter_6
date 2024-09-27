################################################################################
# Script: Beta Diversity
# March 2022 & revised August 2024
################################################################################




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
library(dplyr)
library(tidyverse)
library(vegan) 
library(betapart) 
library(cluster)
library(hablar)
library(harrietr)
library(ape)
library(FD)
library(picante)
library(data.table)
library(psych)
library(cowplot)
library(gridExtra)
library(data.table)
library(pairwiseAdonis)
library(lme4)
library(lmerTest)
#library(permute) # for ms
library(emmeans)
library(ggpubr)
library(cowplot)
library(MuMIn)


rm(list=ls())
myd <- getwd()


path_plots <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/Plots"
path_tables <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/Tables"




#===============================================================================
# Read data:--------------------------------------------------------------------
#===============================================================================
path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits"
path_phylo <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/PhyloObjects"


load(paste0(path_lists, "/l67.RData"))    # main
load(paste0(path_lists, "/l53.RData"))    # excludes assemblages with no data in the contemporary period
load(paste0(path_lists, "/NDT67.RData"))  # for checks (rawdata list)

load(paste0(path_traits, "/dist_mat1.RData"))      # main trait distance matrix (inherited from trait redundancy chapter)
load(paste0(path_phylo, "/fishtree_modU.RData"))   # trees for phylogenetic distance




#===============================================================================
# Beta Diversity Functions: ----------------------------------------------------
#===============================================================================
source(paste0(myd,"/Functions/Functions_Ricotta_et_al_2016.txt")) # For Jac
rm(DJac, PADDis)                                                  # Remove outdated versions
source(paste0(myd,"/Functions/DJac.R"))                           # For DJac (UPDATED by S Pavoine)
source(paste0(myd,"/Functions/PADDis.R"))                         # For obtaining all separate components (UPDATED by S Pavoine)




#===============================================================================
# Phylogenetic Distance Matrix: ------------------------------------------------
#===============================================================================
pd <- lapply(fishtree_modU, function(x) {cophenetic.phylo(x)})
pd_m <- do.call(cbind, pd)
pd_m <- array(pd_m, dim=c(dim(pd[[1]]), length(pd)))
pd_m <- as.matrix(apply(pd_m, c(1, 2), mean, na.rm = TRUE))     # Mean phylo distance across trees (i.e., 100 objects)
colnames(pd_m) <- colnames(pd[[1]])
rownames(pd_m) <- colnames(pd[[1]])
pd_m <- pd_m[,order(colnames(pd_m))]
pd_m <- pd_m[order(rownames(pd_m)),] 
pd_m <- (pd_m - min(pd_m))/(max(pd_m) - min(pd_m))              # Re-scale between 0 & 1
range(pd_m)




#===============================================================================
# Ensure names coincide: -------------------------------------------------------
#===============================================================================

#### trait matrix ####
identical(colnames(dist_mat1), rownames(dist_mat1)) # T
colnames(dist_mat1) <- gsub(" ", "_", colnames(dist_mat1))
rownames(dist_mat1) <- gsub(" ", "_", rownames(dist_mat1))
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]


#### phylo matrix ####
identical(colnames(pd_m), rownames(pd_m)) # T
colnames(pd_m) <- plyr::revalue(colnames(pd_m), c("Agonostomus_monticola"= "Dajaus_monticola",
                                                  "Cichlasoma_beani"= "Mayaheros_beani",
                                                  "Cichlasoma_istlanum"= "Amphilophus_istlanus",
                                                  "Tampichthys_dichromus"= "Tampichthys_dichroma"))
rownames(pd_m) <- plyr::revalue(rownames(pd_m), c("Agonostomus_monticola"= "Dajaus_monticola",
                                                  "Cichlasoma_beani"= "Mayaheros_beani",
                                                  "Cichlasoma_istlanum"= "Amphilophus_istlanus",
                                                  "Tampichthys_dichromus"= "Tampichthys_dichroma"))
pd_m <- pd_m[order(rownames(pd_m)), order(colnames(pd_m))]


#### check with 67 site subset ####
names(l67) <- gsub(" ", "_", names(l67))
l67 <- l67[,order(names(l67))]
l67 <- l67 %>% relocate(c(SiteNameE, DrainageBasinE, Period), .before=Algansea_aphanea)

identical(colnames(pd_m), colnames(as.matrix(l67[,c(4:ncol(l67))])))         # TRUE, equally ordered
identical(colnames(dist_mat1), colnames(as.matrix(l67[,c(4:ncol(l67))])))    # TRUE, equally ordered




#===============================================================================
# Format comm data for Jac: ----------------------------------------------------
#===============================================================================

# The interest is in spatial dissimilarities,
# and on how these compare between Historical & Contemporary


## Subsets: --------------------------------------------------------------------

## 1) 67 localities, even samples and 0s:---------------------------------------
l67 <- l67[!l67$Period=="ContE",]
sum(rowSums(l67[,-c(1:3)]) == 0)   # 35
table(l67$Period)

rownames(l67) <- paste0(l67$SiteNameE, "_", l67$Period)
l67_HNC <- subset(l67, l67$Period %in% c("HNC", "ContN", "ContAll"))   
l67_HNC <- within(l67_HNC, rm(SiteNameE, DrainageBasinE, Period))
sum(colSums(l67_HNC)==0) # 1

l67_HNB <- subset(l67, l67$Period %in% c("HNB", "ContN", "ContAll")) 
l67_HNB <- within(l67_HNB, rm(SiteNameE, DrainageBasinE, Period))
sum(colSums(l67_HNB)==0) # 0

groups67 <- list("HNC"=l67_HNC,
                 "HNB"=l67_HNB)
identical(colnames(dist_mat1), names(l67_HNC)) # T
identical(colnames(dist_mat1), names(l67_HNB)) # T


## 2) 53 localities, excludes no-fish locs:-------------------------------------
rm <- setdiff(unique(l67$SiteNameE), unique(l53$SiteNameE))

l67 <- l67[!l67$SiteNameE %in% rm,]

l53_HNC <- subset(l67, l67$Period %in% c("HNC", "ContN", "ContAll"))   
l53_HNC <- within(l53_HNC, rm(SiteNameE, DrainageBasinE, Period))
sum(colSums(l53_HNC)==0) # 3

l53_HNB <- subset(l67, l67$Period %in% c("HNB", "ContN", "ContAll")) 
l53_HNB <- within(l53_HNB, rm(SiteNameE, DrainageBasinE, Period))
sum(colSums(l53_HNB)==0) # 2

groups53 <- list("HNC"=l53_HNC,
                 "HNB"=l53_HNB)
identical(colnames(dist_mat1), names(l53_HNC)) # T
identical(colnames(dist_mat1), names(l53_HNB)) # T


## 3) Uneven samples, but no 0s:------------------------------------------------
sum(rowSums(l67[,-c(1:3)]) == 0)   # 7 additional ones
lunev <- l67[!rowSums(l67[,-c(1:3)]) == 0,]

#sum(rowSums(lunev[,-c(1:3)]) == 1) # 37
#lunev <- lunev[!rowSums(lunev[,-c(1:3)]) == 1,]
table(lunev$Period)

rownames(lunev) <- paste0(lunev$SiteNameE, "_", lunev$Period)
lunev_HNC <- subset(lunev, lunev$Period %in% c("HNC", "ContN", "ContAll"))   
lunev_HNC <- within(lunev_HNC, rm(SiteNameE, DrainageBasinE, Period))
lunev_HNB <- subset(lunev, lunev$Period %in% c("HNB", "ContN", "ContAll")) 
lunev_HNB <- within(lunev_HNB, rm(SiteNameE, DrainageBasinE, Period))
groups_unev <- list("HNC"=lunev_HNC,
                    "HNB"=lunev_HNB)
identical(names(lunev_HNC), names(lunev_HNB)) # T




#===============================================================================
# Jaccard dissimilarities: -----------------------------------------------------
#===============================================================================

# Seems no need to rm sps with no appearence for Jac functions

## Taxonomic: ------------------------------------------------------------------
Tax_unev <- lapply(groups_unev, function(x) {Jac(x)})
Tax_unevJ <- lapply(Tax_unev, function(x) {as.matrix(x$J)})

Tax53 <- lapply(groups53, function(x) {Jac(x)})
Tax53J <- lapply(Tax53, function(x) {as.matrix(x$J)})

Tax67 <- lapply(groups67, function(x) {Jac(x)})
Tax67J <- lapply(Tax67, function(x) {as.matrix(x$J)})


## Phylogenetic: ---------------------------------------------------------------
Phy_unev <- lapply(groups_unev, function(x) {DJac(x, as.dist(pd_m))})
Phy_unevJ <- lapply(Phy_unev, function(x) {as.matrix(x$J)})

Phy53 <- lapply(groups53, function(x) {DJac(x, as.dist(pd_m))})          # warnings here OK 
Phy53J <- lapply(Phy53, function(x) {as.matrix(x$J)})

Phy67 <- lapply(groups67, function(x) {DJac(x, as.dist(pd_m))})          # warnings here OK 
Phy67J <- lapply(Phy67, function(x) {as.matrix(x$J)})


## Trait: ----------------------------------------------------------------------
Trait_unev <- lapply(groups_unev, function(x) {DJac(x, as.dist(dist_mat1))})
Trait_unevJ <- lapply(Trait_unev, function(x) {as.matrix(x$J)})

Trait53 <- lapply(groups53, function(x) {DJac(x, as.dist(dist_mat1))})   # warnings here OK 
Trait53J <- lapply(Trait53, function(x) {as.matrix(x$J)})

Trait67 <- lapply(groups67, function(x) {DJac(x, as.dist(dist_mat1))})   # warnings here OK 
Trait67J <- lapply(Trait67, function(x) {as.matrix(x$J)})



## Test the influence of ColSums 0s:--------------------------------------------
#p <- matrix(ncol=4, nrow = 4, byrow = T, 
#            c(1, 0, 0, 0, 
#              0, 0, 1, 1, 
#              1, 0, 0, 1,
#              1, 0, 1, 1))
#colnames(p) <- c("sp1", "sp2", "sp3", "sp4")
#data <- data.frame(
#  sp1 = c(1.2, 2.3, 3.4, 4.5),
#  sp2 = c(2.1, 3.4, 1.6, 5.2),
#  sp3 = c(3.1, 2.5, 2.9, 3.9),
#  sp4 = c(1.8, 2.9, 3.1, 4.2)
#)

#data <- t(data)
#distance_matrix <- dist(data, method = "euclidean")
#distance_matrix <- as.matrix(distance_matrix)
#min_val <- min(distance_matrix)
#max_val <- max(distance_matrix)
#normalised_matrix <- (distance_matrix - min_val) / (max_val - min_val)

#DJac(p, as.dist(normalised_matrix))

#p2 <- p[, -2]
#normalised_matrix2 <- normalised_matrix[!rownames(normalised_matrix) == "sp2",!colnames(normalised_matrix) == "sp2"]
#DJac(p2, as.dist(normalised_matrix2)) # idem result to above, no issues

# NOTE #########################################################################
# Jac & DJac can be used with colSums = 0 with no issues
################################################################################




#===============================================================================
# Deal with NAs & NaNs: --------------------------------------------------------
#===============================================================================


## Taxonomic: ------------------------------------------------------------------
lapply(Tax_unevJ, function(x) {sum(is.na(x))})       # 0
lapply(Tax_unevJ, function(x) {sum(is.nan(x))})      # 0
lapply(Tax_unevJ, function(x) {sum(is.infinite(x))}) # 0

lapply(Tax53J, function(x) {sum(is.na(x))})       # 42
lapply(Tax53J, function(x) {sum(is.nan(x))})      # 42
lapply(Tax53J, function(x) {sum(is.infinite(x))}) # 0

lapply(Tax67J, function(x) {sum(is.na(x))})       # 1190 (these are between sites pairs when both or 1 have richness = 0)
lapply(Tax67J, function(x) {sum(is.nan(x))})      # 1190 (these are between sites pairs when both or 1 have richness = 0)
lapply(Tax67J, function(x) {sum(is.infinite(x))}) # 0

Tax53J <- lapply(Tax53J, function(x) {x[is.na(x)] <- 0;x}) # assign 0s to the NAs & NaNs above
Tax67J <- lapply(Tax67J, function(x) {x[is.na(x)] <- 0;x}) # assign 0s to the NAs & NaNs above

save(Tax_unevJ, file="Tax_unevJ.RData")
save(Tax53J, file = "Tax53J.RData")
save(Tax67J, file = "Tax67J.RData")


## Phylogenetic: ---------------------------------------------------------------
lapply(Phy_unevJ, function(x) {sum(is.na(x))})       # 0
lapply(Phy_unevJ, function(x) {sum(is.nan(x))})      # 0
lapply(Phy_unevJ, function(x) {sum(is.infinite(x))}) # 0

lapply(Phy53J, function(x) {sum(is.na(x))})       # 42
lapply(Phy53J, function(x) {sum(is.nan(x))})      # 42
lapply(Phy53J, function(x) {sum(is.infinite(x))}) # 2128

lapply(Phy67J, function(x) {sum(is.na(x))})       # 1190
lapply(Phy67J, function(x) {sum(is.nan(x))})      # 1190
lapply(Phy67J, function(x) {sum(is.infinite(x))}) # 11620


Phy53J <- lapply(Phy53J, function(x) {x[is.na(x)] <- 0; x})       # assign 0s 
Phy53J <- lapply(Phy53J, function(x) {x[is.infinite(x)] <- 1; x}) # these are maximum differences (species exist vs species don't exist) 

Phy67J <- lapply(Phy67J, function(x) {x[is.na(x)] <- 0; x})       # assign 0s 
Phy67J <- lapply(Phy67J, function(x) {x[is.infinite(x)] <- 1; x}) # these are maximum differences (species exist vs species don't exist) 

save(Phy_unevJ, file="Phy_unevJ.RData")
save(Phy53J, file = "Phy53J.RData")
save(Phy67J, file = "Phy67J.RData")


## Trait: ----------------------------------------------------------------------
lapply(Trait_unevJ, function(x) {sum(is.na(x))})       # 0
lapply(Trait_unevJ, function(x) {sum(is.nan(x))})      # 0
lapply(Trait_unevJ, function(x) {sum(is.infinite(x))}) # 0

lapply(Trait53J, function(x) {sum(is.na(x))})       # idem phy53
lapply(Trait53J, function(x) {sum(is.nan(x))})      # idem phy53
lapply(Trait53J, function(x) {sum(is.infinite(x))}) # idem phy53

lapply(Trait67J, function(x) {sum(is.na(x))})       # idem phy67
lapply(Trait67J, function(x) {sum(is.nan(x))})      # idem phy67
lapply(Trait67J, function(x) {sum(is.infinite(x))}) # idem phy67

Trait53J <- lapply(Trait53J, function(x) {x[is.na(x)] <- 0;x})        # assign 0s 
Trait53J <- lapply(Trait53J, function(x) {x[is.infinite(x)] <- 1;x})  # these are maximum differences (species exist vs species don't exist) 

Trait67J <- lapply(Trait67J, function(x) {x[is.na(x)] <- 0;x})        # assign 0s 
Trait67J <- lapply(Trait67J, function(x) {x[is.infinite(x)] <- 1;x})  # these are maximum differences (species exist vs species don't exist) 


save(Trait_unevJ, file="Trait_unevJ.RData")
save(Trait53J, file = "Trait53J.RData")
save(Trait67J, file = "Trait67J.RData")




#===============================================================================
# Betadisper: ------------------------------------------------------------------
#===============================================================================
g_unev <- factor(c(rep(1, 53), rep(2,46), rep(3, 53)), labels=c("Historical", "Contemporary_Native", "Contemporary_All"))
g53 <- factor(c(rep(1, 53), rep(2,53), rep(3, 53)), labels=c("Historical", "Contemporary_Native", "Contemporary_All"))
g67 <- factor(c(rep(1, 67), rep(2,67), rep(3, 67)), labels=c("Historical", "Contemporary_Native", "Contemporary_All"))

#### Uneven HNC sites ####
tax_disper_unev_hnc <- betadisper(as.dist(Tax_unevJ$HNC), g_unev)
phy_disper_unev_hnc <- betadisper(as.dist(Phy_unevJ$HNC), g_unev)  
trait_disper_unev_hnc <- betadisper(as.dist(Trait_unevJ$HNC), g_unev)

#### Uneven HNB sites ####
tax_disper_unev_hnb <- betadisper(as.dist(Tax_unevJ$HNB), g_unev)
phy_disper_unev_hnb <- betadisper(as.dist(Phy_unevJ$HNB), g_unev)  
trait_disper_unev_hnb <- betadisper(as.dist(Trait_unevJ$HNB), g_unev)

#### 53 sites ####
tax_disper53 <- lapply(Tax53J, function(x) {betadisper(as.dist(x), g53)})  
phy_disper53 <- lapply(Phy53J, function(x) {betadisper(as.dist(x), g53)})  
trait_disper53 <- lapply(Trait53J, function(x) {betadisper(as.dist(x), g53)})  

#### 67 sites ####
tax_disper67 <- lapply(Tax67J, function(x) {betadisper(as.dist(x), g67)})  
phy_disper67 <- lapply(Phy67J, function(x) {betadisper(as.dist(x), g67)})  
trait_disper67 <- lapply(Trait67J, function(x) {betadisper(as.dist(x), g67)})  


## Preliminary tests: ----------------------------------------------------------]
anova(tax_disper_unev_hnc)   # Taxonomic (all idem results if tested with permutest, OK)
anova(phy_disper_unev_hnc)   # Phylogenetic
anova(trait_disper_unev_hnc) # Trait

anova(tax_disper_unev_hnb)   # Taxonomic
anova(phy_disper_unev_hnb)   # Phylogenetic
anova(trait_disper_unev_hnb) # Trait


# NOTE #########################################################################
# ANOVAs show no differences in groups dispersions, 
# yet this does not account for the repeated measures structure of the data!
################################################################################


# NOTE #########################################################################
# Betadisper warnings: obtained for Phylogenetic models

# In betadisper: some squared distances are negative and changed to zero.
# https://rdrr.io/cran/vegan/man/betadisper.html
# Non-metric dissimilarity coefficients can produce principal coordinate axes 
# that have negative Eigenvalues. 
# These correspond to the imaginary, non-metric part of the distance between objects. 
# If negative Eigenvalues are produced, we must correct for these imaginary distances.

# The distance to its centroid of a point is: 
# z[ij]^c = sqrt(Delta^2(u[ij]^+, c[i]^+) - Delta^2(u[ij]^-, c[i]^-)),
# where Delta^2 is the squared Euclidean distance between u[ij], 
# the principal coordinate for the jth point in the ith group, and c[i], 
# the coordinate of the centroid for the ith group. 
# The super-scripted ‘+’ and ‘-’ indicate the real and imaginary parts respectively.
# This is equation (3) in Anderson (2006). 
# If the imaginary part is greater in magnitude than the real part, 
# then we would be taking the square root of a negative value, resulting in NaN, 
# and these cases are changed to zero distances (with a warning). 
# This is in line with the behaviour of Marti Anderson's PERMDISP2 programme.
################################################################################




#===============================================================================
# Statistical differences CENTROIDS: -------------------------------------------
#===============================================================================

# Double-check that the method works: ------------------------------------------
#HNB67 <- subset(l67, l67$Period %in% c("HNB"))                 
#HNB67 <- within(HNB67, rm(SiteNameE, DrainageBasinE, Period))
#HNB672 <- HNB67
#rownames(HNB672) <- gsub("HNB", "HNB2", rownames(HNB672))
#HNB67 <- as.data.frame(rbind(HNB67, HNB67))
#Jac_test <- Jac(HNB67)
#Jac_testJ <- as.matrix(Jac_test$J)
#g67_test <- factor(c(rep(1, 67), rep(2,67)), labels=c("Historical_Test_1", "Historical_Test_2"))
#Sites <- rownames(Tax67J$HNC)
#Sites <- str_split_fixed(Sites, "_", 2)[,1]
#Sites_Test <- Sites[c(1:134)]
#adonis2(Jac_testJ ~ g67_test, permutations = 1000, strata = Sites_Test) # OK, identical, as it should.


# NOTE: ########################################################################
# If there are unequal dispersions (see below) a significant result of a permanova 
# may indicate differences in groups dispersions rather than in centroids! 
# Interpret cautiously.
################################################################################



## Strata: ---------------------------------------------------------------------
# NOTE using strata is idem to using permute with sites as blocks.
Sites_unev <- rownames(Tax_unevJ$HNC)
Sites_unev <- str_split_fixed(Sites_unev, "_", 2)[,1]   # for strata (non-independent observations)

Sites53 <- rownames(Tax53J$HNC)
Sites53 <- str_split_fixed(Sites53, "_", 2)[,1]         # for strata (non-independent observations)

Sites67 <- rownames(Tax67J$HNC)
Sites67 <- str_split_fixed(Sites67, "_", 2)[,1]         # for strata (non-independent observations)


locations <- function(l_matrix=NULL, factorN = NULL, Strata = NULL){
  overall_loc_diff <- lapply(l_matrix, function(x) {adonis2(x ~ factorN, permutations = 1000, strata = Strata)})
  return(overall_loc_diff)
}     # permanova with strata

locations_no_strata <- function(l_matrix=NULL, factorN = NULL){
  overall_loc_diff <- lapply(l_matrix, function(x) {adonis2(x ~ factorN, permutations = 1000)})
  return(overall_loc_diff)
}     # permanova no strata

pairwise_locations <- function(l_matrix=NULL, factorN = NULL, Strata = NULL, range=NULL){
  
  factorN <- droplevels(factorN[range])
  Sites <- Strata[range]
  
  pair <- lapply(l_matrix, function(x) {x[range,range]})
  pairwiseadonis2 <- lapply(pair, function(x) {adonis2(x ~ factorN, permutations = 1000, strata = Sites)})
  
  return(pairwiseadonis2)
} # function for PERMANOVA between group pairs

pairwise_locations_no_strata <- function(l_matrix=NULL, factorN = NULL, range=NULL){
  
  factorN <- droplevels(factorN[range])
  
  pair <- lapply(l_matrix, function(x) {x[range,range]})
  pairwiseadonis2 <- lapply(pair, function(x) {adonis2(x ~ factorN, permutations = 1000)})
  
  return(pairwiseadonis2)
} # function for PERMANOVA between group pairs



## main (no 0s): ---------------------------------------------------------------
Tax_unev_adonis <- locations(Tax_unevJ, g_unev, Sites_unev)
Tax_unev_adonis  

Phy_unev_adonis <- locations(Phy_unevJ, g_unev, Sites_unev)
Phy_unev_adonis 

Trait_unev_adonis <- locations(Trait_unevJ, g_unev, Sites_unev)
Trait_unev_adonis 


range_H_ContN <- c(1:99)                 # between historical (HNC/HNB) & ContN
range_H_ContAll <- c(1:53, 100:152)      # between historical (HNC/HNB) & ContAll

main1 <- rbind(as.data.frame(Tax_unev_adonis), 
              as.data.frame(Phy_unev_adonis), 
              as.data.frame(Trait_unev_adonis))


#### strata ####
Tax_H_ContN <- pairwise_locations(Tax_unevJ, g_unev, Sites_unev, range_H_ContN)
Tax_H_ContN
Tax_H_ContAll <- pairwise_locations(Tax_unevJ, g_unev, Sites_unev, range_H_ContAll)
Tax_H_ContAll


Phy_H_ContN <- pairwise_locations(Phy_unevJ, g_unev, Sites_unev, range_H_ContN)
Phy_H_ContN
Phy_H_ContAll <- pairwise_locations(Phy_unevJ, g_unev, Sites_unev, range_H_ContAll)
Phy_H_ContAll


Trait_H_ContN <- pairwise_locations(Trait_unevJ, g_unev, Sites_unev, range_H_ContN)
Trait_H_ContN
Trait_H_ContAll <- pairwise_locations(Trait_unevJ, g_unev, Sites_unev, range_H_ContAll)
Trait_H_ContAll

main2_HNC <- rbind(as.data.frame(Tax_H_ContN$HNC), 
               as.data.frame(Tax_H_ContAll$HNC),
               as.data.frame(Trait_H_ContN$HNC),
               as.data.frame(Trait_H_ContAll$HNC),
               as.data.frame(Phy_H_ContN$HNC),
               as.data.frame(Phy_H_ContAll$HNC))
main2_HNC$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
main2_HNC$Period <- rep(c("ContN", "ContAll"), each=3)
main2_HNC <- main2_HNC %>% relocate(c(Dim, Period), .before=Df)

main2_HNB <- rbind(as.data.frame(Tax_H_ContN$HNB), 
               as.data.frame(Tax_H_ContAll$HNB),
               as.data.frame(Trait_H_ContN$HNB),
               as.data.frame(Trait_H_ContAll$HNB),
               as.data.frame(Phy_H_ContN$HNB),
               as.data.frame(Phy_H_ContAll$HNB))
main2_HNB$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
main2_HNB$Period <- rep(c("ContN", "ContAll"), each=3)
#p.adjust(rep(0.000999001,times= 6), method = "BH")
main2_HNB <- main2_HNB %>% relocate(c(Dim, Period), .before=Df)

write.csv(main2_HNC, file = paste0(path_tables, "/SummaryStatsPERMANOVA_STRATA_HNC.csv"), row.names = F)
write.csv(main2_HNB, file = paste0(path_tables, "/SummaryStatsPERMANOVA_STRATA_HNB.csv"), row.names = F)


#### no strata ####
Tax_H_ContN <- pairwise_locations_no_strata(Tax_unevJ, g_unev, range_H_ContN)
Tax_H_ContN
Tax_H_ContAll <- pairwise_locations_no_strata(Tax_unevJ, g_unev, range_H_ContAll)
Tax_H_ContAll


Phy_H_ContN <- pairwise_locations_no_strata(Phy_unevJ, g_unev, range_H_ContN)
Phy_H_ContN
Phy_H_ContAll <- pairwise_locations_no_strata(Phy_unevJ, g_unev, range_H_ContAll)
class(Phy_H_ContAll)


Trait_H_ContN <- pairwise_locations_no_strata(Trait_unevJ, g_unev, range_H_ContN)
Trait_H_ContN
Trait_H_ContAll <- pairwise_locations_no_strata(Trait_unevJ, g_unev, range_H_ContAll)
Trait_H_ContAll


main2_HNC <- rbind(as.data.frame(Tax_H_ContN$HNC), 
                   as.data.frame(Tax_H_ContAll$HNC),
                   as.data.frame(Trait_H_ContN$HNC),
                   as.data.frame(Trait_H_ContAll$HNC),
                   as.data.frame(Phy_H_ContN$HNC),
                   as.data.frame(Phy_H_ContAll$HNC))
main2_HNC$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
main2_HNC$Period <- rep(c("ContN", "ContAll"), each=3)
main2_HNC <- main2_HNC %>% relocate(c(Dim, Period), .before=Df)

main2_HNB <- rbind(as.data.frame(Tax_H_ContN$HNB), 
                   as.data.frame(Tax_H_ContAll$HNB),
                   as.data.frame(Trait_H_ContN$HNB),
                   as.data.frame(Trait_H_ContAll$HNB),
                   as.data.frame(Phy_H_ContN$HNB),
                   as.data.frame(Phy_H_ContAll$HNB))
main2_HNB$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
main2_HNB$Period <- rep(c("ContN", "ContAll"), each=3)
#p.adjust(rep(0.000999001,times= 6), method = "BH")
main2_HNB <- main2_HNB %>% relocate(c(Dim, Period), .before=Df)

write.csv(main2_HNC, file = paste0(path_tables, "/SummaryStatsPERMANOVA_NO_STRATA_HNC.csv"), row.names = F)
write.csv(main2_HNB, file = paste0(path_tables, "/SummaryStatsPERMANOVA_NO_STRATA_HNB.csv"), row.names = F)



## 67 sites:--------------------------------------------------------------------
Tax67adonis <- locations(Tax67J, g67, Sites67)
Tax67adonis  

Phy67adonis <- locations(Phy67J, g67, Sites67)
Phy67adonis 

Trait67adonis <- locations(Trait67J, g67, Sites67)
Trait67adonis 


range_H_ContN_67 <- c(1:134)                # between historical (HNC/HNB) & ContN
range_H_ContAll_67 <- c(1:67, 135:201)      # between historical (HNC/HNB) & ContAll

Tax67_H_ContN_67 <- pairwise_locations(Tax67J, g67, Sites67, range_H_ContN_67)
Tax67_H_ContN_67
Tax67_H_ContAll_67 <- pairwise_locations(Tax67J, g67, Sites67, range_H_ContAll_67)
Tax67_H_ContAll_67


Phy67_H_ContN_67 <- pairwise_locations(Phy67J, g67, Sites67, range_H_ContN_67)
Phy67_H_ContN_67
Phy67_H_ContAll_67 <- pairwise_locations(Phy67J, g67, Sites67, range_H_ContAll_67)
Phy67_H_ContAll_67


Trait67_H_ContN_67 <- pairwise_locations(Trait67J, g67, Sites67, range_H_ContN_67)
Trait67_H_ContN_67
Trait67_H_ContAll_67 <- pairwise_locations(Trait67J, g67, Sites67, range_H_ContAll_67)
Trait67_H_ContAll_67

HNC_67 <- rbind(as.data.frame(Tax67_H_ContN_67$HNC), 
                   as.data.frame(Tax67_H_ContAll_67$HNC),
                   as.data.frame(Trait67_H_ContN_67$HNC),
                   as.data.frame(Trait67_H_ContAll_67$HNC),
                   as.data.frame(Phy67_H_ContN_67$HNC),
                   as.data.frame(Phy67_H_ContAll_67$HNC))
HNC_67$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
HNC_67$Period <- rep(c("ContN", "ContAll"), each=3)
HNC_67 <- HNC_67 %>% relocate(c(Dim, Period), .before=Df)

HNB_67 <- rbind(as.data.frame(Tax67_H_ContN_67$HNB), 
                   as.data.frame(Tax67_H_ContAll_67$HNB),
                   as.data.frame(Trait67_H_ContN_67$HNB),
                   as.data.frame(Trait67_H_ContAll_67$HNB),
                   as.data.frame(Phy67_H_ContN_67$HNB),
                   as.data.frame(Phy67_H_ContAll_67$HNB))
HNB_67$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
HNB_67$Period <- rep(c("ContN", "ContAll"), each=3)
#p.adjust(rep(0.000999001,times= 6), method = "BH")
HNB_67 <- HNB_67 %>% relocate(c(Dim, Period), .before=Df)

write.csv(HNC_67, file = paste0(path_tables, "/SummaryStatsPERMANOVA_STRATA_HNC_67.csv"), row.names = F)
write.csv(HNB_67, file = paste0(path_tables, "/SummaryStatsPERMANOVA_STRATA_HNB_67.csv"), row.names = F)



#### no strata ####
Tax67_H_ContN_67 <- pairwise_locations_no_strata(Tax67J, g67, range_H_ContN_67)
Tax67_H_ContN_67
Tax67_H_ContAll_67 <- pairwise_locations_no_strata(Tax67J, g67, range_H_ContAll_67)
Tax67_H_ContAll_67


Phy67_H_ContN_67 <- pairwise_locations_no_strata(Phy67J, g67, range_H_ContN_67)
Phy67_H_ContN_67
Phy67_H_ContAll_67 <- pairwise_locations_no_strata(Phy67J, g67, range_H_ContAll_67)
Phy67_H_ContAll_67


Trait67_H_ContN_67 <- pairwise_locations_no_strata(Trait67J, g67, range_H_ContN_67)
Trait67_H_ContN_67
Trait67_H_ContAll_67 <- pairwise_locations_no_strata(Trait67J, g67, range_H_ContAll_67)
Trait67_H_ContAll_67


HNC_67 <- rbind(as.data.frame(Tax67_H_ContN_67$HNC), 
                as.data.frame(Tax67_H_ContAll_67$HNC),
                as.data.frame(Trait67_H_ContN_67$HNC),
                as.data.frame(Trait67_H_ContAll_67$HNC),
                as.data.frame(Phy67_H_ContN_67$HNC),
                as.data.frame(Phy67_H_ContAll_67$HNC))
HNC_67$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
HNC_67$Period <- rep(c("ContN", "ContAll"), each=3)
HNC_67 <- HNC_67 %>% relocate(c(Dim, Period), .before=Df)

HNB_67 <- rbind(as.data.frame(Tax67_H_ContN_67$HNB), 
                as.data.frame(Tax67_H_ContAll_67$HNB),
                as.data.frame(Trait67_H_ContN_67$HNB),
                as.data.frame(Trait67_H_ContAll_67$HNB),
                as.data.frame(Phy67_H_ContN_67$HNB),
                as.data.frame(Phy67_H_ContAll_67$HNB))
HNB_67$Dim <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
HNB_67$Period <- rep(c("ContN", "ContAll"), each=3)
#p.adjust(rep(0.000999001,times= 6), method = "BH")
HNB_67 <- HNB_67 %>% relocate(c(Dim, Period), .before=Df)

write.csv(HNC_67, file = paste0(path_tables, "/SummaryStatsPERMANOVA_NO_STRATA_HNC_67.csv"), row.names = F)
write.csv(HNB_67, file = paste0(path_tables, "/SummaryStatsPERMANOVA_NO_STRATA_HNB_67.csv"), row.names = F)



#===============================================================================
# Statistical differences DISPERSION: ------------------------------------------
#===============================================================================

dispersions <- function(l_disper=NULL, factorN = NULL, Strata = NULL, type=NULL){
  
  if(type=="list"){
   
    l_distances <- lapply(l_disper, function(x) {data.frame(cbind("Dist"=x$distances,
                                                                  "Period"=factorN,
                                                                  "Sites"=Strata))})
    
    l_distances <- lapply(l_distances, function(x) {x$Dist <- as.numeric(x$Dist);x
    x$Period <- as.factor(x$Period);x
    x$Sites <- as.factor(x$Sites);x})
    
    l_test <- lapply(l_distances, function(x) {lmer(Dist ~ Period + (1 | Sites), data = x)})
    l_null <- lapply(l_distances, function(x) {lmer(Dist ~ 1 + (1 | Sites), data = x)})
    
    s_test <- lapply(l_test, function(x) {summary(x)})
    s_coef <- lapply(s_test, function(x) {x$coefficients})
    s_R2 <- lapply(l_test, function(x) {r.squaredGLMM(x)})
    s_R2_random <- lapply(l_null, function(x) {r.squaredGLMM(x)})
    
    em <- lapply(l_test, function(x) {emmeans(x, ~ Period)})
    em2 <- mapply(function(x,y) {as.data.frame(cbind(as.data.frame(x), as.data.frame(y[,colnames(y) == "Pr(>|t|)"])))}, em, s_coef, SIMPLIFY=F)
    em2 <- lapply(em2, function(x) {names(x)[names(x)=='y[, colnames(y) == "Pr(>|t|)"]'] <- "P";x})
    em2 <- lapply(em2, function(x) {x %>% mutate_if(is.numeric, round, digits=3)})
    #aov_test <- lapply(l_test, function(x) {anova(x)})
    
    em2 <- mapply(function(x,y) {x$R2m <- rep(y[1,1], nrow(x));x}, em2, s_R2, SIMPLIFY=F)
    em2 <- mapply(function(x,y) {x$R2c1 <- rep(y[1,2], nrow(x));x}, em2, s_R2_random, SIMPLIFY=F)
    em2 <- mapply(function(x,y) {x$R2c2 <- rep(y[1,2], nrow(x));x}, em2, s_R2, SIMPLIFY=F)
    
    l_results <- list("MODEL"=l_test,"SUMMARY"=s_test, "EEMEANS"= em2) #"ANOVA"=aov_test
    
  }else{
    
    
    l_distances <- data.frame(cbind("Dist"=l_disper$distances,"Period"=factorN,"Sites"=Strata))
    
    l_distances$Dist <- as.numeric(l_distances$Dist)
    l_distances$Period <- as.factor(l_distances$Period)
    l_distances$Sites <- as.factor(l_distances$Sites)
    
    l_test <- lmer(Dist ~ Period + (1 | Sites), data = l_distances)
    l_null <- lmer(Dist ~ 1 + (1 | Sites), data = l_distances)
    
    s_test <- summary(l_test)
    s_coef <- s_test$coefficients
    s_R2 <- r.squaredGLMM(l_test)
    s_R2_random <- r.squaredGLMM(l_null)
    
    em <- emmeans(l_test, ~ Period)
    em2 <- as.data.frame(cbind(as.data.frame(em), as.data.frame(s_coef[,colnames(s_coef) == "Pr(>|t|)"])))
    names(em2)[names(em2)=='s_coef[,colnames(s_coef) == "Pr(>|t|)"]'] <- "P"
    em2$R2m <- rep(s_R2[1,1], nrow(em2))  # variance explained by period
    em2$R2random <- rep(s_R2_random[1,2], nrow(em2)) # variance explained by site
    em2$R2c <- rep(s_R2[1,2], nrow(em2))  # variance explained by whole model
    em2 <- em2 %>% mutate_if(is.numeric, round, digits=3)
    
    l_results <- list("MODEL"=l_test,"SUMMARY"=s_test, "EEMEANS"= em2) 
    
  }
  
  return(l_results)
}


#### Taxonomic ####
Tax_unev_hnc_disp <- dispersions(tax_disper_unev_hnc, factorN = g_unev, Strata = Sites_unev, type = "dt")
Tax_unev_hnb_disp <- dispersions(tax_disper_unev_hnb, factorN = g_unev, Strata = Sites_unev, type = "dt")
Tax67_disp <- dispersions(tax_disper67, factorN = g67, Strata = Sites67, type = "list")


#### Trait ####
Trait_unev_hnc_disp <- dispersions(trait_disper_unev_hnc, factorN = g_unev, Strata = Sites_unev, type = "dt")
Trait_unev_hnb_disp <- dispersions(trait_disper_unev_hnb, factorN = g_unev, Strata = Sites_unev, type = "dt")
Trait67_disp <- dispersions(trait_disper67, factorN = g67, Strata = Sites67, type = "list")


#### Phylogenetic ####
Phy_unev_hnc_disp <- dispersions(phy_disper_unev_hnc, factorN = g_unev, Strata = Sites_unev, type = "dt")
Phy_unev_hnb_disp <- dispersions(phy_disper_unev_hnb, factorN = g_unev, Strata = Sites_unev, type = "dt")
Phy67_disp <- dispersions(phy_disper67, factorN = g67, Strata = Sites67, type = "list")


unev <- rbind(Tax_unev_hnc_disp$EEMEANS, Tax_unev_hnb_disp$EEMEANS,
              Trait_unev_hnc_disp$EEMEANS, Trait_unev_hnb_disp$EEMEANS,
              Phy_unev_hnc_disp$EEMEANS, Phy_unev_hnb_disp$EEMEANS)
names(unev)[7] <- "P"
unev$Dimension <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
unev$Baseline <- rep(c("HC", "HB", 
                       "HC", "HB",
                       "HC", "HB"), each=3)
unev <- unev %>% relocate(c(Dimension, Baseline), .before=Period)
write.csv(unev, file = paste0(path_tables, "/SummaryStatsDispersions.csv"), row.names = F)


summary67 <- rbind(Tax67_disp$EEMEANS$HNC, Tax67_disp$EEMEANS$HNB,
                 Trait67_disp$EEMEANS$HNC, Trait67_disp$EEMEANS$HNB,
                 Phy67_disp$EEMEANS$HNC, Phy67_disp$EEMEANS$HNB)
names(summary67)[7] <- "P"
summary67$Dimension <- rep(c("Taxonomic", "Trait", "Phylogenetic"), each=6)
summary67$Baseline <- rep(c("HC", "HB", 
                       "HC", "HB",
                       "HC", "HB"), each=3)
summary67 <- summary67 %>% relocate(c(Dimension, Baseline), .before=Period)
write.csv(summary67, file = paste0(path_tables, "/SummaryStats_67Sites_Dispersions.csv"), row.names = F)



## Normality of residuals: -----------------------------------------------------
png(paste0(path_plots, "/Normality_ResidualsUnevDispersion.png"), width = 700, height = 600)
par(mfrow = c(2, 3))
hist(resid(Tax_unev_hnc_disp$MODEL, type="pearson"), main = "Taxonomic (HC)", xlab = "residuals")
hist(resid(Trait_unev_hnc_disp$MODEL, type="pearson"), main="Trait (HC)", xlab = "residuals")
hist(resid(Phy_unev_hnc_disp$MODEL, type="pearson"), main = "Phylogenetic (HC)", xlab = "residuals")

hist(resid(Tax_unev_hnb_disp$MODEL, type="pearson"), main = "Taxonomic (HB)", xlab = "residuals")
hist(resid(Trait_unev_hnb_disp$MODEL, type="pearson"), main="Trait (HB)", xlab = "residuals")
hist(resid(Phy_unev_hnb_disp$MODEL, type="pearson"), main = "Phylogenetic (HB)", xlab = "residuals")
dev.off()


png(paste0(path_plots, "/Normality_Residuals67sitesDispersion.png"), width = 700, height = 600)
par(mfrow = c(2, 3))
hist(resid(Tax67_disp$MODEL$HNC, type="pearson"), main = "Taxonomic (HC) [67 sites]", xlab = "residuals")
hist(resid(Trait67_disp$MODEL$HNC, type="pearson"), main="Trait (HC) [67 sites]", xlab = "residuals")
hist(resid(Phy67_disp$MODEL$HNC, type="pearson"), main = "Phylogenetic (HC) [67 sites]", xlab = "residuals")

hist(resid(Tax67_disp$MODEL$HNB, type="pearson"), main = "Taxonomic (HB) [67 sites]", xlab = "residuals")
hist(resid(Trait67_disp$MODEL$HNB, type="pearson"), main="Trait (HB) [67 sites]", xlab = "residuals")
hist(resid(Phy67_disp$MODEL$HNB, type="pearson"), main = "Phylogenetic (HB) [67 sites]", xlab = "residuals")
dev.off()




## Homoscedasticity: -----------------------------------------------------------
## uneven samples (no 0s, HNC)

png(paste0(path_plots, "/HomoscedasticityDispersion.png"), width = 700, height = 600)
par(mfrow = c(2, 3))
plot(resid(Tax_unev_hnc_disp$MODEL, type="pearson")~fitted(Tax_unev_hnc_disp$MODEL), main = "Taxonomic (HC)", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Trait_unev_hnc_disp$MODEL, type="pearson")~fitted(Trait_unev_hnc_disp$MODEL), main = "Trait (HC)", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Phy_unev_hnc_disp$MODEL, type="pearson")~fitted(Phy_unev_hnc_disp$MODEL), main = "Phylogenetic (HC)", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")

## uneven samples (no 0s, HNB)
plot(resid(Tax_unev_hnb_disp$MODEL, type="pearson")~fitted(Tax_unev_hnb_disp$MODEL), main = "Taxonomic (HB)", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Trait_unev_hnb_disp$MODEL, type="pearson")~fitted(Trait_unev_hnb_disp$MODEL), main = "Trait (HB)", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Phy_unev_hnb_disp$MODEL, type="pearson")~fitted(Phy_unev_hnb_disp$MODEL), main = "Phylogenetic (HB)", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
dev.off()

png(paste0(path_plots, "/Homoscedasticity_67SitesDispersion.png"), width = 700, height = 600)
par(mfrow = c(3, 2))  
plot(resid(Tax67_disp$MODEL$HNC, type="pearson")~fitted(Tax67_disp$MODEL$HNC), main = "Taxonomic (HC) [67 sites]", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Trait67_disp$MODEL$HNC, type="pearson")~fitted(Trait67_disp$MODEL$HNC), main = "Trait (HC) [67 sites]", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Phy67_disp$MODEL$HNC, type="pearson")~fitted(Phy67_disp$MODEL$HNC), main = "Phylogenetic (HC) [67 sites]", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")

## uneven samples (no 0s, HNB)
plot(resid(Tax67_disp$MODEL$HNB, type="pearson")~fitted(Tax67_disp$MODEL$HNB), main = "Taxonomic (HB) [67 sites]", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Trait67_disp$MODEL$HNB, type="pearson")~fitted(Trait67_disp$MODEL$HNB), main = "Trait (HB) [67 sites]", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
plot(resid(Phy67_disp$MODEL$HNB, type="pearson")~fitted(Phy67_disp$MODEL$HNB), main = "Phylogenetic (HB) [67 sites]", ylab = "residuals", xlab = "fitted")
abline(h=0, col="red")
dev.off()




#===============================================================================
# Betadisper Plots: ------------------------------------------------------------
#===============================================================================
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/
# https://mattsigal.github.io/eqcov_supp/betadisp-ex.html

## Palette: --------------------------------------------------------------------
c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC")

## Extract plot components: ----------------------------------------------------

extract_comp <- function(l, type=NULL){
  
  if(type=="list"){
    
    ce <- lapply(l, function(x) {data.frame(grps=rownames(x$centroids),data.frame(x$centroids))})  # centroids in space
    
    ve <- lapply(l, function(x) {data.frame(group=x$group, data.frame(x$vectors))})                # vectors in space
    
    seg_data <- mapply(function(x,y){
      data.frame(cbind(y[,1:3],x[rep(1:nrow(x),as.data.frame(table(y$group))$Freq),2:3]))
    }, ce, ve, SIMPLIFY = F)                                                                       # coordinates PCoA1 & PCoA2 per group
    
    seg_data <- lapply(seg_data , function(x) {names(x) <- c("group","PCoA1","PCoA2","centPCoA1","centPCoA2");x})
    
    grp1.hull <- lapply(seg_data, function(x){
      x[x$group=="Historical",1:3][chull(x[x$group=="Historical",2:3]),]
    })                                                                                             # convex hull historical
    
    grp2.hull <- lapply(seg_data, function(x){
      x[x$group=="Contemporary_Native",1:3][chull(x[x$group=="Contemporary_Native",2:3]),]
    }) 
    
    grp3.hull <- lapply(seg_data, function(x){
      x[x$group=="Contemporary_All",1:3][chull(x[x$group=="Contemporary_All",2:3]),]
    })   
    
    all.hull <- mapply(function(x,y,z) {rbind(x, y, z)}, grp1.hull, grp2.hull, grp3.hull, SIMPLIFY = F) 
    all.hull <- lapply(all.hull, function(x) {names(x) <- c("Period", "PCoA1", "PCoA2");x})
    
    output <- list("Centroids"=ce, "data"=seg_data, "hulls"=all.hull)
  }else{
    
    ce <- data.frame(grps=rownames(l$centroids),data.frame(l$centroids))
    ve <- data.frame(group=l$group, data.frame(l$vectors))
    seg_data <- data.frame(cbind(ve[,1:3],ce[rep(1:nrow(ce),as.data.frame(table(ve$group))$Freq),2:3]))
    names(seg_data) <- c("group","PCoA1","PCoA2","centPCoA1","centPCoA2")
    
    grp1.hull <- seg_data[seg_data$group=="Historical",1:3][chull(seg_data[seg_data$group=="Historical",2:3]),]
    grp2.hull <- seg_data[seg_data$group=="Contemporary_Native",1:3][chull(seg_data[seg_data$group=="Contemporary_Native",2:3]),]
    grp3.hull <- seg_data[seg_data$group=="Contemporary_All",1:3][chull(seg_data[seg_data$group=="Contemporary_All",2:3]),]
    
    all.hull <- rbind(grp1.hull, grp2.hull, grp3.hull)
    names(all.hull) <- c("Period", "PCoA1", "PCoA2")
    
    output <- list("Centroids"=ce, "data"=seg_data, "hulls"=all.hull)
    
  }
  
  return(output) #returns centroids and coordinates because these are enough for the final plot layout
}

#### uneven sample (no 0s) HNC ####
plot_tax_disper_unev_hnc <- extract_comp(tax_disper_unev_hnc, type="dt")
plot_phy_disper_unev_hnc <- extract_comp(phy_disper_unev_hnc, type="dt")
plot_trait_disper_unev_hnc <- extract_comp(trait_disper_unev_hnc, type="dt")

#### uneven sample (no 0s) HNB ####
plot_tax_disper_unev_hnb <- extract_comp(tax_disper_unev_hnb, type="dt")
plot_phy_disper_unev_hnb <- extract_comp(phy_disper_unev_hnb, type="dt")
plot_trait_disper_unev_hnb <- extract_comp(trait_disper_unev_hnb, type="dt")

#### 67 sites ####
plot_tax_disper67 <- extract_comp(tax_disper67, type="list")
plot_phy_disper67 <- extract_comp(phy_disper67, type="list")
plot_trait_disper67 <- extract_comp(trait_disper67, type="list")


## Plotting: -------------------------------------------------------------------
betadisper_plot <- function(ce, seg.data, hull, col.h, dimension) {
  pmain <- ggplot() +
    geom_segment(data=seg.data,aes(x=centPCoA1,xend=PCoA1,y=centPCoA2,yend=PCoA2, color=group),alpha=0.2) + 
    geom_point(data=ce, aes(x=PCoA1,y=PCoA2, color=grps),size=4.5,shape=18, alpha=0.6) + 
    geom_point(data=seg.data, aes(x=PCoA1,y=PCoA2, color=group),size=1,shape=17, alpha=0.6) +
    geom_polygon(data=hull, aes(x=PCoA1,y=PCoA2, fill=Period, color=Period), alpha=0.1)+
    labs(title=paste0(dimension, " Beta Diversity"),x="PCoA1",y="PCoA2", color="Period") +
    scale_color_manual(values=c(col.h, "#F5C710", "#CD0BBC"), 
                       breaks=c("Historical", "Contemporary_Native", "Contemporary_All"), 
                       labels=c("Historical", "Contemporary Native", "Contemporary Native + Introduced"))+
    scale_fill_manual(values=c(col.h, "#F5C710", "#CD0BBC"), 
                       breaks=c("Historical", "Contemporary_Native", "Contemporary_All"), 
                       labels=c("Historical", "Contemporary Native", "Contemporary Native + Introduced"))+
    theme_classic() + 
    theme(legend.position = "none")
  
  # Marginal densities along x axis:
  (xdens <- axis_canvas(pmain, axis = "x")+
      geom_density(data = seg.data, aes(x = PCoA1, fill = group),
                   alpha = 0.3, size = 0.2)+
      scale_fill_manual(values=c(col.h, "#F5C710", "#CD0BBC"), name="Period"))
  
  # Marginal densities along y axis:
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  (ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
      geom_density(data = seg.data, aes(x = PCoA2, fill = group),
                   alpha = 0.3, size = 0.2)+
      coord_flip()+
      scale_fill_manual(values=c(col.h, "#F5C710", "#CD0BBC"), name="Period"))
  
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(p2)
  
} 




#===============================================================================
# Save Plots: ------------------------------------------------------------------
#===============================================================================

# 4 Figures of 3 panels each (HNC unev, HNB unev, HNC 67, HNB 67)


#### HNC unev (supp) ####
(tax_hnc_unev <- betadisper_plot(plot_tax_disper_unev_hnc$Centroids,
                                 plot_tax_disper_unev_hnc$data,
                                 plot_tax_disper_unev_hnc$hulls, "#28E2E5", "Taxonomic"))
(trait_hnc_unev <- betadisper_plot(plot_trait_disper_unev_hnc$Centroids,
                                   plot_trait_disper_unev_hnc$data,
                                   plot_trait_disper_unev_hnc$hulls, "#28E2E5", "Trait")) 
(phy_hnc_unev <- betadisper_plot(plot_phy_disper_unev_hnc$Centroids,
                                 plot_phy_disper_unev_hnc$data,
                                 plot_phy_disper_unev_hnc$hulls, "#28E2E5", "Phylogenetic"))   
(spatial_beta_unev_hnc <- ggarrange(tax_hnc_unev,
                          trait_hnc_unev,
                          phy_hnc_unev,
                          labels = c("A)", "B)", "C)"),
                          align="hv",
                          font.label = list(size = 16, color = "black", face = "bold", family = NULL, position = "top"))) #legends added manually
ggsave(spatial_beta_unev_hnc, filename= paste0(path_plots, "/spatial_beta_unev_hnc.jpg"), width = 8, height = 8)

#### HNB unev ####  [MAIN TEXT]
(tax_hnb_unev <- betadisper_plot(plot_tax_disper_unev_hnb$Centroids,
                                 plot_tax_disper_unev_hnb$data,
                                 plot_tax_disper_unev_hnb$hulls, "#61D04F", "Taxonomic"))
(trait_hnb_unev <- betadisper_plot(plot_trait_disper_unev_hnb$Centroids,
                                   plot_trait_disper_unev_hnb$data,
                                   plot_trait_disper_unev_hnb$hulls, "#61D04F", "Trait")) 
(phy_hnb_unev <- betadisper_plot(plot_phy_disper_unev_hnb$Centroids,
                                 plot_phy_disper_unev_hnb$data,
                                 plot_phy_disper_unev_hnb$hulls, "#61D04F", "Phylogenetic"))   
(spatial_beta_unev_hnb <- ggarrange(tax_hnb_unev,
                                    trait_hnb_unev,
                                    phy_hnb_unev,
                                    labels = c("A)", "B)", "C)"),
                                    align="hv",
                                    font.label = list(size = 16, color = "black", face = "bold", family = NULL, position = "top"))) #legends added manually
ggsave(spatial_beta_unev_hnb, filename= paste0(path_plots, "/spatial_beta_unev_hnb.jpg"), width = 8, height = 8)

#### HNC unev (supp) #### [67 SITES]
(tax_hnc_67 <- betadisper_plot(plot_tax_disper67$Centroids$HNC,
                               plot_tax_disper67$data$HNC,
                               plot_tax_disper67$hulls$HNC, "#28E2E5", "Taxonomic"))
(trait_hnc_67 <- betadisper_plot(plot_trait_disper67$Centroids$HNC,
                                 plot_trait_disper67$data$HNC,
                                 plot_trait_disper67$hulls$HNC,"#28E2E5", "Trait")) 
(phy_hnc_67 <- betadisper_plot(plot_phy_disper67$Centroids$HNC,
                               plot_phy_disper67$data$HNC,
                               plot_phy_disper67$hulls$HNC, "#28E2E5", "Phylogenetic"))   
(spatial_beta_67_hnc <- ggarrange(tax_hnc_67,
                                    trait_hnc_67,
                                    phy_hnc_67,
                                    labels = c("A)", "B)", "C)"),
                                    align="hv",
                                    font.label = list(size = 16, color = "black", face = "bold", family = NULL, position = "top"))) #legends added manually
ggsave(spatial_beta_67_hnc, filename= paste0(path_plots, "/spatial_beta_67_hnc.jpg"), width = 8, height = 8)


#### HNB unev (supp) ####  [67 SITES]
(tax_hnb_67 <- betadisper_plot(plot_tax_disper67$Centroids$HNB,
                               plot_tax_disper67$data$HNB,
                               plot_tax_disper67$hulls$HNB,"#61D04F", "Taxonomic"))
(trait_hnb_67 <- betadisper_plot(plot_trait_disper67$Centroids$HNB,
                                 plot_trait_disper67$data$HNB,
                                 plot_trait_disper67$hulls$HNB, "#61D04F", "Trait")) 
(phy_hnb_67 <- betadisper_plot(plot_phy_disper67$Centroids$HNB,
                              plot_phy_disper67$data$HNB,
                              plot_phy_disper67$hulls$HNB, "#61D04F", "Phylogenetic"))   
(spatial_beta_67_hnb <- ggarrange(tax_hnb_67,
                                    trait_hnb_67,
                                    phy_hnb_67,
                                    labels = c("A)", "B)", "C)"),
                                    align="hv",
                                    font.label = list(size = 16, color = "black", face = "bold", family = NULL, position = "top"))) #legends added manually
ggsave(spatial_beta_67_hnb, filename= paste0(path_plots, "/spatial_beta_67_hnb.jpg"), width = 8, height = 8)


#### Legends ####
(pmain_hnc <- ggplot() +
   geom_segment(data=plot_phy_disper67$data$HNC,aes(x=centPCoA1,xend=PCoA1,y=centPCoA2,yend=PCoA2, color=group),alpha=0.2) + 
   geom_point(data=plot_phy_disper67$data$HNC, aes(x=PCoA1,y=PCoA2, color=group, fill=group),size=1,shape=17, alpha=0.6) +
   geom_polygon(data=plot_phy_disper67$hulls$HNC, aes(x=PCoA1,y=PCoA2, fill=Period, color=Period), alpha=0.1)+
   labs(fill="Period", color="Period")+
   scale_color_manual(values=c("#28E2E5", "#F5C710", "#CD0BBC"), 
                     breaks=c("Historical", "Contemporary_Native", "Contemporary_All"), 
                     labels=c("Historical", "Contemporary Native", "Contemporary Native + Introduced"))+
   scale_fill_manual(values=c("#28E2E5", "#F5C710", "#CD0BBC"), 
                    breaks=c("Historical", "Contemporary_Native", "Contemporary_All"), 
                    labels=c("Historical", "Contemporary Native", "Contemporary Native + Introduced"))+
   theme_classic() + 
   theme(legend.position = "right")) # # meaningless, just for legend HNC

(pmain_hnb <- ggplot() +
    geom_segment(data=plot_phy_disper67$data$HNB,aes(x=centPCoA1,xend=PCoA1,y=centPCoA2,yend=PCoA2, color=group),alpha=0.2) + 
    geom_point(data=plot_phy_disper67$data$HNB, aes(x=PCoA1,y=PCoA2, color=group, fill=group),size=1,shape=17, alpha=0.6) +
    geom_polygon(data=plot_phy_disper67$hulls$HNB, aes(x=PCoA1,y=PCoA2, fill=Period, color=Period), alpha=0.1)+
    labs(fill="Period", color="Period")+
    scale_color_manual(values=c("#61D04F", "#F5C710", "#CD0BBC"), 
                       breaks=c("Historical", "Contemporary_Native", "Contemporary_All"), 
                       labels=c("Historical", "Contemporary Native", "Contemporary Native + Introduced"))+
    scale_fill_manual(values=c("#61D04F", "#F5C710", "#CD0BBC"), 
                      breaks=c("Historical", "Contemporary_Native", "Contemporary_All"), 
                      labels=c("Historical", "Contemporary Native", "Contemporary Native + Introduced"))+
    theme_classic() + 
    theme(legend.position = "right")) # meaningless, just for legend HNB


ggsave(pmain_hnc, filename= paste0(path_plots, "/legend_hnc_betadisper.jpg"), width = 5, height = 5) 
ggsave(pmain_hnb, filename= paste0(path_plots, "/legend_hnb_betadisper.jpg"), width = 5, height = 5) 





# End of script ################################################################
