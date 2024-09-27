################################################################################
# Script: Tanglegrams & Ternary Plots
# March 2022 & revised August 2024
################################################################################




#===============================================================================
# Libraries: -------------------------------------------------------------------
#===============================================================================
library(dplyr)         # Data formatting & subsetting
library(tidyverse)     # Data processing
library(dendextend)    # Tanglegrams & computation of cophenetic correlations
library(corrplot)      # Plotting
library(ggplot2)       # Plotting
library(harrietr)      # melt_dist() to create data.frames from distance matrices
library(gridGraphics)



rm(list=ls())
myd <- getwd()




#===============================================================================
# Read Data: -------------------------------------------------------------------
#===============================================================================
load("Tax53J.RData")
load("Tax67J.RData")   # taxonomic

load("Phy53J.RData")
load("Phy67J.RData")   # phylo

load("Trait53J.RData")
load("Trait67J.RData") # trait




#===============================================================================
# Retrieve Hist, ContN & Cont All: ---------------------------------------------
#===============================================================================


## 53 Sites: -------------------------------------------------------------------
lTax53J_H <- lapply(Tax53J, function(x) {x[c(1:53), c(1:53)]})        # retrieves Historical 
lTax53J_C <- lapply(Tax53J, function(x) {x[c(54:106), c(54:106)]})    # retrieves ContN
lTax53J_CA <- lapply(Tax53J, function(x) {x[c(107:159), c(107:159)]}) # retrieves ContAll

lPhy53J_H <- lapply(Phy53J, function(x) {x[c(1:53), c(1:53)]})        # retrieves Historical 
lPhy53J_C <- lapply(Phy53J, function(x) {x[c(54:106), c(54:106)]})    # retrieves ContN
lPhy53J_CA <- lapply(Phy53J, function(x) {x[c(107:159), c(107:159)]}) # retrieves ContAll

lTrait53J_H <- lapply(Trait53J, function(x) {x[c(1:53), c(1:53)]})        # retrieves Historical 
lTrait53J_C <- lapply(Trait53J, function(x) {x[c(54:106), c(54:106)]})    # retrieves ContN
lTrait53J_CA <- lapply(Trait53J, function(x) {x[c(107:159), c(107:159)]}) # retrieves ContAll

lTax53 <- list("HNC"=lTax53J_H[[1]], "HNB"=lTax53J_H[[2]],
               "ContN"=lTax53J_C[[1]], "ContAll"=lTax53J_CA[[1]])

lPhy53 <- list("HNC"=lPhy53J_H[[1]], "HNB"=lPhy53J_H[[2]],
               "ContN"=lPhy53J_C[[1]], "ContAll"=lPhy53J_CA[[1]])

lTrait53 <- list("HNC"=lTrait53J_H[[1]], "HNB"=lTrait53J_H[[2]],
                 "ContN"=lTrait53J_C[[1]], "ContAll"=lTrait53J_CA[[1]])



## 67 Sites: -------------------------------------------------------------------
lTax67J_H <- lapply(Tax67J, function(x) {x[c(1:67), c(1:67)]})        # retrieves Historical 
lTax67J_C <- lapply(Tax67J, function(x) {x[c(68:134), c(68:134)]})    # retrieves ContN
identical(lTax67J_C[[1]], lTax67J_C[[2]])
lTax67J_CA <- lapply(Tax67J, function(x) {x[c(135:201), c(135:201)]}) # retrieves ContAll

lPhy67J_H <- lapply(Phy67J, function(x) {x[c(1:67), c(1:67)]})     
lPhy67J_C <- lapply(Phy67J, function(x) {x[c(68:134), c(68:134)]})    
lPhy67J_CA <- lapply(Phy67J, function(x) {x[c(135:201), c(135:201)]}) 

lTrait67J_H <- lapply(Trait67J, function(x) {x[c(1:67), c(1:67)]})      
lTrait67J_C <- lapply(Trait67J, function(x) {x[c(68:134), c(68:134)]})   
lTrait67J_CA <- lapply(Trait67J, function(x) {x[c(135:201), c(135:201)]}) 


lTax67 <- list("HNC"=lTax67J_H[[1]], "HNB"=lTax67J_H[[2]],
               "ContN"=lTax67J_C[[1]], "ContAll"=lTax67J_CA[[1]])

lPhy67 <- list("HNC"=lPhy67J_H[[1]], "HNB"=lPhy67J_H[[2]],
               "ContN"=lPhy67J_C[[1]], "ContAll"=lPhy67J_CA[[1]])

lTrait67 <- list("HNC"=lTrait67J_H[[1]], "HNB"=lTrait67J_H[[2]],
                 "ContN"=lTrait67J_C[[1]], "ContAll"=lTrait67J_CA[[1]])




#===============================================================================
# Dendograms: ------------------------------------------------------------------
#===============================================================================


dTax53 <- lapply(lTax53, function(x) {as.dist(x) %>% hclust(method="ward.D2") %>% as.dendrogram %>% highlight_branches_lwd})
dPhy53 <- lapply(lPhy53, function(x) {as.dist(x) %>% hclust(method="ward.D2") %>% as.dendrogram %>% highlight_branches_lwd})
dTrait53 <- lapply(lTrait53, function(x) {as.dist(x) %>% hclust(method="ward.D2") %>% as.dendrogram %>% highlight_branches_lwd})

dTax67 <- lapply(lTax67, function(x) {as.dist(x) %>% hclust(method="ward.D2") %>% as.dendrogram %>% highlight_branches_lwd})
dPhy67 <- lapply(lPhy67, function(x) {as.dist(x) %>% hclust(method="ward.D2") %>% as.dendrogram %>% highlight_branches_lwd})
dTrait67 <- lapply(lTrait67, function(x) {as.dist(x) %>% hclust(method="ward.D2") %>% as.dendrogram %>% highlight_branches_lwd})


#### HNC ####
Tax_Phy_HNC_53 <- dendlist(dTax53[[1]], dPhy53[[1]])         # taxo - phylo HNC
Tax_Trait_HNC_53 <- dendlist(dTax53[[1]], dTrait53[[1]])     # taxo - trait HNC
Phy_Trait_HNC_53 <- dendlist(dPhy53[[1]], dTrait53[[1]])     # phylo - trait HNC


Tax_Phy_HNC_67 <- dendlist(dTax67[[1]], dPhy67[[1]])         # taxo - phylo HNC
Tax_Trait_HNC_67 <- dendlist(dTax67[[1]], dTrait67[[1]])     # taxo - trait HNC
Phy_Trait_HNC_67 <- dendlist(dPhy67[[1]], dTrait67[[1]])     # phylo - trait HNC


#### HNB ####
Tax_Phy_HNB_53 <- dendlist(dTax53[[2]], dPhy53[[2]])         # taxo - phylo HNB
Tax_Trait_HNB_53 <- dendlist(dTax53[[2]], dTrait53[[2]])     # taxo - trait HNB
Phy_Trait_HNB_53 <- dendlist(dPhy53[[2]], dTrait53[[2]])     # phylo - trait HNB


Tax_Phy_HNB_67 <- dendlist(dTax67[[2]], dPhy67[[2]])         # taxo - phylo HNB
Tax_Trait_HNB_67 <- dendlist(dTax67[[2]], dTrait67[[2]])     # taxo - trait HNB
Phy_Trait_HNB_67 <- dendlist(dPhy67[[2]], dTrait67[[2]])     # phylo - trait HNB


#### ContN ####
Tax_Phy_ContN_53 <- dendlist(dTax53[[3]], dPhy53[[3]])         # taxo - phylo ContN
Tax_Trait_ContN_53 <- dendlist(dTax53[[3]], dTrait53[[3]])     # taxo - trait ContN
Phy_Trait_ContN_53 <- dendlist(dPhy53[[3]], dTrait53[[3]])     # phylo - trait ContN


Tax_Phy_ContN_67 <- dendlist(dTax67[[3]], dPhy67[[3]])         # taxo - phylo ContN
Tax_Trait_ContN_67 <- dendlist(dTax67[[3]], dTrait67[[3]])     # taxo - trait ContN
Phy_Trait_ContN_67 <- dendlist(dPhy67[[3]], dTrait67[[3]])     # phylo - trait ContN


#### ContAll ####
Tax_Phy_ContAll_53 <- dendlist(dTax53[[4]], dPhy53[[4]])         # taxo - phylo ContAll
Tax_Trait_ContAll_53 <- dendlist(dTax53[[4]], dTrait53[[4]])     # taxo - trait ContAll
Phy_Trait_ContAll_53 <- dendlist(dPhy53[[4]], dTrait53[[4]])     # phylo - trait ContAll


Tax_Phy_ContAll_67 <- dendlist(dTax67[[4]], dPhy67[[4]])         # taxo - phylo ContAll
Tax_Trait_ContAll_67 <- dendlist(dTax67[[4]], dTrait67[[4]])     # taxo - trait ContAll
Phy_Trait_ContAll_67 <- dendlist(dPhy67[[4]], dTrait67[[4]])     # phylo - trait ContAll




#===============================================================================
# Tanglegrams: -----------------------------------------------------------------
#===============================================================================

drawTanglegram <- function(x, a, b){
  y <- tanglegram(x,
                  common_subtrees_color_lines = T,
                  common_subtrees_color_branches = F,
                  highlight_distinct_edges  = F,
                  highlight_branches_lwd = F, 
                  margin_inner = 0.5, 
                  main_left=a,
                  main_right =b, 
                  cex_main_left=2,
                  cex_main_right=2)
  y <- y %>% untangle(method = "step2side")

  return(y) 
} 
# For drawing tanglegrams between pairs of dissimilarity matrices (not done for
# 67 sites because sites with 0s div align perfectly // no point in drawing 67 
# site tanglegrams)


## HNC : -----------------------------------------------------------------------
Tax_Phy_HNC_53plot <- drawTanglegram(Tax_Phy_HNC_53, a="Taxonomic", b="Phylogenetic")
Tax_Trait_HNC_53plot <- drawTanglegram(Tax_Trait_HNC_53, a="Taxonomic", b="Trait")
Phy_Trait_HNC_53plot <- drawTanglegram(Phy_Trait_HNC_53, a="Phylogenetic", b="Trait")


## HNB : -----------------------------------------------------------------------
Tax_Phy_HNB_53plot <- drawTanglegram(Tax_Phy_HNB_53, a="Taxonomic", b="Phylogenetic")
Tax_Trait_HNB_53plot <- drawTanglegram(Tax_Trait_HNB_53, a="Taxonomic", b="Trait")
Phy_Trait_HNB_53plot <- drawTanglegram(Phy_Trait_HNB_53, a="Phylogenetic", b="Trait")


## ContN : ---------------------------------------------------------------------
Tax_Phy_ContN_53plot <- drawTanglegram(Tax_Phy_ContN_53, a="Taxonomic", b="Phylogenetic")
Tax_Trait_ContN_53plot <- drawTanglegram(Tax_Trait_ContN_53, a="Taxonomic", b="Trait")
Phy_Trait_ContN_53plot <- drawTanglegram(Phy_Trait_ContN_53, a="Phylogenetic", b="Trait")


## ContAll : -------------------------------------------------------------------
Tax_Phy_ContAll_53plot <- drawTanglegram(Tax_Phy_ContAll_53, a="Taxonomic", b="Phylogenetic")
Tax_Trait_ContAll_53plot <- drawTanglegram(Tax_Trait_ContAll_53, a="Taxonomic", b="Trait")
Phy_Trait_ContAll_53plot <- drawTanglegram(Phy_Trait_ContAll_53, a="Phylogenetic", b="Trait")




#===============================================================================
# Null Models: -----------------------------------------------------------------
#===============================================================================
#cor_bakers_gamma(Tax_Phy_HNC_53)
#cor_bakers_gamma(dTax53[[1]], dPhy53[[1]])  # idem, OK
#cor_bakers_gamma(dTax53[[1]], dTax53[[1]])  # 1, OK


baker_gamma <- function(x,y,t1,t2,s){
  
  cor_bakers_gamma_results <- numeric(R)
  dend_mixed <- x # 1st dendogram in the pair
  
  for(i in 1:R) {
    dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)    # reshuffle
    cor_bakers_gamma_results[i] <- cor_bakers_gamma(x, dend_mixed)  # cor dendogram with itself but reshuffled
  }
    
  plot(density(cor_bakers_gamma_results),
       main = paste0("Baker's gamma distribution under H0 ", t1, "-", t2, " ", s),
       xlim = c(-1,1))
  abline(v = 0, lty = 2)
  
  cor1 <- cor_bakers_gamma(x, x) # dendogram with itself (always 1)
  abline(v = cor1, lty = 2)
  
  cor2 <- cor_bakers_gamma(x, y) # cor between tree pair
  abline(v = cor2, lty = 1, col = 2)
  
  round(sum(cor2 < cor_bakers_gamma_results)/ R, 4) # p-value for x-y
 
  title(sub = paste("One sided p-value =",
                    round(sum(cor2 < cor_bakers_gamma_results)/ R, 4)
  ))
  
  p <- round(sum(cor2 < cor_bakers_gamma_results)/ R, 4)
  
  output <- as.list(c(cor2,p))
  return(output)
}


Tax_Trait_53 <- mapply(function(a,b) {baker_gamma(a, b, "Tax", "Trait", "53")}, dTax53, dTrait53, SIMPLIFY = F) # added manually to tanglegram figure
Tax_Phy_53 <- mapply(function(a,b) {baker_gamma(a, b, "Tax", "Phy", "53")}, dTax53, dPhy53, SIMPLIFY = F) 
Trait_Phy_53 <- mapply(function(a,b) {baker_gamma(a, b, "Trait", "Phy", "53")}, dTrait53, dPhy53, SIMPLIFY = F) 

Tax_Trait_67 <- mapply(function(a,b) {baker_gamma(a, b, "Tax", "Trait", "67")}, dTax67, dTrait67, SIMPLIFY = F) 
Tax_Phy_67 <- mapply(function(a,b) {baker_gamma(a, b, "Tax", "Phy", "67")}, dTax67, dPhy67, SIMPLIFY = F) 
Trait_Phy_67 <- mapply(function(a,b) {baker_gamma(a, b, "Trait", "Phy", "67")}, dTrait67, dPhy67, SIMPLIFY = F) 




# End of script ################################################################

