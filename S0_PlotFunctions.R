
################################################################################################
# Script: Plot functions (to source)
# AFE
# March 2022
################################################################################################


# Libraries:------------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggtern)
library(ggrepel)
library(factoextra)
library(dendextend)
library(gridGraphics)
library(gridExtra)

################################################################################################
# A) Betadisper plot functions:-----------------------------------------------------------------

################################################################################################
# 1) Representation on the multidimensional space:----------------------------------------------

betadisper_plot <- function(x) {
  main <- ggplot() +
  geom_segment(data=seg.data[1:67,],aes(x=centPCoA1,xend=PCoA1,y=centPCoA2,yend=PCoA2, color="Historical"),alpha=0.3) + 
  geom_segment(data=seg.data[68:134,],aes(x=centPCoA1,xend=PCoA1,y=centPCoA2,yend=PCoA2, color="Contemporary"),alpha=0.3) + 
  geom_point(data=ce[1,1:3], aes(x=PCoA1,y=PCoA2, color="Historical"),size=4,shape=16) + 
  geom_point(data=ce[2,1:3], aes(x=PCoA1,y=PCoA2, color="Contemporary"),size=4,shape=17) + 
  geom_point(data=seg.data[1:67,], aes(x=PCoA1,y=PCoA2, color="Historical"),size=2,shape=16, alpha=0.6) +
  geom_point(data=seg.data[68:134,], aes(x=PCoA1,y=PCoA2, color="Contemporary"),size=2,shape=17, alpha=0.6) +
  labs(title="Taxonomic Beta Diversity",x="PCoA1",y="PCoA2") +
  theme_classic() + 
  ggpubr::color_palette("jco")+
  theme(legend.position="bottom")
  
  # Marginal densities along x axis:
  (xdens <- axis_canvas(pmain, axis = "x")+
      geom_density(data = seg.data, aes(x = PCoA1, fill = group),
                   alpha = 0.7, size = 0.2)+
      scale_fill_manual(values=c("#EFC000FF", "#0073C2FF"), name="Period"))
  
  # Marginal densities along y axis:
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  (ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
      geom_density(data = seg.data, aes(x = PCoA2, fill = group),
                   alpha = 0.7, size = 0.2)+
      coord_flip()+
      scale_fill_manual(values=c("#EFC000FF", "#0073C2FF"), name="Period"))
  
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(p2)
    
}







################################################################################################
# 2) Violin plots:------------------------------------------------------------------------------
violin_plot <- function(df, dc, group, title){

pviolin <- ggplot(data=df, aes(x=group,y=dc, colour=group, fill=group)) +
     geom_violin(alpha=0.5) +
     scale_fill_manual(values=c("#f04546", "#3591d1")) +
     labs(title=title, x="Period", y="Distance to centroid", color='Period', fill="Period")+
     geom_point(data = df %>% group_by(group) %>% summarise_all(mean), size = 3, color = "black")+ 
     theme_bw()+
     scale_y_continuous(limits=c(0,1))+
     theme(legend.position = "bottom")

return(pviolin)
}

# OR:
(pviolin <- ggplot(data=dfv,aes(x=Group,y=Distance_to_centroid, colour=g, fill=g)) +
    geom_violin(alpha=0.5) +
    scale_fill_manual(values=c("#f04546", "#3591d1")) +
    labs(title="", x="", y="Distance to centroid", color='Period', fill="Period")+
    stat_summary(fun.data=mean_sdl,geom="pointrange", color="black") +
    theme_classic()+
    scale_y_continuous(limits=c(0,1))+
    theme(legend.position = "bottom"))

################################################################################################
# B) Tanglegram representation:--------------------------------------------------------------------
draw_tanglegram <- function(x, a, b){ # x is a dendlist object, a is the left title to display, b is the right title to display
  tang <- tanglegram(x,
                  common_subtrees_color_lines = FALSE,
                  common_subtrees_color_branches = FALSE,
                  highlight_distinct_edges  = FALSE,
                  highlight_branches_lwd = FALSE, 
                  margin_inner = 0.5, 
                  main_left = a,
                  main_right = b, 
                  cex_main_left=2,
                  cex_main_right=2)
  return(tang) 
} # Tanglegrams between pairs of dissimilarity matrices



################################################################################################
# C) Null model:-----------------------------------------------------------------------------------
# http://talgalili.github.io/dendextend/articles/dendextend.html

# Below: Function to permute over the labels of one tree many times, 
# calculating the distribution under the null hypothesis (keeping the trees topologies constant)

permute_tree <- function(x){ # x is a tree to permute
set.seed(23235)
P <- 1000
cor_cophenetic_results <- numeric(P)
dend_mixed <- x
for(i in 1:P) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_cophenetic_results[i] <- cor_cophenetic(x, dend_mixed)
}
res <- cor_cophenetic_results
return(res)
} 


################################################################################################
# Below: distribution of cor_cophenetic under the null hypothesis (assuming fixed tree topologies). 
# Results: x vs x (after shuffling x labels) & tree x vs tree y:
  
plot_null <- function(x,y,a,b){ # x is the permuted tree, y is the non permuted tree, a is the result from permute_tree, b is the title to display

corxx <- cor_cophenetic(x,x)  # Correlation of a tree with itself
corxy <- cor_cophenetic(x,y)  # Correlation between tree x(permuted) and y
P <- 1000

null <- plot(density(a),
        main = paste0("Cophenetic correlation under H0", " ", b),
        xlim = c(-1,1))
        abline(v = 0, lty = 2)
        abline(v = corxx, lty = 2, col = 2)     # When comparing a tree with itself (max correlation of 1)
        abline(v = corxy, lty = 2, col = 4)     # Observed correlation           
        legend("topleft", legend = c("Maximum", "Observed"), fill = c(2,4))
        round(sum(corxy < a)/ P, 4)

        title(sub = paste("One sided p-value:",
                  "Maximum =",  round(sum(corxx < a)/ P, 4),
                  " ; Observed =",  round(sum(corxy < a)/ P, 4)))

return (null)
}


################################################################################################
# Below: bootstrap confidence interval.
bootstrap_int <- function(x,y){
set.seed(23801)
P <- 1000
x_labels <- labels(x)
y_labels <- labels(y)
cor_cophenetic_results <- numeric(P)

for(i in 1:P) {
  sampled_labels <- sample(x_labels, replace = TRUE)
  # members needs to be fixed since it will be later used in nleaves
  dend_mixed1 <- sample.dendrogram(x, 
                                   dend_labels=x_labels,
                                   fix_members=TRUE,fix_order=TRUE,fix_midpoint=FALSE,
                                   replace = TRUE, sampled_labels=sampled_labels
  )
  dend_mixed2 <- sample.dendrogram(y, dend_labels=y_labels,
                                   fix_members=TRUE,fix_order=TRUE,fix_midpoint=FALSE,
                                   replace = TRUE, sampled_labels=sampled_labels
  )                                    
  cor_cophenetic_results[i] <- cor_cophenetic(dend_mixed1, dend_mixed2, warn = FALSE)
}
res <- cor_cophenetic_results
return(res)
} # Returns warning (standard deviation is 0)


################################################################################################
# Below: plotting the confidence interval obtained by bootstrapping between two trees.
plot_bootstrap_int <- function(x,y,a,b){ # x is the permuted tree, y is the non permuted tree, a is the result from permute_tree, b is the title to display
CI95 <- quantile(a, probs=c(.025,.975))
cor_cophenetic(x, y)
par(mfrow = c(1,1))
plot(density(a),
     main = paste0("Cor cophenetic bootstrap distribution", " ", b),
     xlim = c(-1,1))
abline(v = CI95, lty = 2, col = 3)
abline(v = corxy, lty = 2, col = 2)
legend("topleft", legend =c("95% CI", "Observed"), fill = c(3,2))
}


################################################################################################
# D) ternary plots:-----------------------------------------------------------------------------

ternary_plot <- function(df, S, R, D, title){
lines <- data.frame(x = c(0.5, 0, 0.5), 
                    y = c(0.5, 0.5, 0), 
                    z = c(0, 0.5, 0.5), 
                    xen = c(1, 1, 1)/3, 
                    yen = c(1, 1, 1)/3, 
                    zen = c(1, 1, 1)/3)

tplot <- ggtern(data = df, aes(S, R, D, color=Period)) + 
    geom_point() + 
    geom_segment(data = lines, 
                 aes(x, y, z, 
                     xend = xen, yend = yen, zend = zen), 
                 color = 'red', size = 1)+
    ggtitle(title)+
    scale_color_manual(values=c("#F8766D","#00BFC4"), breaks=c("H", "C"), labels=c("Historical", "Contemporary"))+
    theme(legend.position = "bottom")

return(tplot)
}

