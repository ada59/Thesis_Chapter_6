# Thesis Chapter 6

## Extirpations & introductions as drivers of freshwater fish diversity change in an endemism hotspot

This repository contains the data and scripts used for the analyses presented in Chapter 6 of my thesis. 


### Folders:

-Plots: all figures, main & supplementary. 
-Tables: tables of summary stats for different analyses. 


### Scripts:
-**S1_SourcePhylogeny:** souring phylogenetic information from the Fish Tree of Life API. The main outputs are resolved trees for phylogenetic diversity change analyses.

-**S2_GammaDiversity:** plotting gamma diversity based on status group and obtaining regional-level PD & FRic values.

-**S3_AlphaDiversity:** alpha diversity quantification (phylogenetic, whereas taxonomic & trait are inherited from Chapter 7) and analysis of change (Friedman & Wilcoxon rank tests). This script also includes code for generating main and supplementary raster-map figures of alpha diversity change.

-**S4_PairwiseBetaDiversity:** calculating pairwise taxonomic, trait and phylogenetic Jaccard dissimilarities, running betadisper plus subsequent analyses of differences in spatial beta diversity among periods.

-**S5_Tanglegrams:** hierarchical clustering analysis. 


### Please note that:
-Some aesthetic changes in the figures were made using Power Point.

-Some tables were filled in manually.

-The code for the functions Jac & DJac for computing dissimilarities & FD_MLE() for computing trait diversity are available in the original sources (see Chapter 6 for references).

-Analyses associated with Chapter 7 in the thesis were completed first. Hence, scripts for the standardisation of the central Mexico surveys can be found in the repo: Thesis_Chapter7. Some objects are also inherited from this repo (e.g., fish morphological trait data matrix).
