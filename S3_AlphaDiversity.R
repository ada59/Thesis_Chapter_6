################################################################################
# Script: Alpha diversity
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
library(data.table)
library(rstatix)


rm(list=ls())
myd <- getwd()


path_plots <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/Plots"
path_tables <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_DataChapter_2/Tables"



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

load(paste0(path_traits, "/divT0.RData"))        # assemblage taxo diversity
load(paste0(path_traits, "/divF0.RData"))        # assemblage trait diversity


load(paste0(path_phylo, "/fishtree_modU.RData")) # phylo tree for plotting




#===============================================================================
# Alpha PD:---------------------------------------------------------------------
#===============================================================================
#comm = dummy = FD::dummy$abun
#tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
#hill_phylo(comm, tree, q = 0)
#tree_com1 <- keep.tip(tree, c("sp1", "sp2", "sp5", "sp6"))
#sum(tree_com1$edge.length) # idem, OK

rownames(l67) <- paste0(l67$SiteNameE, "_", l67$Period)
l67 <- within(l67, rm(SiteNameE, DrainageBasinE, Period))

# names:
names(l67) <- paste0(substr(names(l67),1,1), "_", str_split_fixed(names(l67), " ", 2)[,2])
fishtree_modU <- lapply(fishtree_modU, function(x) {
  x$tip.label[x$tip.label=="Agonostomus_monticola"] <- "Dajaus_monticola";x
  x$tip.label[x$tip.label=="Cichlasoma_beani"] <- "Mayaheros_beani";x
  x$tip.label[x$tip.label=="Cichlasoma_istlanum"] <- "Amphilophus_istlanus";x
  x$tip.label[x$tip.label=="Tampichthys_dichromus"] <- "Tampichthys_dichroma";x
})

fishtree_modU <- lapply(fishtree_modU, function(x) {x$tip.label <- paste0(substr(x$tip.label,1,1), "_", str_split_fixed(x$tip.label, "_", 2)[,2]);x})
identical(sort(unique(names(l67))), sort(unique(fishtree_modU[[1]]$tip.label))) # TRUE

# subset trees:
lassemb <- split(l67, f=rownames(l67))
ltips <- lapply(lassemb, function(x) {x[, x == 1, drop=F]})
ltips <- lapply(ltips, function(x) {names(x)})

divP0 <- list()
for (i in 1:length(ltips)){
  chara <- ltips[[i]]
  tree_chara <- lapply(fishtree_modU, function(x) {keep.tip(x, chara)})
  tree_chara <- do.call(rbind, lapply(tree_chara, function(x) {sum(x$edge.length)}))
  tree_chara_mean <- mean(tree_chara)
  divP0[[i]] <- tree_chara_mean
} # OK, warnings OK, comms with no sps
names(divP0) <- names(ltips)


# NOTE #########################################################################
# Alpha Taxonomic & Trait inherited from Script: AssemblageLevel1 (Chapter 7)
################################################################################




#===============================================================================
# ALPHA DIV DATAFRAME: ---------------------------------------------------------
#===============================================================================

divT0 <- as.data.frame(do.call(rbind, divT0))
names(divT0) <- "Taxonomic"
divF0 <- as.data.frame(do.call(rbind, divF0))
names(divF0) <- "Trait"
divP0 <- as.data.frame(do.call(rbind, divP0))
names(divP0) <- "Phylo"

Adf <- merge(divT0, divF0, by=0, all=TRUE)           # merge taxonomic & trait by rownames (i.e., by 0)
divP0$Row.names <- rownames(divP0)
Adf <- merge(Adf, divP0, by="Row.names", all=TRUE)   # merge previous & phylogenetic 

sum(is.na(Adf)) # 0
Adf <- Adf[!Adf$Row.names %like% "ContE",]

sum(Adf$Taxonomic==0) # 35 (always 0 for each facet)
sum(Adf$Taxonomic==1) # 43 (phylo does give a branch length...)
Adf$Phylo2 <- Adf$Phylo
Adf$Phylo2[Adf$Taxonomic==1] <- 1 # set phylo div of comms with only 1 sps to 1 following T0 & F0


Adf$SiteNameE <- str_split_fixed(Adf$Row.names, "_", 2)[,1] # re-separate site & period info
Adf$Period <- str_split_fixed(Adf$Row.names, "_", 2)[,2]
Adf <- within(Adf, rm(Row.names))

Adf$Latitude <- NDT67$Latitude[match(Adf$SiteNameE, NDT67$SiteNameE)]  # add coords for plot
Adf$Longitude <- NDT67$Longitude[match(Adf$SiteNameE, NDT67$SiteNameE)]

Adf <- Adf %>% relocate(c(SiteNameE, Period, Latitude, Longitude), .before=Taxonomic)

#rm <- setdiff(sort(unique(Adf$SiteNameE)), sort(unique(l53$SiteNameE))) # skewed results
#Adf <- subset(Adf, !Adf$SiteNameE %in% rm)
#table(Adf$Period) # [Added for Friedman tests]



#===============================================================================
# Friedman Tests: --------------------------------------------------------------
#===============================================================================
#https://www.datanovia.com/en/lessons/friedman-test-in-r/


## Explore data: ---------------------------------------------------------------

Adf$Period <- factor(Adf$Period, levels=c("HNC", "HNB", "ContN", "ContAll"))

#### Taxonomic ####
Adf_taxo <- Adf %>%
  group_by(Period) %>%
  get_summary_stats(Taxonomic, type = "common") 
taxo <- ggboxplot(Adf[, names(Adf) %in% c("SiteNameE", "Period", "Taxonomic")], x = "Period", y = "Taxonomic", add = "jitter", color = "Period")+
  scale_x_discrete(labels = c("HNC" = "HC", "HNB" = "HB", "ContN" = "CN", "ContAll"= "CN+I"))+
  scale_colour_manual(values = c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC"))+
  xlab("")+
  theme(legend.position = "none")
write.csv(Adf_taxo, file=paste0(path_tables, "/Adf_taxo.csv"), row.names=FALSE)

#### Trait ####
Adf_trait <- Adf %>%
  group_by(Period) %>%
  get_summary_stats(Trait, type = "common") 
trait <- ggboxplot(Adf[, names(Adf) %in% c("SiteNameE", "Period", "Trait")], x = "Period", y = "Trait", add = "jitter", color = "Period")+
  scale_x_discrete(labels = c("HNC" = "HC", "HNB" = "HB", "ContN" = "CN", "ContAll"= "CN+I"))+
  xlab("")+
  scale_colour_manual(values = c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC"))+
  theme(legend.position = "none")
write.csv(Adf_trait, file=paste0(path_tables, "/Adf_trait.csv"), row.names=FALSE)

#### Phylogenetic ####
#Adf %>%
#  group_by(Period) %>%
#  get_summary_stats(Phylo2, type = "common") # phylo when assemblages with 1 sps have phylo = 1
#phylo <- ggboxplot(Adf[, names(Adf) %in% c("SiteNameE", "Period", "Phylo2")], x = "Period", y = "Phylo2", add = "jitter", color = "Period")+
#  scale_x_discrete(labels = c("HNC" = "HC", "HNB" = "HB", "ContN" = "CN", "ContAll"= "C"))+
#  labs(y="Phylogenetic2")+
#  scale_colour_manual(values = c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC"))+
#  theme(legend.position = "none")

Adf_phylo <- Adf %>%
  group_by(Period) %>%
  get_summary_stats(Phylo, type = "common") # phylo when assemblages with 1 have the branch length for the taxon
phylo2 <- ggboxplot(Adf[, names(Adf) %in% c("SiteNameE", "Period", "Phylo")], x = "Period", y = "Phylo", add = "jitter", color = "Period")+
  scale_x_discrete(labels = c("HNC" = "HC", "HNB" = "HB", "ContN" = "CN", "ContAll"= "CN+I"))+
  labs(y="Phylogenetic")+
  scale_colour_manual(values = c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC"))+
  theme(legend.position = "none")
write.csv(Adf_phylo, file=paste0(path_tables, "/Adf_phylo.csv"), row.names=FALSE)

supp_alpha <- ggarrange(taxo,
          trait,
          phylo2,
          common.legend = T, 
          legend="none",
          labels = c("A)", "B)", "C)"),
          align="hv",
          ncol=1,
          font.label = list(size = 16, color = "black", face = "bold", family = NULL, position = "top")) # composite taxo, trait & phylo

ggsave(supp_alpha, filename= paste0(path_plots, "/Alpha_Appendix.jpg"), width = 6, height = 10)



## Friedman Tests: -------------------------------------------------------------
# Summary tables generated manually

#### Taxonomic ####

res_Tax <- Adf %>% friedman_test(Taxonomic ~ Period|SiteNameE)     # perform test
res_Tax
eff_Tax <- Adf %>% friedman_effsize(Taxonomic ~ Period|SiteNameE, ci=T, ci.type = "norm") # compute effect size (boot=1000)
eff_Tax
pwc_Tax <- Adf %>%                                                 # pairwise comparisons
  wilcox_test(Taxonomic ~ Period, paired = TRUE, p.adjust.method = "BH")
as.data.frame(pwc_Tax)


#### Trait ####

res_Trait <- Adf %>% friedman_test(Trait ~ Period|SiteNameE)     # perform test
res_Trait
eff_Trait <- Adf %>% friedman_effsize(Trait ~ Period|SiteNameE, ci=T, ci.type = "norm") # compute effect size
eff_Trait
pwc_Trait <- Adf %>%                                             # pairwise comparisons
  wilcox_test(Trait ~ Period, paired = TRUE, p.adjust.method = "BH")
as.data.frame(pwc_Trait)


#### Phylo ####

#res_Phylo <- Adf %>% friedman_test(Phylo ~ Period |SiteNameE)    # perform test
#res_Phylo
#eff_Phylo <- Adf %>% friedman_effsize(Phylo ~ Period|SiteNameE, ci=T, ci.type = "norm") # compute effect size
#eff_Phylo
#pwc_Phylo <- Adf %>%                                             # pairwise comparisons
#  wilcox_test(Phylo ~ Period, paired = TRUE, p.adjust.method = "BH")
#as.data.frame(pwc_Phylo)

res_Phylo2 <- Adf %>% friedman_test(Phylo2 ~ Period|SiteNameE)    # perform test
res_Phylo2
eff_Phylo2 <- Adf %>% friedman_effsize(Phylo2 ~ Period|SiteNameE, ci=T, ci.type = "norm") # compute effect size
eff_Phylo2
pwc_Phylo2 <- Adf %>%                                             # pairwise comparisons
  wilcox_test(Phylo2 ~ Period, paired = TRUE, p.adjust.method = "BH")
as.data.frame(pwc_Phylo2)




#===============================================================================
# Alpha diversity (PLOT): ------------------------------------------------------
#===============================================================================

## Prep Data:-------------------------------------------------------------------

Adf$TD <- (Adf$Taxonomic - min(Adf$Taxonomic))/(max(Adf$Taxonomic) - min(Adf$Taxonomic))
Adf$PD <- (Adf$Phylo - min(Adf$Phylo))/(max(Adf$Phylo) - min(Adf$Phylo))
Adf$PD2 <- (Adf$Phylo2 - min(Adf$Phylo2))/(max(Adf$Phylo2) - min(Adf$Phylo2))
Adf$FD <- (Adf$Trait - min(Adf$Trait))/(max(Adf$Trait) - min(Adf$Trait))


## Horizontal format: ----------------------------------------------------------

#### Taxonomic ####
Adf_tax <- Adf[, c(1:4, 9)]
Adf_tax_h <- spread(Adf_tax, key="Period", value="TD")
Adf_tax_h$diff_HNC_ContN <- Adf_tax_h$HNC - Adf_tax_h$ContN
Adf_tax_h$diff_HNC_ContAll <- Adf_tax_h$HNC - Adf_tax_h$ContAll
Adf_tax_h$diff_HNB_ContN <- Adf_tax_h$HNB - Adf_tax_h$ContN
Adf_tax_h$diff_HNB_ContAll <- Adf_tax_h$HNB - Adf_tax_h$ContAll


#### Trait ####
Adf_trait <- Adf[, c(1:4, 12)]
Adf_trait_h <- spread(Adf_trait, key="Period", value="FD")
Adf_trait_h$diff_HNC_ContN <- Adf_trait_h$HNC - Adf_trait_h$ContN
Adf_trait_h$diff_HNC_ContAll <- Adf_trait_h$HNC - Adf_trait_h$ContAll
Adf_trait_h$diff_HNB_ContN <- Adf_trait_h$HNB - Adf_trait_h$ContN
Adf_trait_h$diff_HNB_ContAll <- Adf_trait_h$HNB - Adf_trait_h$ContAll


#### Phylo ####
Adf_phylo <- Adf[, c(1:4, 10)]
Adf_phylo_h <- spread(Adf_phylo, key="Period", value="PD")
Adf_phylo_h$diff_HNC_ContN <- Adf_phylo_h$HNC - Adf_phylo_h$ContN
Adf_phylo_h$diff_HNC_ContAll <- Adf_phylo_h$HNC - Adf_phylo_h$ContAll
Adf_phylo_h$diff_HNB_ContN <- Adf_phylo_h$HNB - Adf_phylo_h$ContN
Adf_phylo_h$diff_HNB_ContAll <- Adf_phylo_h$HNB - Adf_phylo_h$ContAll


#### Phylo2 ####
Adf_phylo2 <- Adf[, c(1:4, 11)]
Adf_phylo2_h <- spread(Adf_phylo2, key="Period", value="PD2")
Adf_phylo2_h$diff_HNC_ContN <- Adf_phylo2_h$HNC - Adf_phylo2_h$ContN
Adf_phylo2_h$diff_HNC_ContAll <- Adf_phylo2_h$HNC - Adf_phylo2_h$ContAll
Adf_phylo2_h$diff_HNB_ContN <- Adf_phylo2_h$HNB - Adf_phylo2_h$ContN
Adf_phylo2_h$diff_HNB_ContAll <- Adf_phylo2_h$HNB - Adf_phylo2_h$ContAll



## Rasterise: ------------------------------------------------------------------
names(Adf_tax_h)
emptras = raster(ext=extent(-107,-97, 16.5, 28),  # (xmin, xmax, ymin, ymax) central Mexico
                 res=c(0.5,0.5))                  # in degrees of x and y (55.5km x 55.5km)

crs(emptras) <- CRS("+init=epsg:4326")
p <- as(emptras@extent, "SpatialPolygons")

comparison <- "diff_HNC_ContAll"                  # [SET COMPARISON HERE]

coordinates(Adf_tax_h) <- ~Longitude + Latitude
meanrTDd <- rasterize(Adf_tax_h, emptras, comparison, fun=mean)

coordinates(Adf_trait_h) <- ~Longitude + Latitude
meanrFDd <- rasterize(Adf_trait_h, emptras, comparison, fun=mean)

coordinates(Adf_phylo_h) <- ~Longitude + Latitude
meanrPDd <- rasterize(Adf_phylo_h, emptras, comparison, fun=mean)   # phylo 

coordinates(Adf_phylo2_h) <- ~Longitude + Latitude
meanrPDd2 <- rasterize(Adf_phylo2_h, emptras, comparison, fun=mean) # phylo when branch length of sites with 1 sps is set to 1

rasdfTDd <- as.data.frame(meanrTDd, xy=TRUE)
rasdfFDd <- as.data.frame(meanrFDd, xy=TRUE)
rasdfPDd <- as.data.frame(meanrPDd, xy=TRUE)
rasdfPDd2 <- as.data.frame(meanrPDd2, xy=TRUE)

rasdfTDd <- rasdfTDd[!is.na(rasdfTDd$layer),]
names(rasdfTDd)[names(rasdfTDd)=="layer"] <- "Taxonomic"
rasdfFDd <- rasdfFDd[!is.na(rasdfFDd$layer),]
names(rasdfFDd)[names(rasdfFDd)=="layer"] <- "Trait"
rasdfPDd <- rasdfPDd[!is.na(rasdfPDd$layer),]
names(rasdfPDd)[names(rasdfPDd)=="layer"] <- "Phylogenetic"
rasdfPDd2 <- rasdfPDd2[!is.na(rasdfPDd2$layer),]
names(rasdfPDd2)[names(rasdfPDd2)=="layer"] <- "Phylogenetic2"

ras_df <- merge(rasdfTDd, rasdfFDd, by=c("x", "y"))
ras_df <- merge(ras_df, rasdfPDd, by=c("x", "y"))
ras_df <- merge(ras_df, rasdfPDd2, by=c("x", "y"))
ras_df <- gather(ras_df, key="Dimension", value="Value", -c(1:2))

ras_df$Dimension <-factor(ras_df$Dimension,
                          levels=c("Taxonomic", "Trait", 
                                   "Phylogenetic", "Phylogenetic2"))
ras_df <- ras_df[ras_df$Dimension %in% c("Taxonomic", "Trait", "Phylogenetic2"),]
ras_df$Dimension[ras_df$Dimension =="Phylogenetic2"] <- "Phylogenetic"


## Check basins: ---------------------------------------------------------------
# https://rmacroecology.netlify.app/2018/02/12/mapping-plant-richness/
world_map <- map_data('world', c('Mexico'), xlim=c(-115, -90), ylim=c(14,28))
world_mapII <- world_map[world_map$long > -107 & world_map$long < -97,]
world_mapII <- world_mapII[world_mapII$lat < 28 & world_mapII$lat > 16.5,]

#NDT67_sub <- subset(NDT67, NDT67$DrainageBasinE %in% "Balsas River")
#(g_ba<-ggplot()+coord_fixed()+xlab("")+ylab(""))
#(g_ba<-g_ba+geom_polygon(data=world_mapII, aes(x=long, y=lat, group=group), colour="black", fill="white"))
#(g_ba<-g_ba+theme_classic())
#(g_ba<-g_ba+geom_point(data= NDT67_sub, aes(x = Longitude, y = Latitude, 
#                                    color = SiteNameE), size=3))


## Maps:------------------------------------------------------------------------


(g2<-ggplot()+coord_fixed()+xlab("")+ylab(""))
(g2<-g2+geom_polygon(data=world_mapII, aes(x=long, y=lat, group=group), colour="black", fill="gray60"))
(g2<-g2+theme_classic())
(g2<-g2+geom_tile(data= ras_df, aes(x = x, y = y, 
                                      fill = Value))+
    labs(fill="α diversity change") +
    scale_fill_gradient2(low = 'darkblue', mid = 'white', high = "darkred")+
    ggtitle("")+
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 12)) +
    facet_wrap(~Dimension, nrow = 1))

ggsave(g2, file=paste0(path_plots, "/AlphaDiv_HNC_ContAll.jpg"), width=12, height=8) # save based on comparison (assigned above)




# End of script ################################################################
