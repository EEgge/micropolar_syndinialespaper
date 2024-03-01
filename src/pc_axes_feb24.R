library(lubridate)
library(mvabund)
library(compositions)
library(vegan)
library(magrittr)
library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)
library(metagMisc)
library(plotly)
library(ragg)
library(here)
library(ggpubr)
library(cowplot)




#### Read ASV tables prepared in script "prep_syndi_asvtable.R" ####
synditab_043_rn_num <- read_delim(here("data", "synditab_043_rn_num.txt"), delim = "\t")
synditab_043_rn_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")

synditab_3200_rn_num <- read_delim(here("data", "synditab_3200_rn_num_fromsubsamp.txt"), delim = "\t")
synditab_3200_rn_tax <- read_delim(here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")

#### Read env_table ####
env_table0 <- read_delim(here::here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y')) %>% 
  mutate(doy = yday(Dato))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

#### Read cluster grouping df ####
grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")


#### Join env table and group table (one for each size fraction) ####
env_table_wgrp_043 <- left_join(env_table0, grp_df_043, by = "env_sample")
env_table_wgrp_3200 <- left_join(env_table0, grp_df_3200, by = "env_sample")


#### Function for transforming ASV matrices to phyloseq objects #### (consider moving to own script and source instead)
phylosize <- function(asvtab, asvtax, samples) {
  asvtab_mat <- as.matrix(asvtab)
  asvtax_mat <- as.matrix(asvtax)
  rownames(asvtab_mat) <- asvtax$asv_id
  rownames(asvtax_mat) <- asvtax$asv_id
  asv_otutab <- otu_table(asvtab_mat, taxa_are_rows = TRUE)
  asv_taxtab <- tax_table(asvtax_mat)
  samples <- sample_data(samples)
  asv_ps <- phyloseq(asv_otutab, asv_taxtab, samples)
  return(asv_ps)
}


cols04 <- c("#999999", "#E69F00", "#56B4E9", "purple", "darkblue", "darkgreen")
cols3 <- c("#999999", "#E69F00", "#56B4E9", "purple", "darkblue", "darkgreen")


#### 0.4-3 ####
#### centered-log transformation 0.4-3###
asv043_ps <- phylosize(synditab_043_rn_num, synditab_043_rn_tax, env_table)
syndiclr_043 <- otu_table(microbiome::transform(asv043_ps, transform = "clr", target = "OTU"))

#### For mvabund-analyses, keep only ASVs where CLR value was above 0 in at least one sample within a give size fraction ####
syndiclr_abund_043 <- syndiclr_043
syndiclr_abund_043[syndiclr_abund_043 < 0] <- 0 # Set all values lower than the chosen limit to 0
syndiclrabunddf_043 <- phyloseq_to_df(syndiclr_abund_043, addtax =F) # Transform to data frame to be able to take rowsums
syndiclrabunddf_043$rowsums <- rowSums(syndiclr_abund_043)
syndiclrabunddf_043$asv_id <- rownames(syndiclr_043)
asv_abund_043 <- syndiclrabunddf_043 %>% filter(rowsums > 0) %>% dplyr::pull(asv_id) # Keep only asvs where clr value was above the limit in at least one sample
syndiclrdf_043 <- phyloseq_to_df(syndiclr_043, addtax =F) # Transform original clr-transformed to data frame (this function creates a column named 'OTU')
syndiclrdf_sub_043 <- syndiclrdf_043 %>% filter(OTU %in% asv_abund_043) # Keep OTUs/ASVs that passed the criterion
dim(syndiclrdf_sub_043) #All Asvs had clr>0 in at least one sample

adist_043 <- vegdist(t(syndiclrdf_sub_043[,-1]), method = "euclidean")

pca_043 <- prcomp(adist_043)
propvar1_043 <- round(summary(pca_043)$importance[2,1], 2)*100
propvar2_043 <- round(summary(pca_043)$importance[2,2], 2)*100

cluster_order_043 <- c("JanMar_Atl", "JanMar_Arc",  "May_epi", "Aug_epi_N03", "MayAuNo_meso", "1000m")

comptab0_043 <- cbind.data.frame("pc1" = pca_043$x[,1], "pc2" = pca_043$x[,2], "env_sample" = rownames(pca_043$x))
comptab_043 <- left_join(comptab0_043, env_table_wgrp_043, by = "env_sample") %>% 
  mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T)) %>% 
mutate(cluster = factor(pcaclust, levels = cluster_order_043, ordered = TRUE))

## Convex hull
comptab_043_chull <- comptab_043 %>% 
  group_by(cluster) %>% 
  slice(chull(pc1, pc2))

## Plot regular PCA plot with two first axes ##
pcaplot_043 <- ggplot(data = comptab_043, aes(x = pc1, y = pc2, color = month, shape = depthbin))+
  xlab(paste0("PC1, ", propvar1_043, "%"))+
  ylab(paste0("PC2, ", propvar2_043, "%"))+
  theme_cowplot(font_size = 12)+
  theme(plot.margin = unit(c(3,0,6,2), "pt"))+
  labs(color = "Month", shape = "Depth zone", fill = "Cluster")+
  scale_shape_discrete(name = "Depth zone", labels = c("Epipelagic", "Mesopelagic"))+
  geom_polygon(data = comptab_043_chull, alpha = 0.6, aes(x = pc1, y = pc2, fill = cluster), inherit.aes = F)+
  scale_fill_manual(values = cols04)+
  geom_point(size = 2)
pcaplot_043

#### 3-200 ####
#### centered-log transformation 0.4-3###
asv3200_ps <- phylosize(synditab_3200_rn_num, synditab_3200_rn_tax, env_table)
syndiclr_3200 <- otu_table(microbiome::transform(asv3200_ps, transform = "clr", target = "OTU"))

#### For mvabund-analyses, keep only ASVs where CLR value was above 0 in at least one sample within a give size fraction ####
syndiclr_abund_3200 <- syndiclr_3200
syndiclr_abund_3200[syndiclr_abund_3200 < 0] <- 0 # Set all values lower than the chosen limit to 0
syndiclrabunddf_3200 <- phyloseq_to_df(syndiclr_abund_3200, addtax =F) # Transform to data frame to be able to take rowsums
syndiclrabunddf_3200$rowsums <- rowSums(syndiclr_abund_3200)
syndiclrabunddf_3200$asv_id <- rownames(syndiclr_3200)
asv_abund_3200 <- syndiclrabunddf_3200 %>% filter(rowsums > 0) %>% dplyr::pull(asv_id) # Keep only asvs where clr value was above the limit in at least one sample
syndiclrdf_3200 <- phyloseq_to_df(syndiclr_3200, addtax =F) # Transform original clr-transformed to data frame (this function creates a column named 'OTU')
syndiclrdf_sub_3200 <- syndiclrdf_3200 %>% filter(OTU %in% asv_abund_3200) # Keep OTUs/ASVs that passed the criterion
dim(syndiclrdf_sub_3200)

adist_3200 <- vegdist(t(syndiclrdf_sub_3200[,-1]), method = "euclidean")

cluster_order_3200 <- c("JanMar_Atl", "JanMar_Arc",   "May_epi", "Aug_epi", "AuNo_meso", "1000m")

pca_3200 <- prcomp(adist_3200)
propvar1_3200 <- round(summary(pca_3200)$importance[2,1], 2)*100
propvar2_3200 <- round(summary(pca_3200)$importance[2,2], 2)*100

comptab0_3200 <- cbind.data.frame("pc1" = pca_3200$x[,1], "pc2" = pca_3200$x[,2], "env_sample" = rownames(pca_3200$x))
comptab_3200 <- left_join(comptab0_3200, env_table_wgrp_3200, by = "env_sample") %>% 
  mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T)) %>% 
  mutate(cluster = factor(pcaclust, levels = cluster_order_3200, ordered = TRUE))

## Convex hull
comptab_3200_chull <- comptab_3200 %>% 
  group_by(cluster) %>% 
  slice(chull(pc1, pc2))

## Plot regular PCA plot with two first axes ##
pcaplot_3200 <- ggplot(data = comptab_3200, aes(x = pc1, y = pc2, color = month, shape = depthbin, text = env_sample))+
  #geom_point(size = 2)+
  xlab(paste0("PC1, ", propvar1_3200, "%"))+
  ylab(paste0("PC2, ", propvar2_3200, "%"))+
  theme_cowplot(font_size = 12)+
  #scale_color_manual(values = cols3)+
  theme(plot.margin = unit(c(3,0,6,2), "pt"))+
  labs(color = "Month", shape = "Depth zone")+
  scale_shape_discrete(name = "Depth zone", labels = c("Epi-\npelagic", "Meso-\npelagic"))+
  geom_polygon(data = comptab_3200_chull, alpha = 0.6, aes(x = pc1, y = pc2, fill = cluster), inherit.aes = F)+
  scale_fill_manual(values = cols3)+
  geom_point(size = 2)
pcaplot_3200


  