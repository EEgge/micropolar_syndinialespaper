library(phyloseq)
library(ape)
library(ggtree)
library(stringr)
phylosize_wtree <- function(asvtab, asvtax, samples, tree) {
  asvtab_mat <- as.matrix(asvtab)
  asvtax_mat <- as.matrix(asvtax)
  rownames(asvtab_mat) <- asvtax$asv_id
  rownames(asvtax_mat) <- asvtax$asv_id
  asv_otutab <- otu_table(asvtab_mat, taxa_are_rows = TRUE)
  asv_taxtab <- tax_table(asvtax_mat)
  samples <- sample_data(samples)
  asv_ps <- phyloseq(asv_otutab, asv_taxtab, samples, tree)
  return(asv_ps)
}

#### Read asv tables and taxonomy (from script "prep_syndi_asvtable.R" ####
synditab_043_rn_num <- read_delim(here("data", "synditab_043_rn_num.txt"), delim = "\t") %>% as.data.frame(.)
synditab_043_rn_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t") %>% as.data.frame(.)
rownames(synditab_043_rn_num) <- synditab_043_rn_tax$asv_id
rownames(synditab_043_rn_tax) <- synditab_043_rn_tax$asv_id
synditab_043_prop_num <- sweep(synditab_043_rn_num, 2, colSums(synditab_043_rn_num), FUN = "/")
rownames(synditab_043_prop_num) <- synditab_043_rn_tax$asv_id


#### Read file with 200 most abundant ASVs ####
topseqs <- read_delim(here("data", "top200syndiasvs.txt"), delim = "\t")


#### Read env_table ####
env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

#### Read cluster grouping df, from script "cluster_barplot_figures_clr_new.R" ####
grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")

cluster_order_043 <- c("JanMar_Atl", "JanMar_Arc", "May_epi", "Aug_epi_N03", "MayAuNo_meso", "1000m")
cluster_order_3200 <- c("JanMar_Atl", "JanMar_Arc",   "May_epi", "Aug_epi", "AuNo_meso", "1000m")

#### Join env table and group table (one for each size fraction) ####
env_table_wgrp_043 <- left_join(env_table0, grp_df_043, by = "env_sample") %>% 
  tibble::column_to_rownames("env_sample") 
env_table_wgrp_3200 <- left_join(env_table0, grp_df_3200, by = "env_sample")  %>% 
  tibble::column_to_rownames("env_sample") 

#### Filter out 200 most abundant and get samples in correct order####
synditab_043_tax200 <- synditab_043_rn_tax[topseqs$asv_id,]
synditab_043_rn_num200 <- synditab_043_rn_num[topseqs$asv_id,rownames(env_table_wgrp_043)]
synditab_043_prop_num200 <- synditab_043_prop_num[topseqs$asv_id,rownames(env_table_wgrp_043)]


#syndiraxml_file <- system.file("RAxML_mafft_Syndiniales", "RAxML_bipartitionsBranchLabels.top200", package="treeio")
#treeio::read.raxml(syndiraxml_file)
#syndi200raxml <- treeio::read.raxml("RAxML_mafft_Syndiniales/RAxML_bipartitionsBranchLabels.top200")
#syndi200raxml <- treeio::read.tree("RAxML_mafft_Syndiniales/RAxML_bestTree.top200")
syndi200raxml <- treeio::read.tree("RAxML_mafft_Syndiniales/RAxML_bipartitions.top200")


#syndi200nw <- treeio::read.newick("syndi200_NJmafft.nwk") #if use treeio - have to read phylo-object in to phyloseq -- nono, depends on function
#syndi200nw$tip.label <- sub(".*\\_ASV", "", syndi200nw$tip.label)  
#syndi200nw$tip.label <- gsub(".*^", "ASV", syndi200nw$tip.label)  



syndi200raxml$tip.label <- str_sub(syndi200raxml$tip.label, end=-7)
syndi200raxml$tip.label <- gsub("Clade", "C", syndi200raxml$tip.label)

#syndi200raxml@phylo$tip.label <- gsub(".*^", "ASV", syndi200raxml@phylo$tip.label) 

syndi200raxml$tip.label

# 
# syndi200 <- ape::read.tree("syndi200_NJmafft.nh")
# syndi200$tip.label <- sub(".*\\_ASV", "", syndi200$tip.label)  
# syndi200$tip.label <- gsub(".*^", "ASV", syndi200$tip.label)  

ggtree(syndi200raxml)

syndi200phylo <- syndi200raxml
#phylo_top200_043 <- phylosize_wtree(synditab_043_rn_num200, synditab_043_tax200, env_table_wgrp_043, syndi200@phylo)
synditab_043_prop_num200

synditab_043_tax200$asv_id <- gsub("Clade", "C", synditab_043_tax200$asv_id)
synditab_043_tax200$asv_id
rownames(synditab_043_tax200) <- synditab_043_tax200$asv_id
rownames(synditab_043_prop_num200) <- synditab_043_tax200$asv_id
otu <- otu_table(synditab_043_prop_num200, taxa_are_rows = T)

tax <- tax_table(as.matrix(synditab_043_tax200)) 
env_table_wgrp_043 <- as.data.frame(env_table_wgrp_043)
#row.names(env_table_wgrp_0432) <- env_table_wgrp_043$env_sample
samples <- sample_data(env_table_wgrp_043)
#sample_names(samples) <- env_table_wgrp_043$env_sample

phylo_top200_043prop <- phyloseq(otu, tax, samples, syndi200phylo)


DG_II_C7 <- subset_taxa(phylo_top200_043prop, family == "Dino-Group-II-Clade-7")
plot_bar(DG_II_C7, fill = "asv_id")

cols04 <- c("#999999", "#E69F00", "#56B4E9", "purple", "darkblue", "darkgreen")
cols3 <- c("#999999", "#E69F00", "#56B4E9", "purple", "darkblue", "darkgreen")

names(cols3) <- cluster_order_3200

ragg::agg_png("DGIIC7_tree_jan24_2.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
phyloseq::plot_tree(DG_II_C7, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = 0.05, sizebase = 5)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))
dev.off()

DG_II_C6_043 <- subset_taxa(phylo_top200_043prop, family == "Dino-Group-II-Clade-6")

ragg::agg_png("DGIIC6_tree_mar24.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_II_C6_043, color="cluster", label.tips = "asv_id", size = "abundance", base.spacing = 0.05)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))
dev.off()


DG_I_C1_043 <- subset_taxa(phylo_top200_043prop, family == "Dino-Group-I-Clade-1")

ragg::agg_png("DGIC1_tree_jan24.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_I_C1, color="cluster", label.tips = "asv_id", size = "abundance", base.spacing = 0.05)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))
dev.off()

DG_I_C5 <- subset_taxa(phylo_top200_043prop, family == "Dino-Group-I-Clade-5")
ragg::agg_png("DGIC5_tree.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_I_C5, color="cluster", label.tips = "asv_id", size = "abundance", base.spacing = 0.05)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))
dev.off()

bs <- 0.07
dg2c7tre <- plot_tree(DG_II_C7, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))
dg2c6tre_04 <- plot_tree(DG_II_C6_043, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))
dg2c6tre_04

dg1c1tre_04 <- plot_tree(DG_I_C1_043, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = c(get(cluster_order_043[1]) ="#999999",
                                                  get(cluster_order_043[2]) = "#E69F00",
                                                  cluster_order_043[3] = "#56B4E9",
                                                  cluster_order_043[4] = "purple",
                                                  cluster_order_043[5] = "darkblue",
                                                  cluster_order_043[6] = "darkgreen")))

dg1c1tre_04 <- plot_tree(DG_I_C1_043, color="pcaclust", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = cols04)

dg1c1tre_04

dg2c6tre_04 <- plot_tree(DG_II_C6_043, color="pcaclust", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = cols04)

dg2c6tre_04


dg1c5tre <- plot_tree(DG_I_C5, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(100,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow"))

trefig_043 <- ggarrange(dg1c1tre, dg1c5tre, dg2c6tre, dg2c7tre, common.legend = T, labels = "AUTO")

ragg::agg_png("trefig_043_jan24.png", width = 25, height = 10, units = "cm", res = 300, scaling = 0.6)
trefig_043
dev.off()

asv_122 <- otu_table(subset_taxa(phylo_top200_043prop, asv_id == "ASV_94_DG-I-Clade-1"))

#### nano-micro fraction ####
#### Read asv tables and taxonomy (from script "prep_syndi_asvtable.R" ####
synditab_3200_rn_num <- read_delim(here("data", "synditab_3200_rn_num.txt"), delim = "\t") %>% as.data.frame(.)
synditab_3200_rn_tax <- read_delim(here("data", "synditab_3200_rn_tax.txt"), delim = "\t") %>% as.data.frame(.)

rownames(synditab_3200_rn_num) <- synditab_3200_rn_tax$asv_id
rownames(synditab_3200_rn_tax) <- synditab_3200_rn_tax$asv_id
synditab_3200_prop_num <- sweep(synditab_3200_rn_num, 2, colSums(synditab_3200_rn_num), FUN = "/")
rownames(synditab_3200_prop_num) <- synditab_3200_rn_tax$asv_id


#### Read file with 200 most abundant ASVs ####
topseqs <- read_delim(here("data", "top200syndiasvs.txt"), delim = "\t")


#### Read env_table ####
env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

#### Read cluster grouping df, from script "cluster_barplot_figures_clr_new.R" ####
#grp_df_3200 <- read_delim(here("data", "grp_3200_df.txt"), delim = "\t")
#grp_df_3200 <- read_delim(here("data", "grp_df_3200new.txt"), delim = "\t")

#### Join env table and group table (one for each size fraction) ####
env_table_wgrp_3200 <- left_join(env_table0, grp_df_3200, by = "env_sample") %>% 
  tibble::column_to_rownames("env_sample") 
env_table_wgrp_3200 <- left_join(env_table0, grp_df_3200, by = "env_sample")  %>% 
  tibble::column_to_rownames("env_sample") 

#### Filter out 200 most abundant and get samples in correct order####
synditab_3200_tax200 <- synditab_3200_rn_tax[topseqs$asv_id,]
synditab_3200_rn_num200 <- synditab_3200_rn_num[topseqs$asv_id,rownames(env_table_wgrp_3200)]
synditab_3200_prop_num200 <- synditab_3200_prop_num[topseqs$asv_id,rownames(env_table_wgrp_3200)]




#phylo_top200_3200 <- phylosize_wtree(synditab_3200_rn_num200, synditab_3200_tax200, env_table_wgrp_3200, syndi200@phylo)
synditab_3200_tax200$asv_id <- gsub("Clade", "C", synditab_3200_tax200$asv_id)
synditab_3200_tax200$asv_id
rownames(synditab_3200_tax200) <- synditab_3200_tax200$asv_id
rownames(synditab_3200_prop_num200) <- synditab_3200_tax200$asv_id

otu <- otu_table(synditab_3200_prop_num200, taxa_are_rows = T)
tax <- tax_table(as.matrix(synditab_3200_tax200))
env_table_wgrp_3200 <- as.data.frame(env_table_wgrp_3200)
samples <- sample_data(env_table_wgrp_3200)

phylo_top200_3200prop <- phyloseq(otu, tax, samples, syndi200raxml)

DG_II_C7 <- subset_taxa(phylo_top200_3200prop, family == "Dino-Group-II-Clade-7")
plot_bar(DG_II_C7, fill = "asv_id")

ragg::agg_png("DGIIC7_tree.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_II_C7, color="cluster", label.tips = "asv_id", size = "abundance", base.spacing = 0.05, sizebase = 5)
dev.off()

DG_II_C6_3200 <- subset_taxa(phylo_top200_3200prop, family == "Dino-Group-II-Clade-6")

ragg::agg_png("DGIIC6_tree.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_II_C6_3200, color="cluster", label.tips = "asv_id", size = "abundance", base.spacing = 0.05)
dev.off()


DG_I_C1_3200 <- subset_taxa(phylo_top200_3200prop, family == "Dino-Group-I-Clade-1")

ragg::agg_png("DGIC1_tree_3200_jan24.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_I_C1, color="cluster", label.tips = "asv_id", size = "abundance", base.spacing = 0.05)+
  scale_color_manual(name = "Cluster", values = c("dark_jan_B08" ="#999999",
                                                  "dark_1000m" = "#E69F00",
                                                  "dark_above_1000m" = "#56B4E9",
                                                  "light_epi_sum" = "yellow",
                                                  "dark_may" = "darkblue"))
dev.off()

DG_I_C5 <- subset_taxa(phylo_top200_3200prop, family == "Dino-Group-I-Clade-5")
ragg::agg_png("DGIC5_tree.png", width = 18, height = 10, units = "cm", res = 300, scaling = 0.7)
plot_tree(DG_I_C5, color="cluster", label.tips = "asv_id", size = "abundance",  nodelabf=nodeplotboot(100,50,4),base.spacing = 0.05)
dev.off()

bs <- 0.07
#dg2c7tre <- plot_tree(DG_II_C7, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(80,50,4, 1), size = "abundance", base.spacing = bs)
dg2c6tre_3200 <- plot_tree(DG_II_C6_3200, color="pcaclust", label.tips = "asv_id", nodelabf=nodeplotboot(80,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = cols3)
dg2c6tre_3200

dg1c1tre_3200 <- plot_tree(DG_I_C1_3200, color="pcaclust", label.tips = "asv_id", nodelabf=nodeplotboot(80,50,4, 1), size = "abundance", base.spacing = bs)+
  scale_color_manual(name = "Cluster", values = cols3)
#dg1c5tre <- plot_tree(DG_I_C5, color="cluster", label.tips = "asv_id", nodelabf=nodeplotboot(80,50,4, 1), size = "abundance", base.spacing = bs)
dg1c1tre_3200

trefig_3200 <- ggarrange(dg1c1tre, dg1c5tre, dg2c6tre, dg2c7tre, common.legend = T, labels = "AUTO")

treleg <- get

trefig_both <- ggarrange(dg1c1tre_04, dg2c6tre_04, dg1c1tre_3200, dg2c6tre_3200, common.legend = T, labels = "AUTO")
 
trefig_both <- dg1c1tre_04+theme(legend.position = "none")+
  dg2c6tre_04+theme(legend.position = "none")+
  dg1c1tre_3200+theme(legend.position = "none")+
  dg2c6tre_3200+theme(legend.text = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size=4)))+
  plot_layout(ncol = 2, guides = "collect")+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

trefig_both

ragg::agg_png("trefig_both_mar24.png", width = 45, height = 20, units = "cm", res = 300, scaling = 1)
trefig_both
dev.off()



