# A lot of packages you may or may not need
library(phyloseq)
library(microbiome)
library(vegan)
library(metagMisc)
library(ggdendro)
library(here)
library(dplyr)
library(magrittr)
library(cowplot)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(ragg)
library(grid)
library(ggcharts)
library(scales)
library(ggtrace)
library(readr)
library(readxl)

#### Read ASV tables prepared in script "prep_syndi_asvtable.R" ####
synditab_043_rn_num <- read_delim(here("data", "synditab_043_rn_num.txt"), delim = "\t")
synditab_043_rn_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")

synditab_3200_rn_num <- read_delim(here("data", "synditab_3200_rn_num_fromsubsamp.txt"), delim = "\t")
synditab_3200_rn_tax <- read_delim(here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")

#### Create proportions ####
synditab_043_prop_num <- sweep(synditab_043_rn_num, 2, colSums(synditab_043_rn_num), FUN = "/")
synditab_043_prop <- cbind.data.frame(synditab_043_prop_num, synditab_043_rn_tax)

synditab_3200_prop_num <- sweep(synditab_3200_rn_num, 2, colSums(synditab_3200_rn_num), FUN = "/")
synditab_3200_prop <- cbind.data.frame(synditab_3200_prop_num, synditab_3200_rn_tax)

grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")

synditab_043_prop_long_wclust <- pivot_longer(synditab_043_prop, cols = grp_df_043$env_sample, names_to = "env_sample", values_to = "prop") %>% 
  left_join(., grp_df_043, by = "env_sample")

synditab_prop_clades_clust <- synditab_043_prop_long_wclust %>% 
  group_by(family, pcaclust, env_sample) %>% 
  reframe(propfam = sum(prop)) %>%
  group_by(family, pcaclust) %>% 
  reframe(meanprop = mean(propfam)) %>% 
  group_by(pcaclust) %>% 
  slice_max(meanprop, n = 5)

synditab_prop_clades_clust %>% mutate(cladeprop = paste(gsub("Dino-Group", "DG", family))) %>% 
  #rowwise() %>% 
  mutate(propp = paste0("(",as.character(round(meanprop*100, 1)),")")) %>% 
  mutate(prop2 = paste(cladeprop, propp)) %>%
  print(n = 30)
  
  
synditab_3200_prop_long_wclust <- pivot_longer(synditab_3200_prop, cols = grp_df_3200$env_sample, names_to = "env_sample", values_to = "prop") %>% 
    left_join(., grp_df_3200, by = "env_sample")
  
  synditab_prop_clades_clust <- synditab_3200_prop_long_wclust %>% 
    group_by(family, pcaclust, env_sample) %>% 
    reframe(propfam = sum(prop)) %>%
    group_by(family, pcaclust) %>% 
    reframe(meanprop = mean(propfam)) %>% 
    group_by(pcaclust) %>% 
    slice_max(meanprop, n = 5)
  
  synditab_prop_clades_clust %>% mutate(cladeprop = paste(gsub("Dino-Group", "DG", family))) %>% 
    #rowwise() %>% 
    mutate(propp = paste0("(",as.character(round(meanprop*100, 1)),")")) %>% 
    mutate(prop2 = paste(cladeprop, propp)) %>%
    print(n = 30)
  
#### Read env_table (from Egge et al. 2021) ####
env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

#### Read files with indication of significantly differential abundance of ASVs (from script "mvabund_syndiniales.R) ####
signmulti_043 <- read_delim(here("data", "signmulti_043_feb24.txt"), delim = "\t")
signmulti_3200 <- read_delim(here("data", "signmulti_3200_feb24.txt"), delim = "\t")


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


#### centered-log transformation 0.4-3 ####
asv043_ps <- phylosize(synditab_043_rn_num, synditab_043_rn_tax, env_table)

syndiclr_043 <- otu_table(microbiome::transform(asv043_ps, transform = "clr", target = "OTU"))

syndiclr_043_df <- phyloseq_to_df(syndiclr_043, addtax =F)

#write_delim(syndiclr_043_df, here("data", "syndiclr_043_df.txt"), delim = "\t")

#### Aitchison distance = euclidean between clr-transformed samples ### c.f. Gloor et al. 2017
adist_043 <- vegdist(t(syndiclr_043), method = "euclidean")

#### Hierarchical clustering with complete linkage 0.4-3 ####
asvclust0_043 <- hclust(adist_043, method = "complete")


#### clr, 3-200 ####
asv3200_ps <- phylosize(synditab_3200_rn_num, synditab_3200_rn_tax, env_table)

syndiclr_3200 <- otu_table(microbiome::transform(asv3200_ps, transform = "clr", target = "OTU"))

syndiclr_3200_df <- phyloseq_to_df(syndiclr_3200, addtax =F)

#write_delim(syndiclr_3200_df, here("data", "syndiclr_3200_df.txt"), delim = "\t")

#### Aitchison distance = euclidean between clr-transformed samples ### c.f. Gloor et al. 2017
adist3200 <- vegdist(t(syndiclr_3200), method = "euclidean")

#### Hierarchical clustering with complete linkage, 3-200 ####
asvclust0_3200 <- hclust(adist3200, method = "complete")


sample_order_clustfig_043 <- grp_df_043$env_sample

writeLines(sample_order_clustfig_043, here("data", "sample_order_clustfig_043.txt"))
  

#### Heatmap with clr-transformed values 0.4-3 ####
#Include only ASVs with CLR>=10 in at least one sample in heatmap

syndiclr_abund_043 <- syndiclr_043
syndiclr_abund_043[syndiclr_abund_043<10] <- 0 #Exactly 20 ASVs

syndiclrabunddf_043 <- phyloseq_to_df(syndiclr_abund_043, addtax =F)
syndiclrabunddf_043$rowsums <- rowSums(syndiclr_abund_043)
syndiclrabunddf_043$asv_id <- rownames(syndiclr_043)
asv_abund_043 <- syndiclrabunddf_043 %>% filter(rowsums > 0) %>% pull(asv_id)
asv_abund_043
syndiclrdf_043 <- phyloseq_to_df(syndiclr_043, addtax =F)
syndiclrdf_sub_043 <- syndiclrdf_043 %>% filter(OTU %in% asv_abund_043)
syndiclrdf_sub_long_043 <- pivot_longer(syndiclrdf_sub_043, cols = - OTU, names_to = "env_sample", values_to = "clr_value")
str(syndiclrdf_sub_long_043)
tax_abund_043 <- phyloseq_to_df(asv043_ps, addtax = T) %>% filter(OTU %in% asv_abund_043) %>% dplyr::select(1:12) #Taxonomy file for the abundant ASVs

cluster_order_043 <- c("JanMar_Atl", "JanMar_Arc",  "May_epi", "Aug_epi_N03", "MayAuNo_meso", "1000m")

syndiclrdf_sub_long_forplot_043 <- left_join(syndiclrdf_sub_long_043, grp_df_043, by = "env_sample") %>% 
  mutate(cluster = factor(pcaclust, levels = cluster_order_043, ordered = T)) %>% left_join(., tax_abund_043, by = c("OTU" = "asv_id")) %>% 
  arrange(.data$family, .data$order, .data$OTU) %>%
  left_join(., signmulti_043, by = c("OTU" = "asv_id")) %>% 
  mutate(asv_id_star = case_when(!is.na(pval) ~ asv_id_star,
                                 TRUE ~ OTU)) %>% 
  mutate(asv_id_star  = str_replace(asv_id_star, "Clade", "C")) %>% 
  mutate(asv_id_star  = str_replace(asv_id_star, "Hematodinium-Group", "Hem-Gr")) %>% 
  mutate(env_sample = factor(env_sample, levels = sample_order_clustfig_043, ordered = T)) %>% 
  mutate(asv_id_star = factor(asv_id_star, levels = unique(asv_id_star))) %>% 
  mutate(cluster = recode(cluster, JanMar_Atl = "JanMar_\nAtl"))


#View(syndiclrdf_sub_long_forplot_043)

heatmap_asvs_043 <- ggplot(data = syndiclrdf_sub_long_forplot_043, aes(y = asv_id_star, x = env_sample, fill = clr_value))+
  geom_tile()+
  facet_grid(cols = vars(cluster), scales = "free", space = "free")+
  scale_y_discrete(limits = rev)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "right",
        panel.border = element_blank(),  
              # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
              # Remove panel background
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.01, 0.08, 0.1, "cm"),
        axis.text = element_text(size = 8),
        strip.background = element_rect(color = "black", fill = "white", linewidth = 0.2))+
  ylab("ASV")+
  xlab("Sample")+
  scale_fill_gradient2("CLR", low = "darkred", high = "darkblue", limits = c(-2, 12.5), breaks = c(-2,0,2,4,6,8,10,12))
heatmap_asvs_043
hm_legend <- get_legend(heatmap_asvs_043)

#### Barplot clades ####
bpcol <- read_delim(here("data","col_syndi_asv_both_fin.txt"), col_names = T, delim = "\t", comment = "")
bpcol <- bpcol %>% mutate_at(vars(family), list( ~factor(., levels = unique(bpcol$family), ordered = T)))
bpcolvec <- as.character(bpcol$color)
names(bpcolvec) <- bpcol$family


### Barplot at family level, 04-3 ####
taxlevel <- "family"

# Group all clades with less than 10% reads in all samples into "other"
lim <- 10
limfun <- function(x) {
  ifelse(x>=lim/100,1,0)
}

# Group by family
taxgroupspre_043 <-  synditab_043_prop %>% group_by(family) %>% summarise_if(is.numeric, sum)

# Apply "limfun" 
taxgroupspre_mat_043 <- as.matrix(taxgroupspre_043[,-1, drop = FALSE])
taxgroupspre_bin_043 <- apply(taxgroupspre_mat_043, 2, limfun)

if (is.vector(taxgroupspre_bin_043)) {
  taxgroupspre_bin2_043 <- as.data.frame(as.list(taxgroupspre_bin_043))
} else {
  taxgroupspre_bin2_043 <- as.data.frame(taxgroupspre_bin_043)}

# Filter the taxgroups accordingly
taxgroupspre_bin_sums_043 <- taxgroupspre_bin2_043 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
taxgroupspre_bin_sums_043$taxgroups <- taxgroupspre_043 %>% pull(taxlevel)
taxgroupspre_bin_sums_yes_043 <- filter(taxgroupspre_bin_sums_043, rowsumm >0) %>% dplyr::select(taxgroups)

# Create new data frames with taxgroups with > 10% reads in at least one sample, and the others
taxgroups_select_043 <- taxgroupspre_043 %>% filter(.data[[taxlevel]] %in% taxgroupspre_bin_sums_yes_043$taxgroups)
taxgroups_other_043 <- taxgroupspre_043 %>% filter(!.data[[taxlevel]] %in% taxgroupspre_bin_sums_yes_043$taxgroups)

# Create sum of low-abundant taxgroups, and vector of taxonomic group, including category "other"
if (dim(taxgroupspre_bin_sums_yes_043)[1] < dim(taxgroupspre_bin_sums_043)[1]) {
  Taxonomic_group_043 <- c(taxgroupspre_bin_sums_yes_043$taxgroups,"Other")
  othersum_043 <- colSums(taxgroups_other_043[,-1])
} else {
  Taxonomic_group_043 <- c(taxgroupspre_bin_sums_yes_043$taxgroups)
  othersum_043 <- NULL
}

#rbind into new data frame
taxgroups_select2_043 <- rbind(taxgroups_select_043[,-1],othersum_043)

#Add column with taxonomic group (i.e., "family"), + category "other"
taxgroups_select31_043 <- cbind.data.frame(Taxonomic_group_043,taxgroups_select2_043)

#pivot longer, mutate into factors where necessary, make taxon names shorter
taxgroups_select_longer_043 <- pivot_longer(taxgroups_select31_043, cols = -Taxonomic_group_043, names_to = "env_sample", values_to = "proportion") %>% 
  left_join(., grp_df_043, by = "env_sample") %>% 
  mutate(env_sample = factor(env_sample, levels = sample_order_clustfig_043, ordered = T)) %>%
  mutate(cluster = factor(pcaclust, levels = cluster_order_043, ordered = T)) %>% 
  mutate(Taxonomic_group_043 = gsub("Dino-Group", "DG", Taxonomic_group_043)) %>% 
  mutate_at(vars(Taxonomic_group_043), list( ~factor(., levels = unique(bpcol$family), ordered = T))) %>% 
  mutate(cluster = recode(cluster, JanMar_Atl = "JanMar_\nAtl"))


#### 3-200 NEW ####
sample_order_clustfig_3200 <- grp_df_3200$env_sample


writeLines(sample_order_clustfig_3200, here("data", "sample_order_clustfig_3200.txt"))


#### Heatmap with clr-transformed values 3-200 ####
#Include only ASVs with CLR>=10 in at least one sample in heatmap


syndiclr_abund_3200 <- syndiclr_3200
syndiclr_abund_3200[syndiclr_abund_3200<9] <- 0

syndiclrabunddf_3200 <- phyloseq_to_df(syndiclr_abund_3200, addtax =F)
syndiclrabunddf_3200$rowsums <- rowSums(syndiclr_abund_3200)
syndiclrabunddf_3200$asv_id <- rownames(syndiclr_3200)
asv_abund_3200 <- syndiclrabunddf_3200 %>% filter(rowsums > 0) %>% pull(asv_id)
asv_abund_3200
syndiclrdf_3200 <- phyloseq_to_df(syndiclr_3200, addtax =F)
syndiclrdf_sub_3200 <- syndiclrdf_3200 %>% filter(OTU %in% asv_abund_3200)
syndiclrdf_sub_long_3200 <- pivot_longer(syndiclrdf_sub_3200, cols = - OTU, names_to = "env_sample", values_to = "clr_value")
str(syndiclrdf_sub_long_3200)
tax_abund_3200 <- phyloseq_to_df(asv3200_ps, addtax = T) %>% filter(OTU %in% asv_abund_3200) %>% dplyr::select(1:12) #Taxonomy file for the abundant ASVs

cluster_order_3200 <- c("JanMar_Atl", "JanMar_Arc",   "May_epi", "Aug_epi", "AuNo_meso", "1000m")

syndiclrdf_sub_long_forplot_3200 <- left_join(syndiclrdf_sub_long_3200, grp_df_3200, by = "env_sample") %>% 
  mutate(cluster = factor(pcaclust, levels = cluster_order_3200)) %>% left_join(., tax_abund_3200, by = c("OTU" = "asv_id")) %>% 
  arrange(.data$family, .data$order, .data$OTU) %>%
  left_join(., signmulti_3200, by = c("OTU" = "asv_id")) %>% 
  mutate(asv_id_star = case_when(!is.na(pval) ~ asv_id_star,
                                 TRUE ~ OTU)) %>% 
  mutate(asv_id_star  = str_replace(asv_id_star, "Clade", "C")) %>% 
  mutate(asv_id_star  = str_replace(asv_id_star, "Hematodinium-Group", "Hem-Gr")) %>% 
  mutate(env_sample = factor(env_sample, levels = sample_order_clustfig_3200, ordered = T)) %>% 
  mutate(asv_id_star = factor(asv_id_star, levels = unique(asv_id_star))) %>% 
  mutate(cluster = recode(cluster, JanMar_Atl = "JanMar_\nAtl"))


#View(syndiclrdf_sub_long_forplot_3200)

heatmap_asvs_3200 <- ggplot(data = syndiclrdf_sub_long_forplot_3200, aes(y = asv_id_star, x = env_sample, fill = clr_value))+
  geom_tile()+
  facet_grid(cols = vars(cluster), scales = "free", space = "free")+
  scale_y_discrete(limits = rev)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing = unit(0.1, "lines"),
        #legend.position = "bottom",
        panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.01, 0.08, 0.1, "cm"),
        axis.text = element_text(size = 8),
        strip.background = element_rect(color = "black", fill = "white", size = 0.2))+
  ylab("ASV")+
  xlab("Sample")+
  scale_fill_gradient2("CLR", low = "darkred", high = "darkblue", limits = c(-2, 12.5))


heatmap_asvs_3200

#### Heatmap both size fractions ####
hm_both <- plot_grid(heatmap_asvs_043+theme(legend.position = "none"), heatmap_asvs_3200+theme(legend.position = "none"), 
                     ncol = 1, labels = c("A", "B"))
hm_both
hm_both_leg <- plot_grid(hm_both, hm_legend, rel_widths = c(1, 0.1))
hm_both_leg


ragg::agg_png("Fig_hm_both_mar24.png", width = 8.1, height = 8, units = "in", res = 300)
hm_both_leg
dev.off()


#### Barplot clades ####
### Barplot at family level, 3-200 ####
taxlevel <- "family"

# Group all clades with less than 10% reads in any sample into "other"
lim <- 10
limfun <- function(x) {
  ifelse(x>=lim/100,1,0)
}

# Group by family
taxgroupspre_3200 <-  synditab_3200_prop %>% group_by(family) %>% summarise_if(is.numeric, sum)

# Apply "limfun" 
taxgroupspre_mat_3200 <- as.matrix(taxgroupspre_3200[,-1, drop = FALSE])
taxgroupspre_bin_3200 <- apply(taxgroupspre_mat_3200, 2, limfun)

if (is.vector(taxgroupspre_bin_3200)) {
  taxgroupspre_bin2_3200 <- as.data.frame(as.list(taxgroupspre_bin_3200))
} else {
  taxgroupspre_bin2_3200 <- as.data.frame(taxgroupspre_bin_3200)}

# Filter the taxgroups accordingly
taxgroupspre_bin_sums_3200 <- taxgroupspre_bin2_3200 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
taxgroupspre_bin_sums_3200$taxgroups <- taxgroupspre_3200 %>% pull(taxlevel)
taxgroupspre_bin_sums_yes_3200 <- filter(taxgroupspre_bin_sums_3200, rowsumm >0) %>% dplyr::select(taxgroups)

# Create new data frames with taxgroups with > 10% reads in at least one sample, and the others
taxgroups_select_3200 <- taxgroupspre_3200 %>% filter(.data[[taxlevel]] %in% taxgroupspre_bin_sums_yes_3200$taxgroups)
taxgroups_other_3200 <- taxgroupspre_3200 %>% filter(!.data[[taxlevel]] %in% taxgroupspre_bin_sums_yes_3200$taxgroups)

# Create sum of low-abundant taxgroups, and vector of taxonomic group, including category "other"
if (dim(taxgroupspre_bin_sums_yes_3200)[1] < dim(taxgroupspre_bin_sums_3200)[1]) {
  Taxonomic_group_3200 <- c(taxgroupspre_bin_sums_yes_3200$taxgroups,"Other")
  othersum_3200 <- colSums(taxgroups_other_3200[,-1])
} else {
  Taxonomic_group_3200 <- c(taxgroupspre_bin_sums_yes_3200$taxgroups)
  othersum_3200 <- NULL
}

#rbind into new data frame
taxgroups_select2_3200 <- rbind(taxgroups_select_3200[,-1],othersum_3200)

#Add column with taxonomic group (i.e., "family"), + category "other"
taxgroups_select31_3200 <- cbind.data.frame(Taxonomic_group_3200,taxgroups_select2_3200)

#pivot longer, mutate into factors where necessary, make taxon names shorter
taxgroups_select_longer_3200 <- pivot_longer(taxgroups_select31_3200, cols = -Taxonomic_group_3200, names_to = "env_sample", values_to = "proportion") %>% 
  left_join(., grp_df_3200, by = "env_sample") %>% 
  mutate(env_sample = factor(env_sample, levels = sample_order_clustfig_3200, ordered = T)) %>%
  mutate(cluster = factor(pcaclust, levels = cluster_order_3200, ordered = T)) %>% 
  mutate(Taxonomic_group_3200 = gsub("Dino-Group", "DG", Taxonomic_group_3200)) %>% 
  mutate_at(vars(Taxonomic_group_3200), list( ~factor(., levels = unique(bpcol$family), ordered = T)))%>% 
  mutate(cluster = recode(cluster, JanMar_Atl = "JanMar_\nAtl"))

#### join plots in figure ####

#### Make legend family prop barpplot ####


bp_leg_plot <- ggplot(bpcol, aes(x = family, fill = family))+
  geom_bar()+
  scale_fill_manual(values = bpcolvec)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill = guide_legend(nrow = 7, byrow = F))+
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"))+
  labs(fill = "Family")
bp_leg_plot
bp_legnd <- get_legend(bp_leg_plot)

  
  barph = 0.48

#### Barplot ####
barplot_clust_043 <- ggplot(taxgroups_select_longer_043, aes(x=env_sample, y = proportion, fill = Taxonomic_group_043))+
  geom_bar(stat = "identity", position = "stack")+
  theme(panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(color = "black", fill = "white", size = 0.2))+
  facet_grid(cols = vars(cluster), scales = "free", space = "free")+
  scale_fill_manual(values = bpcolvec)+
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))+
  theme(strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        strip.text.y = element_text(margin = margin(0, 0.035, 0, 0.035, "cm")))+
  ylab("Proportion of reads")+
  #coord_flip()+
  guides(fill=guide_legend(title="Family"))+
  theme(legend.position = "right")+
  scale_y_continuous(labels = scales::percent)+
  xlab("Sample")#+
  #scale_x_discrete(limits = rev)
barplot_clust_043


barplot_clust_3200 <- ggplot(taxgroups_select_longer_3200, aes(x=env_sample, y = proportion, fill = Taxonomic_group_3200))+
  geom_bar(stat = "identity", position = "stack")+
  theme(panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(color = "black", fill = "white", size = 0.2))+
  facet_grid(cols = vars(cluster), scales = "free", space = "free")+
  scale_fill_manual(values = bpcolvec)+
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))+
  theme(strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        strip.text.y = element_text(margin = margin(0, 0.035, 0, 0.035, "cm")))+
  ylab("Proportion of reads")+
  #coord_flip()+
  guides(fill=guide_legend(title="Family"))+
  theme(legend.position = "right")+
  scale_y_continuous(labels = scales::percent)+
  xlab("Sample")#+
 # scale_x_discrete(limits = rev)
barplot_clust_3200


#### Make figure with pca-plot and barplots - pca plot from pc_axes_feb24.R ####
pca_barpl_plot <- pcaplot_043+pcaplot_3200+theme(legend.position = "none")+
  barplot_clust_043+theme(legend.position = "none")+barplot_clust_3200+theme(legend.position = "bottom")+
  plot_layout(ncol = 2, byrow = F, widths = c(1/3, 2/3))+plot_annotation(tag_levels = "A")
pca_barpl_plot

ragg::agg_png("Fig03_pca_bp_mar24.png", width = 12, height = 10, units = "in", res = 300)
pca_barpl_plot
dev.off()

