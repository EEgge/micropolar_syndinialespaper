library(mvabund)
library(compositions)
library(vegan)
library(magrittr)
library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)
library(metagMisc)

#### Read ASV tables prepared in script "prep_syndi_asvtable.R" ####
synditab_043_rn_num <- read_delim(here("data", "synditab_043_rn_num.txt"), delim = "\t")
synditab_043_rn_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")

synditab_3200_rn_num <- read_delim(here("data", "synditab_3200_rn_num_fromsubsamp.txt"), delim = "\t")
synditab_3200_rn_tax <- read_delim(here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")

#### Read env_table ####
env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

#### Read cluster grouping df, delimited by visual inspection of PCA + hclust ####
grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")


#### Join env table and group table (one for each size fraction) ####
env_table_wgrp_043 <- left_join(env_table0, grp_df_043, by = "env_sample")
env_table_wgrp_3200 <- left_join(env_table0, grp_df_3200, by = "env_sample")


#### Function for transforming ASV matrices to phyloseq objects ####
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


#### 0.4-3 ####
#### centered-log transformation 0.4-3###
asv043_ps <- phylosize(synditab_043_rn_num, synditab_043_rn_tax, env_table)
syndiclr_043 <- otu_table(microbiome::transform(asv043_ps, transform = "clr", target = "OTU"))

#### For mvabund-analyses, keep only ASVs where CLR value was above 0 in at least one sample within a give size fraction ####
syndiclr_abund_043 <- syndiclr_043
syndiclr_abund_043[syndiclr_abund_043 < 10] <- 0 # Set all values lower than the chosen limit to 0
syndiclrabunddf_043 <- phyloseq_to_df(syndiclr_abund_043, addtax =F) # Transform to data frame to be able to take rowsums
syndiclrabunddf_043$rowsums <- rowSums(syndiclr_abund_043)
syndiclrabunddf_043$asv_id <- rownames(syndiclr_043)
asv_abund_043 <- syndiclrabunddf_043 %>% filter(rowsums > 0) %>% dplyr::pull(asv_id) # Keep only asvs where clr value was above the limit in at least one sample
syndiclrdf_043 <- phyloseq_to_df(syndiclr_043, addtax =F) # Transform original clr-transformed to data frame (this function creates a column named 'OTU')
syndiclrdf_sub_043 <- syndiclrdf_043 %>% filter(OTU %in% asv_abund_043) # Keep OTUs/ASVs that passed the criterion
dim(syndiclrdf_sub_043)


syndiclrdft0_043 <- t(syndiclrdf_sub_043[, -1][env_table_wgrp_043$env_sample]) # to get samples in same order as in the env_table
syndiclrdft_043 <- as.data.frame(syndiclrdft0_043)

names(syndiclrdft_043) <- syndiclrdf_sub_043$OTU

syndimvabund_043 <- mvabund(syndiclrdft_043, neg = TRUE)

meanvar.plot(syndimvabund_043)

rownames(syndiclrdft_043) == env_table_wgrp_043$env_sample #check

#### Run manylm 10 times, keep ASVs with p <= 0.05 in 5 runs or more ####
aovlmclusterlist_043 <- list(NULL)
for (i in 1:10) {
  lmcluster_043 <- manylm(syndimvabund_043 ~ env_table_wgrp_043$pcaclust)
  aovlmclusterlist_043[[i]] <- anova(lmcluster_043, p.uni = "adjusted")
}
aovlmclusterpvallist_043 <- list(NULL)
for (i in 1:10) {
  aovlmclusterpvallist_043[[i]] <- data.frame(t(aovlmclusterlist_043[[i]]$uni.p))
  aovlmclusterpvallist_043[[i]]$asv_id <- rownames(aovlmclusterpvallist_043[[i]])
  names(aovlmclusterpvallist_043[[i]]) <- c("interc", "pval", "asv_id")
}

pvals_043 <- dplyr::bind_rows(aovlmclusterpvallist_043) %>% dplyr::filter(pval <= 0.05) %>% 
  group_by(asv_id) %>% summarise_if(is.numeric, mean)

signmulti_043 <- dplyr::bind_rows(aovlmclusterpvallist_043) %>% dplyr::filter(pval <= 0.05) %>% 
  group_by(asv_id) %>% count() %>% dplyr::filter(n >4) %>% 
  left_join(., pvals_043, by = "asv_id") %>%  
  mutate(asv_id = gsub("\\.", "-", asv_id))

View(signmulti_043)

#write_delim(signmulti_043, here("data", "signmulti_043_feb24.txt"), delim = "\t")

signmulti_043 <- read_delim(here("data", "signmulti_043_feb24.txt"), delim = "\t") %>% 
  mutate(asv_id_star = case_when(pval > 0.01 ~ paste0("*", asv_id),
                                 TRUE ~ paste0("**", asv_id))) #Modify ASV ids to indicate if they have significantly different CLR between sample clusters
#write_delim(signmulti_043, here("data", "signmulti_043_feb24.txt"), delim = "\t")


#### 3-200 ####
#### centered-log transformation 3-200###
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


syndiclrdft0_3200 <- t(syndiclrdf_sub_3200[, -1][env_table_wgrp_3200$env_sample]) # to get samples in same order as in the env_table
syndiclrdft_3200 <- as.data.frame(syndiclrdft0_3200)

names(syndiclrdft_3200) <- syndiclrdf_sub_3200$OTU

syndimvabund_3200 <- mvabund(syndiclrdft_3200, neg = TRUE)

meanvar.plot(syndimvabund_3200)

rownames(syndiclrdft_3200) == env_table_wgrp_3200$env_sample

lmcluster_3200 <- manylm(syndimvabund_3200 ~ env_table_wgrp_3200$pcaclust)
aovlmclust_3200 <- anova(lmcluster_3200, p.uni = "adjusted")

aovlmclustdf_3200 <- data.frame(t(aovlmclust_3200$uni.p))
aovlmclustdf_3200$asv_id <- rownames(aovlmclustdf_3200)
names(aovlmclustdf_3200) <- c("interc", "pval", "asv_id")
head(aovlmclustdf_3200)
aovlmclustdfsign_3200 <- aovlmclustdf_3200 %>% dplyr::filter(pval <= 0.05)

aovlmclustdfsign_3200_fin <- aovlmclustdfsign_3200 %>% mutate(asv_id = gsub("\\.", "-", asv_id)) %>% 
  mutate(asv_id_star = case_when(pval > 0.01 ~ paste0("*", asv_id),
                                 TRUE ~ paste0("**", asv_id)))
aovlmclustdfsign_3200_fin

#write_delim(aovlmclustdfsign_3200_fin, "data/mvabund_sign_3200_feb24.txt", delim = "\t")

#### Run manylm 10 times, keep ASVs with p <= 0.05 in 5 runs or more ####
aovlmclusterlist_3200 <- list(NULL)
for (i in 1:10) {
  lmcluster_3200 <- manylm(syndimvabund_3200 ~ env_table_wgrp_3200$pcaclust)
  aovlmclusterlist_3200[[i]] <- anova(lmcluster_3200, p.uni = "adjusted")
}
aovlmclusterpvallist_3200 <- list(NULL)
for (i in 1:10) {
  aovlmclusterpvallist_3200[[i]] <- data.frame(t(aovlmclusterlist_3200[[i]]$uni.p))
  aovlmclusterpvallist_3200[[i]]$asv_id <- rownames(aovlmclusterpvallist_3200[[i]])
  names(aovlmclusterpvallist_3200[[i]]) <- c("interc", "pval", "asv_id")
}

pvals_3200 <- dplyr::bind_rows(aovlmclusterpvallist_3200) %>% dplyr::filter(pval <= 0.05) %>% 
  group_by(asv_id) %>% summarise_if(is.numeric, mean)

signmulti_3200 <- dplyr::bind_rows(aovlmclusterpvallist_3200) %>% dplyr::filter(pval <= 0.05) %>% 
  group_by(asv_id) %>% count() %>% dplyr::filter(n >4) %>% 
  left_join(., pvals_3200, by = "asv_id") %>%  
  mutate(asv_id = gsub("\\.", "-", asv_id)) %>% 
  mutate(asv_id_star = case_when(pval > 0.01 ~ paste0("*", asv_id),
                                 TRUE ~ paste0("**", asv_id)))

View(signmulti_3200)
#write_delim(signmulti_3200, here("data", "signmulti_3200_feb24.txt"), delim = "\t")

#### Taxonomic distribution of differentially abundant ASVs ####

order_clade <- synditab_043_rn_tax %>% select(order, family, asv_id) %>% unique()

signmulti_043 <- read_delim(here("data", "signmulti_043_feb24.txt"), delim = "\t")
dim(signmulti_043)

tax_distr_sign043 <- left_join(signmulti_043, order_clade, by = "asv_id") %>% group_by(family) %>% count() 

signmulti_3200 <- read_delim(here("data", "signmulti_3200_feb24.txt"), delim = "\t")
signmulti_3200

tax_distr_sign3200 <- left_join(signmulti_3200, order_clade, by = "asv_id") %>% group_by(family) %>% count()

tax_distr_sign <- full_join(tax_distr_sign043, tax_distr_sign3200, by = "family") %>% rename(Family = family, Pico = n.x, `Nano-micro` = n.y) %>% 
  mutate(across(where(is.numeric), as.character)) %>% 
  tidyr::replace_na(list(Pico = "NA", `Nano-micro` = "NA"))
xtable::print.xtable(tax_distr_sign)
