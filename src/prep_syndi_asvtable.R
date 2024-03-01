require("here")
require(readr)
require(tidyverse)
require(vegan)
require(purrr)

#### Read asv tables from Egge et al. 2021 ####
asvtab2 <- read_delim(here("data", "asvtab2_merged_readnum.txt"), delim = "\t") # To select the 0.4-3 fraction
asvtab4 <- read_delim(here("data", "asvtab4_merged_subsamp_readnum.txt"), delim = "\t") # Using the subsampled-to-equal read number within size fractions table to get (almost) equal representation of the fractions for each sample.
asvtab5 <- read_delim(here("data", "asvtab5_merged_subsamp_prop.txt"), delim = "\t") #To order Syndiniales ASVs after total abundance

#### Make ASV-ids according to abundance, and with clade name ####
syndiniales_prop <- asvtab5 %>% filter(divisionlong == "Dinoflagellata_Syndiniales")

head(syndiniales_prop)

syndiniales_prop_num <- syndiniales_prop %>% purrr::keep(is.numeric)
syndiniales_prop_notnum <- syndiniales_prop %>% purrr::keep(negate(is.numeric)) %>% 
  purrr::map_df(., ~gsub("Dino-Group", "DG", .x))
# Make new proportions of total Syndiniales
syndiniales_prop_prop <- sweep(syndiniales_prop_num, 2, colSums(syndiniales_prop_num), FUN = "/")
syndiniales_prop_new <- cbind.data.frame(syndiniales_prop_notnum, syndiniales_prop_prop)

syndiniales_prop_new <- syndiniales_prop_new %>% 
  mutate(total=rowSums(syndiniales_prop_new %>% dplyr::select_if(is.numeric))) %>% 
  arrange(desc(total)) %>% mutate(asv_id = paste("ASV", seq(1:dim(syndiniales_prop_new)[1]), family, sep = "_")) %>% 
  dplyr::select(-total)

asv_code_id <- syndiniales_prop_new %>% dplyr::select(asv_code, asv_id)

top100_syndiasvsseq <- syndiniales_prop_new %>% dplyr::select(asv_id, sequence) %>% slice(c(1:200))

#write_delim(top100_syndiasvsseq, here("data", "top200syndiasvs.txt"), delim = "\t")


## meta data for each sequencing_event, file from Egge et al. 2021 ##
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

meta_tab_sample_sf_w3_200 <- meta_tab_sample_sf %>% mutate(new_sf = ifelse(size_fraction %in% c("3_10", "10_50", "50_200"), "3_200", 
                                                                           ifelse(size_fraction == "0.4_3", "0.4_3", "3_180"))) #creating meta table with 3-200 size fraction

## Environmental data files ##
env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

taxlevels <- c("asv_id", "asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")

#### Prepare asv table with read counts for 0.4-3 ####
syndiniales_rn_nonsubs <- asvtab2 %>% filter(divisionlong == "Dinoflagellata_Syndiniales") %>% #Adding column with syndiniales asv_id
  left_join(., asv_code_id, by = "asv_code")

synditab_long <- pivot_longer(syndiniales_rn_nonsubs, cols = -all_of(taxlevels), names_to = "sample_sizefract")
synditab_long_meta <- left_join(synditab_long, meta_tab_sample_sf_w3_200, by = "sample_sizefract")

synditab_043_rn_long <- synditab_long_meta %>% dplyr::filter(size_fraction == "0.4_3")
synditab_043_rn <- pivot_wider(synditab_043_rn_long, id_cols = all_of(taxlevels), names_from = env_sample, values_from = value)
synditab_043_rn_num <- synditab_043_rn %>% purrr::keep(is.numeric) %>% as.data.frame()
rownames(synditab_043_rn_num) <- synditab_043_rn$asv_id
synditab_043_rn_tax <- synditab_043_rn %>% dplyr::select(all_of(taxlevels)) %>% as.data.frame()
rownames(synditab_043_rn_tax) <- synditab_043_rn$asv_id
#write_delim(synditab_043_rn_num, "data/synditab_043_rn_num.txt", delim = "\t")
#write_delim(synditab_043_rn_tax, "data/synditab_043_rn_tax.txt", delim = "\t")

#### ASV table subsampled to equal Syndiniales read number 0.4-3 ####
colSums(synditab_043_rn_num)
min(colSums(synditab_043_rn_num)) #2638

synditab_043_raref <- as.data.frame(t(rrarefy(t(synditab_043_rn_num), 5000)))
#write_delim(synditab_043_raref, here("data", "synditab_043_raref.txt"), delim = "\t")

#### 3-200, readnum, from subsampled to equal read number within each size fraction ####
syndiniales_readnum_for3200 <- asvtab4 %>% filter(divisionlong == "Dinoflagellata_Syndiniales") %>% 
  left_join(., asv_code_id, by = "asv_code")

synditab_long_for3200 <- pivot_longer(syndiniales_readnum_for3200, cols = -all_of(taxlevels), names_to = "sample_sizefract")
synditab_long_meta_for3200 <- left_join(synditab_long_for3200, meta_tab_sample_sf_w3_200, by = "sample_sizefract")

synditab_3200_rn_long <- synditab_long_meta_for3200 %>% dplyr::filter(collection_method == "niskin") %>% 
  group_by(asv_code, kingdom, supergroup, division, divisionlong, class, family, order, genus, species, sequence, asv_id, new_sf, env_sample) %>% 
  summarise_if(is.numeric, sum) %>% dplyr::filter(new_sf %in% c("3_180", "3_200"))
synditab_3200_rn <- pivot_wider(synditab_3200_rn_long, id_cols = all_of(taxlevels), names_from = env_sample, values_from = value)
synditab_3200_rn_num <- synditab_3200_rn %>% purrr::keep(is.numeric) %>% as.data.frame()
rownames(synditab_3200_rn_num) <- synditab_3200_rn$asv_id
synditab_3200_rn_tax <- synditab_3200_rn %>% dplyr::select(all_of(taxlevels)) %>% ungroup() %>% as.data.frame()
rownames(synditab_3200_rn_tax) <- synditab_3200_rn$asv_id
#write_delim(synditab_3200_rn_num, here("data", "synditab_3200_rn_num_fromsubsamp.txt"), delim = "\t")
#write_delim(synditab_3200_rn_tax, here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")


#### ASV table subsampled to equal Syndiniales read number 3-200 ####
colSums(synditab_3200_rn_num)
min(colSums(synditab_3200_rn_num)) #1514 hmmm

synditab_3200_raref <- as.data.frame(t(rrarefy(t(synditab_3200_rn_num), 5000)))
#write_delim(synditab_3200_raref, here("data", "synditab_3200_raref.txt"), delim = "\t")






