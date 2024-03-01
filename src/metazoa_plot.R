library(readxl)
library(tidyr)
library(purrr)
library(here)

#### Phylosize function for metazoa ####
phylosize_metazoa <- function(asvtab, asvtax) {
  asvtab_mat <- as.matrix(asvtab)
  asvtax_mat <- as.matrix(asvtax)
  rownames(asvtab_mat) <- asvtax$asv_id
  rownames(asvtax_mat) <- asvtax$asv_id
  asv_otutab <- otu_table(asvtab_mat, taxa_are_rows = TRUE)
  asv_taxtab <- tax_table(asvtax_mat)
  asv_ps <- phyloseq(asv_otutab, asv_taxtab)
  return(asv_ps)
}


#### Read total ASV table and meta data (technical data) ####
asvtab_wtax <- read_excel(here("data", "metapr2_wide_asv_set_207_208_209_Eukaryota.xlsx"))
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles_forplot2.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

#### Filter Metazoa ####
asvtab_metazoa <- asvtab_wtax %>% dplyr::filter(division == "Metazoa")
dim(asvtab_metazoa)
asvtab_metazoa_num <- asvtab_metazoa %>% purrr::keep(is.numeric)
asvtab_zoo_pivot <- pivot_longer(asvtab_metazoa, cols=colnames(asvtab_metazoa_num), names_to = "seq_event")
dim(asvtab_zoo_pivot)

#### Merge replicate sequencing samples ####
asvtab_zoo_pivot_sampsf <- left_join(asvtab_zoo_pivot, meta_tab_seq_event, by = "seq_event")
asvtab_zoo_pivot_sampsf_merge <- asvtab_zoo_pivot_sampsf %>% group_by(sample_sizefract, asv_code, kingdom, supergroup, division, class, family, order, genus, species, sequence, collection_method, env_sample, fraction_max, fraction_min) %>% 
  summarise_if(is.numeric, sum) %>% ungroup()
dim(asvtab_zoo_pivot_sampsf_merge)

#check sample names
asvtab_zoo_pivot_sampsf_merge %>% dplyr::select(sample_sizefract) %>% unique()

# Select only the largest size fractions
asvtab_zoo_pivot_sampsf_merge_50200 <- asvtab_zoo_pivot_sampsf_merge %>% filter(fraction_max %in% c(200, 180), collection_method == "niskin")
asvtab_zoo_pivot_sampsf_merge_50200 %>% dplyr::select(sample_sizefract) %>% unique()

# Pivot wider, remove Homo sapiens and Capra (goat, reindeer)
merged_zoo_tab0 <- pivot_wider(asvtab_zoo_pivot_sampsf_merge_50200, id_cols = c(asv_code, kingdom, supergroup, division, class, family, order, genus, species, sequence), names_from = env_sample, values_from = value)
merged_zoo_tab <- merged_zoo_tab0 %>% filter(!(genus %in% c("Homo", "Capra")))

merged_zoo_tab_num <- merged_zoo_tab %>% purrr::keep(is.numeric)
merged_zoo_tab_num_asvid <- merged_zoo_tab_num %>% mutate(asv_code = merged_zoo_tab$asv_code)
merged_zoo_tab_num %>% dplyr::select(starts_with("May")) %>% summarise_if(is.numeric, sum)

merged_zoo_tab_notnum <- merged_zoo_tab %>% purrr::keep(negate(is.numeric))
merged_zoo_tab_prop0 <- sweep(merged_zoo_tab_num, 2, colSums(merged_zoo_tab_num), FUN = "/") %>% 
  mutate(rsum = rowSums(.))

# Arranging according to total abundance
merged_zoo_tab_prop <- cbind.data.frame(merged_zoo_tab_notnum, merged_zoo_tab_prop0) %>% arrange(desc(rsum)) %>% dplyr::select(-rsum)
merged_zoo_tab_prop <- merged_zoo_tab_prop %>% 
  mutate(asv_id = paste("ASVzoo", seq(1:dim(merged_zoo_tab_prop)[1]), species, sep = "_"))

# Arrange read number table according to total abundance based on proportions:
merged_zoo_tab_num_arr <- merged_zoo_tab_num_asvid %>% dplyr::arrange(factor(asv_code, levels = merged_zoo_tab_prop$asv_code))
merged_zoo_tab_notnum_arr <- merged_zoo_tab_notnum %>% dplyr::arrange(factor(asv_code, levels = merged_zoo_tab_prop$asv_code))
merged_zoo_tab_num_arrnum <- merged_zoo_tab_num_arr %>% dplyr::select(-asv_code)
#merged_zoo_tab_num_arrnum <- merged_zoo_tab_prop$asv_id
rownames(merged_zoo_tab_notnum_arr) <- merged_zoo_tab_prop$asv_id
merged_zoo_tab_notnum_arr$asv_id <- merged_zoo_tab_prop$asv_id
#write_delim(merged_zoo_tab_prop, "data/merged_zoo_tab_prop.txt", delim = "\t")

#### CLR-transform metazoa data ####
metazoa_ps <- phylosize_metazoa(merged_zoo_tab_num_arrnum, merged_zoo_tab_notnum_arr)

metazoa_clr <- otu_table(microbiome::transform(metazoa_ps, transform = "clr", target = "OTU"))

str(metazoa_ps)
str(metazoa_clr)

envsampnames <- c("Jan_B08_0001", "Jan_B08_0020", "Jan_B08_0500", "Jan_B08_1000", "Jan_B16_0001", "Jan_B16_0020", "Jan_B16_0500", "Jan_B16_1000", "Mar_M02_0001", "Mar_M02_0320",
                  "Mar_M02_1000","Mar_M03_0001", "Mar_M03_0020","Mar_M04_0001","Mar_M04_0020","Mar_M05_0020","Mar_M05_0120","Mar_M06_0020","May_P01_0001","May_P01_0020","May_P01_0417","May_P03_0001","May_P03_0015",
                  "May_P03_0447","May_P04_0001","May_P04_0015","May_P04_0500","May_P04_1000","Aug_P05_0001","Aug_P05_0020","Aug_P05_0213","Aug_P06_0001","Aug_P06_0024","Aug_P06_0500","Aug_P06_1000",
                  "Aug_P07_0001","Aug_P07_0025","Aug_P07_0500","Aug_P07_1000","Nov_N02_0020","Nov_N03_0020","Nov_N03_0300","Nov_N04_0020","Nov_N04_1000")

limfun <- function(x) {
  ifelse(x>=0.05,1,0)
}



taxlev_zoo <- "species"
groupby_taxlevel <- merged_zoo_tab_prop %>% group_by(.data[[taxlev_zoo]]) %>% summarise_if(is.numeric, sum)

taxgroupspre_mat <- as.matrix(groupby_taxlevel[,-1, drop = FALSE])
taxgroupspre_bin <- apply(taxgroupspre_mat, 2, limfun)

if (is.vector(taxgroupspre_bin)) {
  taxgroupspre_bin2 <- as.data.frame(as.list(taxgroupspre_bin))
} else {
  taxgroupspre_bin2 <- as.data.frame(taxgroupspre_bin)}

#### Identify taxa present as < selected % in all samples (low prop. taxa) ####
taxgroupspre_bin_sums <- taxgroupspre_bin2 %>% dplyr::select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
taxgroupspre_bin_sums$taxgroups <- groupby_taxlevel %>% pull(get(taxlev_zoo))
taxgroupspre_bin_sums_yes <- filter(taxgroupspre_bin_sums, rowsumm >0) %>% dplyr::select(taxgroups)

taxgroups_select <- groupby_taxlevel %>% filter(get(taxlev_zoo) %in% taxgroupspre_bin_sums_yes$taxgroups)
taxgroups_other <- groupby_taxlevel %>% filter(!get(taxlev_zoo) %in% taxgroupspre_bin_sums_yes$taxgroups)

if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
  othersum <- colSums(taxgroups_other[,-1])
} else {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
  othersum <- NULL
}


taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)

env_table


zoospec_long <- pivot_longer(taxgroups_select3, cols = -Taxonomic_group, names_to = "env_sample", values_to = "value")
zoospec_long <- zoospec_long %>% mutate(env_sample = factor(env_sample, levels = unique(envsampnames), ordered = T),
                                        month = factor(str_sub(env_sample, start = 1, end = 3), levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T))

zoospec_long2 <- left_join(zoospec_long, env_table, by = "env_sample") %>%  mutate(env_sample = factor(env_sample, levels = unique(envsampnames), ordered = T))


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)

zooplot <- ggplot(zoospec_long2, aes(x=reorder(env_sample, desc(env_sample)), y = value, fill = Taxonomic_group))+
  geom_bar(stat = "identity", position = "stack")+
  facet_grid(rows = vars(month.x), scales = "free")+
  scale_fill_manual(values = c25[c(1:19)])+
  xlab("")+
  ylab("Proportion of metazooan reads")+
  coord_flip()
ggplotly(zooplot)


