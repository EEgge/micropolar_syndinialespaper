require(stringr)
require(formattable)
#install.packages("xtable")
require(xtable)
require(readr)
require(here)
require(dplyr)
require(magrittr)
require(tidyr)

customGreen0 = "#DeF7E9"

customGreen = "#71CA97"

customRed = "#ff7f7f"


taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")

#### Read ASV tables prepared in script "prep_syndi_asvtable.R" ####
synditab_043_rn_num <- read_delim(here("data", "synditab_043_rn_num.txt"), delim = "\t")
synditab_043_rn_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")

synditab_3200_rn_num <- read_delim(here("data", "synditab_3200_rn_num_fromsubsamp.txt"), delim = "\t")
synditab_3200_rn_tax <- read_delim(here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")

#### Read clr-transformed tables #### From script "cluster_barplot_figures_clr_new.R"
syndiclr_043_df <- read_delim(here("data", "syndiclr_043_df.txt"), delim = "\t")
syndiclr_043_wtax <- left_join(syndiclr_043_df, synditab_043_rn_tax, by = c("OTU"="asv_id")) %>% rename("asv_id" = "OTU")

syndiclr_3200_df <- read_delim(here("data", "syndiclr_3200_df.txt"), delim = "\t")
syndiclr_3200_wtax <- left_join(syndiclr_3200_df, synditab_3200_rn_tax, by = c("OTU"="asv_id")) %>% rename("asv_id" = "OTU")

syndiclr_abund <- syndiclr_3200
syndiclr_abund[syndiclr_abund < 0] <- 0 # Set all values lower than the chosen limit to 0
syndiclrabunddf <- phyloseq_to_df(syndiclr_abund, addtax =F) # Transform to data frame to be able to take rowsums
syndiclrabunddf$rowsums <- rowSums(syndiclr_abund)
syndiclrabunddf$asv_id <- rownames(syndiclr_abund)
asv_abund <- syndiclrabunddf %>% filter(rowsums > 0) %>% dplyr::pull(asv_id) # Keep only asvs where clr value was above the limit in at least one sample




#### Read env_table ####
env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

env_table <- env_table0 %>% 
  tibble::column_to_rownames("env_sample") 

#### Read cluster grouping df, from script "cluster_barplot_figures_clr_new.R" ####
grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")

#### Transform to proportions ####
synditab_043_prop <- sweep(synditab_043_rn_num, 2, colSums(synditab_043_rn_num), FUN = "/")
synditab_3200_prop <- sweep(synditab_3200_rn_num, 2, colSums(synditab_3200_rn_num), FUN = "/")

synditab_043_prop_wtax <- cbind.data.frame(synditab_043_prop, synditab_043_rn_tax)

synditab_043_prop_long <- pivot_longer(synditab_043_prop_wtax, cols = -c(all_of(taxlevels), asv_id), names_to = "sample", values_to = "prop")
synditab_043_prop_long_wgrp <- left_join(synditab_043_prop_long, grp_df_043, by = c("sample" = "env_sample")) %>% 
  mutate(cluster = as.factor(pcaclust))

synditab_3200_prop <- sweep(synditab_3200_rn_num, 2, colSums(synditab_3200_rn_num), FUN = "/")
synditab_3200_prop <- sweep(synditab_3200_rn_num, 2, colSums(synditab_3200_rn_num), FUN = "/")

synditab_3200_prop_wtax <- cbind.data.frame(synditab_3200_prop, synditab_3200_rn_tax)

synditab_3200_prop_long <- pivot_longer(synditab_3200_prop_wtax, cols = -c(all_of(taxlevels), asv_id), names_to = "sample", values_to = "prop")
synditab_3200_prop_long_wgrp <- left_join(synditab_3200_prop_long, grp_df_3200, by = c("sample" = "env_sample")) %>% 
  mutate(cluster = as.factor(pcaclust))

#### Make clr-transformed #### 
synditab_043_clr_long <- pivot_longer(syndiclr_043_wtax, cols = -c(all_of(taxlevels), asv_id), names_to = "sample", values_to = "clr")
synditab_043_clr_long_wgrp <- left_join(synditab_043_clr_long, grp_df_043, by = c("sample" = "env_sample")) %>% 
  mutate(cluster = as.factor(pcaclust))

synditab_3200_clr_long <- pivot_longer(syndiclr_3200_wtax, cols = -c(all_of(taxlevels), asv_id), names_to = "sample", values_to = "clr")
synditab_3200_clr_long_wgrp <- left_join(synditab_3200_clr_long, grp_df_3200, by = c("sample" = "env_sample")) %>% 
  mutate(cluster = as.factor(pcaclust))

#### Create table with most abundant clades in each cluster ####
clustertab_043 <- synditab_043_prop_long_wgrp %>% group_by(family, sample, cluster) %>% summarise_if(is.numeric, sum) %>%
  group_by(family, cluster) %>% mutate(maxprop = max(prop), meanprop = mean(prop), minprop = min(prop)) %>% 
  group_by(cluster, family) %>% summarise_if(is.numeric, mean) %>%  
  slice_max(order_by = prop, n = 5) %>%
  mutate(rangeprop = paste(round(meanprop*100,1), round(maxprop*100,1), sep = ", ")) %>%
  mutate(family = stringr::str_replace(family, "Dino-Group", "DG")) %>%
  mutate(famrange = paste0(family, " (", round(meanprop*100,1), ")")) %>% 
  summarise(test = toString(famrange))
clustertab_043
  
#xtable::print.xtable(clustertab_043)

clustertabasv_043 <- synditab_043_prop_long_wgrp %>% group_by(asv_id, sample, cluster) %>% summarise_if(is.numeric, sum) %>%
  group_by(asv_id, cluster) %>% mutate(maxprop = max(prop), meanprop = mean(prop), minprop = min(prop)) %>% 
  group_by(cluster, asv_id) %>% summarise_if(is.numeric, mean) %>%  
  slice_max(order_by = prop, n = 5) %>%
  mutate(rangeprop = paste(round(minprop*100,1), round(maxprop*100,1), sep = ", ")) %>%
  #mutate(family = stringr::str_replace(family, "Dino-Group", "DG")) %>%
  mutate(asvrange = paste0(asv_id, " (", rangeprop, ")")) %>% 
  summarise(test = toString(asvrange))

#xtable::print.xtable(clustertabasv_043)

clustertabasvclr_043 <- synditab_043_clr_long_wgrp %>% group_by(asv_id, sample, cluster) %>% summarise_if(is.numeric, sum) %>%
  group_by(asv_id, cluster) %>% mutate(maxprop = max(clr), meanprop = mean(clr), minprop = min(clr)) %>% 
  group_by(cluster, asv_id) %>% summarise_if(is.numeric, mean) %>%  
  slice_max(order_by = clr, n = 10) %>%
  #dplyr::filter(meanprop >= 8) %>% 
  arrange(desc(clr)) %>% 
  mutate(rangeprop = paste(round(minprop,1), round(maxprop,1), sep = ", ")) %>%
  mutate(asvrange = paste0(asv_id, " (", round(meanprop,1), ")")) %>% 
  summarise(test = toString(asvrange))

xtable::print.xtable(clustertabasvclr_043)

clustertabasvclr_043_wtax <- left_join(clustertabasvclr_043, synditab_043_rn_tax, by = "asv_id")

abundasvs_043 <- clustertabasvclr_043_wtax %>% group_by(clustername, family) %>% count() %>% pivot_wider(names_from = clustername, values_from = n, values_fill = 0) %>% 
  rename(Family = family) %>% dplyr::select(Family, phot_summer, aphot_winter, aphot_1000m, Jan_B08)

unit.scale_abund = function(x) (x/58)
abundasvs_043 <- formattable(abundasvs_043, align = c("l", "r", "r", "r", "r"), list(
  area(col = 2:5) ~ color_bar(customRed, fun = unit.scale_abund)
))

png(filename = "abundasvs_043.png", width = 1500, height = 1000,
    units = "px", pointsize = 8, bg = "white", res = 200,
    restoreConsole = TRUE)
abundasvs_043
dev.off()

clustertab_3200 <- synditab_3200_prop_long_wgrp %>% group_by(family, sample, cluster) %>% summarise_if(is.numeric, sum) %>%
  group_by(family, cluster) %>% mutate(maxprop = max(prop), meanprop = mean(prop), minprop = min(prop)) %>% 
  group_by(cluster, family) %>% summarise_if(is.numeric, mean) %>%  
  slice_max(order_by = prop, n = 5) %>%
  mutate(rangeprop = paste(round(minprop*100,1), round(maxprop*100,1), sep = ", ")) %>%
  mutate(family = stringr::str_replace(family, "Dino-Group", "DG")) %>%
  mutate(famrange = paste0(family, " (", round(meanprop*100,1), ")")) %>% 
  summarise(test = toString(famrange))

xtable::print.xtable(clustertab_3200)

clustertabasvclr_3200 <- synditab_3200_clr_long_wgrp %>% group_by(asv_id, sample, cluster) %>% summarise_if(is.numeric, sum) %>%
  group_by(asv_id, cluster) %>% mutate(maxprop = max(clr), meanprop = mean(clr), minprop = min(clr)) %>% 
  group_by(cluster, asv_id) %>% summarise_if(is.numeric, mean) %>%  
  slice_max(order_by = clr, n = 10) %>%
  #dplyr::filter(maxprop >= 10) %>% 
  arrange(desc(clr)) %>%
  mutate(rangeprop = paste(round(minprop,1), round(maxprop,1), sep = ", ")) %>%
  mutate(asvrange = paste0(asv_id, " (", round(meanprop,1), ")")) %>%
  summarise(test = toString(asvrange))

xtable::print.xtable(clustertabasvclr_3200)

clustertabasvclr_3200_wtax <- left_join(clustertabasvclr_3200, synditab_3200_rn_tax, by = "asv_id")

clustertabasvclr_3200_wtax %>% group_by(clustername, family) %>% count() %>% pivot_wider(names_from = clustername, values_from = n, values_fill = 0)

unit.scale_abund = function(x) (x/5)
formattable(abundasvs_043, align = c("l", "r", "r", "r", "r"), list(
  area(col = 2:5) ~ color_bar(customRed, fun = unit.scale_abund)
))



#### Table with overall mean, max and number of ASVs ####
synditab_043_prop_long_wgrp_pres <- synditab_043_prop_long_wgrp %>% dplyr::filter(prop >0) 
meanmax_043 <- synditab_043_prop_long_wgrp_pres %>% 
  group_by(family, sample) %>% summarise_if(is.numeric, sum) %>% group_by(family) %>% mutate(maxprop_043 = max(prop), meanprop_043 = mean(prop), minprop_043 = min(prop)) %>% 
  group_by(family) %>% summarise_if(is.numeric, mean)
  
  
  
uniquetab_043 <- unique(synditab_043_prop_long_wgrp_pres[c("family", "asv_id")]) %>% mutate(dummy_043 = 1) #%>% group_by(family) %>% summarise_if(is.numeric, sum)

left_join(meanmax_043, uniquetab_043, by = "family")

synditab_3200_prop_long_wgrp_pres <- synditab_3200_prop_long_wgrp %>% dplyr::filter(prop >0) 
meanmax_3200 <- synditab_3200_prop_long_wgrp_pres %>% 
  group_by(family, sample) %>% summarise_if(is.numeric, sum) %>% group_by(family) %>% mutate(maxprop_3200 = max(prop), meanprop_3200 = mean(prop), minprop_3200 = min(prop)) %>% 
  group_by(family) %>% summarise_if(is.numeric, mean)



uniquetab_3200 <- unique(synditab_3200_prop_long_wgrp_pres[c("family", "asv_id")]) %>% mutate(dummy_3200 = 1) #%>% group_by(family) %>% summarise_if(is.numeric, sum)

left_join(meanmax_043, uniquetab_043, by = "family")

asvprestab <- full_join(uniquetab_043, uniquetab_3200) 
asvprestab$dummy_043[is.na(asvprestab$dummy_043)] <- 0
asvprestab$dummy_3200[is.na(asvprestab$dummy_3200)] <- 0
asvprestab <- asvprestab %>% mutate(only043 = case_when(dummy_043 == 1 & dummy_3200 == 0 ~1,
                                                        TRUE ~ 0), 
                                    only3200 = case_when(dummy_043 == 0 & dummy_3200 == 1 ~1,
                                                        TRUE ~ 0)) %>% group_by(family) %>% summarise_if(is.numeric, sum)

meanmaxtab <- full_join(meanmax_043, meanmax_3200, by = "family")


order_clade <- synditab_043_rn_tax %>% select(order, family) %>% unique()
finalcladetable <- full_join(asvprestab, meanmaxtab, by = "family") %>% mutate(Both = dummy_043-only043) %>% mutate(Total = Both+only043+only3200) %>% 
  rename(`Only 0.4-3` = only043, `Only 3-200` = only3200, Family = family, 
         `Mean perc. 0.4-3` = meanprop_043, `Max perc. 0.4-3` = maxprop_043, `Mean perc. 3-200` = meanprop_3200, `Max perc. 3-200` = maxprop_3200) %>% 
  dplyr::select(Family, `Mean perc. 0.4-3`, `Max perc. 0.4-3`, `Mean perc. 3-200`, `Max perc. 3-200`, Both, `Only 0.4-3`, `Only 3-200`, Total) %>% 
  mutate(across(contains("perc"), ~round(.*100, 1))) %>% 
  mutate(minimum = pmax(`Mean perc. 0.4-3`, `Mean perc. 3-200`)) %>% 
  dplyr::filter(minimum >= 0.0 ) %>% dplyr::select(-minimum) #%>% 
  #mutate(across(where(is.numeric), as.character)) %>% 
  #mutate(Family = stringr::str_replace(Family, "Dino-Group", "DG"))
  
asvs_perorder <- left_join(finalcladetable, order_clade, by = c("Family" = "family")) %>% 
  mutate(dummy = 1) %>% 
  group_by(order) %>% summarise_if(is.numeric, sum)

colnames(finalcladetable) = gsub("perc. ", "", colnames(finalcladetable))

xtable::print.xtable(finalcladetable)


write_delim(finalcladetable, here("data", "finalcladetable.txt"), delim = "\t")
