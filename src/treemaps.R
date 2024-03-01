require(stringr)
require(formattable)
#install.packages("xtable")
require(xtable)
require(readr)
require(here)
require(dplyr)
require(magrittr)
require(tidyr)
library(treemapify)
require(ggplot2)
require(ggpubr)
library(cowplot)
require(xtable)


#### Read color palette #### 
bpcol <- read_delim(here("data","col_syndi_asv_both_fin.txt"), col_names = T, delim = "\t", comment = "")
bpcol <- bpcol %>% mutate_at(vars(family), list( ~factor(., levels = unique(bpcol$family), ordered = T)))
bpcolvec <- as.character(bpcol$color)
names(bpcolvec) <- bpcol$family



#### Treemap ASVs ####
# Read clade table made in script "clade_tables.R"
cladetable <- read_delim(here("data", "finalcladetable.txt"), delim = "\t")
order_fam <- read_delim(here("data", "order_family.txt"), delim = "\t")

# Read order table of pr2 reference sequences
pr2_ref_order <- read_delim(here("data", "pr2_syndi_ref_order.txt"), delim = "\t") %>% 
  mutate(`PR2 ref.\nsequences` = paste0(dummy, " (", round(perc*100,1), "%)")) %>% 
  dplyr::select(-dummy, -perc)
pr2_ref_order

# Add order column
cladetable_wtax <- left_join(cladetable, order_fam, by = c("Family" = "family")) %>%
  rename(family = Family) %>% 
  mutate(family = gsub("Dino-Group", "DG", family)) %>% 
  mutate(nASVs043 = Both+`Only 0.4-3`) %>% 
  mutate(nASVs3200 = Both+`Only 3-200`)

names(cladetable_wtax)

fortreemap_numasvs <- cladetable_wtax %>% 
  mutate(abundclades = case_when(family %in% names(bpcolvec)~family,
                                 TRUE ~ "Other")) %>% 
  group_by(order, abundclades) %>% summarise_if(is.numeric, sum) %>% 
  mutate(label043nasvs = paste0(abundclades, "\n", nASVs043)) %>% 
  mutate(label3200nasvs = paste0(abundclades, "\n", nASVs3200)) %>% 
  mutate(label043mean = paste0(abundclades, "\n", `Mean 0.4-3`,"%")) %>% 
  mutate(label3200mean = paste0(abundclades, "\n", `Mean 3-200`,"%")) %>% 
  mutate(across(starts_with("label"),  ~gsub("Hematodinium-Group", "Hem.-Gr.",.)))
  

View(fortreemap_numasvs)

treem_nASVs_043 <- ggplot(fortreemap_numasvs, aes(area = nASVs043, fill = abundclades, label = label043nasvs, subgroup = order))+
  geom_treemap()+
  scale_fill_manual(values = bpcolvec)+
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15,
                    grow = TRUE)+
  geom_treemap_subgroup_border(colour = "white", size = 5)+
  labs(fill = "Family", title = "A") +
  theme(plot.title = element_text(face = "bold", size = 14))
treem_nASVs_043

treem_nASVs_3200 <- ggplot(fortreemap_numasvs, aes(area = nASVs3200, fill = abundclades, label = label3200nasvs, subgroup = order))+
  geom_treemap()+
  scale_fill_manual(values = bpcolvec)+
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15,
                    grow = TRUE)+
  geom_treemap_subgroup_border(colour = "white", size = 5)+
  labs(fill = "Family", title = "B") +
  theme(plot.title = element_text(face = "bold", size = 14))
treem_nASVs_3200

treem_prop_043 <- ggplot(fortreemap_numasvs, aes(area = `Mean 0.4-3`, fill = abundclades, label = label043mean, subgroup = order))+
  geom_treemap()+
  scale_fill_manual(values = bpcolvec)+
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15,
                    grow = TRUE)+
  geom_treemap_subgroup_border(colour = "white", size = 5)+
  labs(fill = "Family", title = "C") +
  theme(plot.title = element_text(face = "bold", size = 14))

treem_prop_3200 <- ggplot(fortreemap_numasvs, aes(area = `Mean 3-200`, fill = abundclades, label = label3200mean, subgroup = order))+
  geom_treemap()+
  scale_fill_manual(values = bpcolvec)+
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15,
                    grow = TRUE)+
  geom_treemap_subgroup_border(colour = "white", size = 5)+
  labs(fill = "Family", title = "D") +
  theme(plot.title = element_text(face = "bold", size = 14))
treem_prop_3200


treemap_grid <- ggarrange(treem_nASVs_043, treem_nASVs_3200, treem_prop_043, treem_prop_3200, ncol = 2, nrow = 2, 
          common.legend = T, legend = "bottom")

# ragg::agg_png("Fig02_treemaps.png", width = 15, height = 12, units = "cm", res = 300, scaling = 0.5)
# treemap_grid
# dev.off()


nasvs_vsprop_043 <- ggplot(data = fortreemap_numasvs, aes(x = nASVs043, y = `Mean 0.4-3`, color = abundclades))+
  geom_point()+
  theme_cowplot(font_size = 10)+
  scale_color_manual(values = bpcolvec)+
  labs(color = "Family")
nasvs_vsprop_043

nasvs_vsprop_3200 <- ggplot(data = fortreemap_numasvs, aes(x = nASVs3200, y = `Mean 3-200`, color = abundclades))+
  geom_point()+
  theme_cowplot(font_size = 10)+
  scale_color_manual(values = bpcolvec)+
  labs(color = "Family")
nasvs_vsprop_3200


nasvnasv <- ggplot(data = fortreemap_numasvs, aes(x = nASVs043, y = nASVs3200, color = abundclades))+
  geom_point()+
  theme_cowplot(font_size = 10)+
  scale_color_manual(values = bpcolvec)+
  labs(color = "Family")
nasvnasv


meanmean <- ggplot(data = fortreemap_numasvs, aes(x = `Mean 0.4-3`, y = `Mean 3-200`, color = abundclades))+
  geom_point()+
  theme_cowplot(font_size = 10)+
  scale_color_manual(values = bpcolvec)+
  labs(color = "Family")
meanmean



fortreemap_numasvs_noother <- fortreemap_numasvs %>% dplyr::filter(abundclades != "Other")

summary(lm(`Mean 0.4-3`~ nASVs043, fortreemap_numasvs_noother))
summary(lm(`Mean 3-200`~ nASVs3200, fortreemap_numasvs_noother))
summary(lm(nASVs043 ~ nASVs3200, fortreemap_numasvs_noother))
summary(lm(`Mean 0.4-3`~ `Mean 3-200`, fortreemap_numasvs_noother))

cor.test(fortreemap_numasvs_noother$`Mean 0.4-3`, fortreemap_numasvs_noother$nASVs043)
summary(lm(`Mean 3-200`~ nASVs3200, fortreemap_numasvs_noother))
cor.test(fortreemap_numasvs_noother$nASVs043, fortreemap_numasvs_noother$nASVs3200)
summary(lm(`Mean 0.4-3`~ `Mean 3-200`, fortreemap_numasvs_noother))

order_tab_num <-  fortreemap_numasvs %>% 
  mutate(order = gsub("-Hematodinium-Group", "", order)) %>%
  mutate(order = gsub("-Syndinium-Group", "", order)) %>%
  group_by(order) %>% summarise_if(is.numeric, sum) %>% 
  dplyr::select(-starts_with("Max")) %>% 
  mutate(nASVs043perc = round((nASVs043/sum(nASVs043))*100, 1)) %>% 
  mutate(nASVs3200perc = round((nASVs3200/sum(nASVs3200))*100, 1)) %>% 
  mutate(totalperc = round((Total/sum(Total))*100, 1))

order_tab <- order_tab_num %>% 
  mutate(`nASVs pico` = paste0(nASVs043," (", nASVs043perc, "%)")) %>% 
  mutate(`nASVs nano-\nmicro` = paste0(nASVs3200," (", nASVs3200perc, "%)")) %>% 
  mutate(`nASVs total` = paste0(Total, " (", totalperc, "%)")) %>% 
  dplyr::select(-starts_with(c("nASVs0", "nASVs3"))) %>% 
  dplyr::select(-Total, -totalperc) %>% 
  rename(`Mean reads (%),\n pico` = `Mean 0.4-3`) %>% 
  rename(`Mean reads (%),\n nano-micro` = `Mean 3-200`) %>% 
  left_join(pr2_ref_order, by = "order") %>% 
  rename(Order = order) %>% 
  mutate(across(.cols = everything(), as.character))

sumvec <- c("Sum", unname(round(colSums(order_tab_num[,c(2:12)]))))
names(sumvec) <- names(order_tab)
  
order_tab2 <- bind_rows(order_tab, sumvec)
View(order_tab2)

write_delim(order_tab2, here("data", "order_tab_wpr2.txt"), delim = "\t") #needs edits in excel
ordertab_fin <- read_delim(here("data", "order_tab_wpr2_edited.txt"), delim = "\t")
View(ordertab_fin)
xtable::print.xtable(ordertab_fin)


family_tab <-  fortreemap_numasvs %>% group_by(abundclades) %>% summarise_if(is.numeric, sum) %>% 
  dplyr::select(-starts_with("Max")) %>% 
  mutate(nASVs043perc = round((nASVs043/sum(nASVs043))*100, 1)) %>% 
  mutate(nASVs3200perc = round((nASVs3200/sum(nASVs3200))*100, 1)) %>% 
  mutate(totalperc = round((Total/sum(Total))*100, 1))
View(family_tab)
