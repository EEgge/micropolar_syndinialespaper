library(grid)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(readr)
library(ggplot2)
library(ggpubr)
library(here)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)

#### Read clr-transformed data made in script "cluster_barplot_figures_clr_new.R" ####
clr_043 <- read_delim(here("data", "syndiclr_043_df.txt"), delim = "\t")
clr_3200 <- read_delim(here("data", "syndiclr_3200_df.txt"), delim = "\t")

dim(clr_043)

synditab_043_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")
synditab_3200_tax <- read_delim(here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")


#### Read table with significant ASVs from "mvabund_syndiniales.R" ####
sign_043 <- read_delim("data/signmulti_043_feb24.txt", delim = "\t")
sign_3200 <- read_delim("data/signmulti_3200_feb24.txt", delim = "\t")

clr_sign_043 <- left_join(sign_043, clr_043, by = c("asv_id" = "OTU")) %>% 
  dplyr::select(-interc, -pval, -n) %>% left_join(., synditab_043_tax, by = "asv_id")

head(clr_sign_043)

clr_sign_3200 <- left_join(sign_3200, clr_3200, by = c("asv_id" = "OTU")) %>% 
  dplyr::select(-interc, -pval, -n) %>% left_join(., synditab_3200_tax, by = "asv_id")

#### Read grouping table from "cluster_barplot_figures_clr_new.R"
samp_ord_043 <- readLines(here("data", "sample_order_clustfig_043.txt"))
samp_ord_3200 <- readLines(here("data", "sample_order_clustfig_3200.txt"))

grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")

cluster_order_043 <- c("JanMar_Atl", "JanMar_Arc", "May_epi", "Aug_epi_N03", "MayAuNo_meso", "1000m")
cluster_order_3200 <- c("JanMar_Atl", "JanMar_Arc",   "May_epi", "Aug_epi", "AuNo_meso", "1000m")

taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")
#### Heatmap with clr-transformed values ####


clr_sign_long_043 <- pivot_longer(clr_sign_043, cols = -c(asv_id, asv_id_star, all_of(taxlevels)), names_to = "env_sample", values_to = "clr_value")

clr_sign_long_wgrp_043 <- left_join(clr_sign_long_043, grp_df_043, by = "env_sample") %>% 
  arrange(.data$family, .data$order, .data$asv_id) %>% mutate(asv_id  = str_replace(asv_id, "Clade", "C")) %>% 
  mutate(pcaclust = factor(pcaclust, levels = cluster_order_043, ordered = T)) %>% 
  mutate(asv_id  = str_replace(asv_id, "Hematodinium-Group", "Hem-Gr")) %>% 
  mutate(asv_id = factor(asv_id, levels = unique(asv_id))) %>% 
  mutate(env_sample = factor(env_sample, levels = samp_ord_043, ordered = T)) %>% 
  mutate(pcaclust = recode(pcaclust, JanMar_Atl = "JanMar_\nAtl"))

heatmap_asvs_043 <- ggplot(data = clr_sign_long_wgrp_043, aes(x = asv_id, y = env_sample, fill = clr_value))+
  geom_tile()+
  facet_grid(rows = vars(pcaclust), scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 0.8),
        panel.spacing = unit(0.1, "lines"))+
  scale_fill_gradient2(low = "darkred", high = "darkblue")


heatmap_asvs_043

heatmap_asvs_043_2 <- ggplot(data = clr_sign_long_wgrp_043, aes(y = asv_id, x = env_sample, fill = clr_value))+
  geom_tile()+
  facet_grid(cols = vars(pcaclust), scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 6))+
  scale_fill_gradient2(low = "darkred", high = "darkblue")+
  scale_y_discrete(limits = rev)+
  ylab("ASV")+
  xlab("Sample")+
  scale_fill_gradient2("CLR", low = "darkred", high = "darkblue", limits = c(-2, 12.5))
heatmap_asvs_043_2

#### 3-200 ####

clr_sign_long_3200 <- pivot_longer(clr_sign_3200, cols = -c(asv_id, asv_id_star, all_of(taxlevels)), names_to = "env_sample", values_to = "clr_value")

clr_sign_long_wgrp_3200 <- left_join(clr_sign_long_3200, grp_df_3200, by = "env_sample") %>% 
  arrange(.data$family, .data$order, .data$asv_id) %>% mutate(asv_id  = str_replace(asv_id, "Clade", "C")) %>% 
  mutate(cluster = factor(pcaclust, levels = cluster_order_3200, ordered = T)) %>% 
  mutate(asv_id  = str_replace(asv_id, "Hematodinium-Group", "Hem-Gr")) %>% 
  mutate(asv_id = factor(asv_id, levels = unique(asv_id))) %>% 
  mutate(env_sample = factor(env_sample, levels = samp_ord_3200, ordered = T)) %>% 
  mutate(cluster = recode(cluster, JanMar_Atl = "JanMar_\nAtl"))

heatmap_asvs_3200 <- ggplot(data = clr_sign_long_wgrp_3200, aes(x = asv_id, y = env_sample, fill = clr_value))+
  geom_tile()+
  facet_grid(rows = vars(cluster), scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 0.8),
        panel.spacing = unit(0.1, "lines"))+
  scale_fill_gradient2(low = "darkred", high = "darkblue")


heatmap_asvs_3200

heatmap_asvs_3200_2 <- ggplot(data = clr_sign_long_wgrp_3200, aes(y = asv_id, x = env_sample, fill = clr_value))+
  geom_tile()+
  facet_grid(cols = vars(cluster), scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        panel.spacing = unit(0.1, "lines"),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        strip.text.x = element_text(size = 8))+
  scale_fill_gradient2(low = "darkred", high = "darkblue")+
  scale_y_discrete(limits = rev)+
  xlab("Sample")+
  scale_fill_gradient2("CLR", low = "darkred", high = "darkblue", limits = c(-2, 12.5))
heatmap_asvs_3200_2

hm_sign_legnd <- get_legend(heatmap_asvs_043_2)


hm_sign_plot <- plot_grid(heatmap_asvs_043_2+theme(legend.position = "none", axis.title.y = element_blank(),
                                   axis.text.y = element_text(size = 5)), 
          heatmap_asvs_3200_2+theme(legend.position = "none", axis.title.y = element_blank(),
                                    axis.text.y = element_text(size = 5)),
          labels = "AUTO")

hm_sign_plot_wleg <- plot_grid(hm_sign_plot, hm_sign_legnd, ncol = 2, rel_widths = c(1,0.06))
hm_sign_plot_wleg

ragg::agg_png("hm_sign_plot_mar24.png", width = 15, height = 15, units = "cm", res = 300, scaling = 0.5)
hm_sign_plot_wleg
dev.off()

