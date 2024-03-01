library(vegan)
library(ggridges)
library(cowplot)
library(readr)
library(here)
library(magrittr)
library(dplyr)
library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)

# Load subsampled tables prepared in script "prep_syndi_asvtable.R"

synditab_043_raref <- read_delim(here("data", "synditab_043_raref.txt"), delim = "\t")
synditab_043_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")

synditab_3200_raref <- read_delim(here("data", "synditab_3200_raref.txt"), delim = "\t")
synditab_3200_tax <- read_delim(here("data", "synditab_3200_rn_tax.txt"), delim = "\t")

# Load sample group tables
grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")

cluster_order_043 <- c("JanMar_Atl", "JanMar_Arc", "May_epi", "Aug_epi_N03", "MayAuNo_meso", "1000m")
cluster_order_3200 <- c("JanMar_Atl", "JanMar_Arc",   "May_epi", "Aug_epi", "AuNo_meso", "1000m")



#### Read env_table ####
env_table <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t") %>% 
  mutate(Dato = as.Date(as.character(`yyyy-mm-dd`), format = '%d/%m/%y'))

rich_043 <- specnumber(t(synditab_043_raref))
shan_043 <- exp(vegan::diversity(t(synditab_043_raref)))
hill_043 <- shan_043/rich_043

rich_3200 <- specnumber(t(synditab_3200_raref))
shan_3200 <- exp(vegan::diversity(t(synditab_3200_raref)))
hill_3200 <- shan_3200/rich_3200


divdf04 <- tibble(env_sample = names(rich_043), numasvs = unname(rich_043), shan = unname(shan_043), hill = unname(hill_043))
divdf04_2 <- tibble(env_sample = rep(names(rich_043), 3), value = c(unname(rich_043), unname(shan_043), unname(hill_043)), 
                  param = c(rep("numasvs", 44), rep("shan", 44), rep("hill", 44)))
difdfenv04_2 <- left_join(divdf04_2, env_table, by = "env_sample") %>% 
  dplyr::select(env_sample, value, param, month, depthbin, station) %>%
  left_join(., grp_df_043, by = "env_sample") %>%
  mutate(cluster = factor(pcaclust, levels = cluster_order_043, ordered = TRUE))


alpha_summary_043 <- difdfenv04_2 %>% group_by(pcaclust, param) %>%
  #mutate(percval = value/max(value)) %>% 
  summarise(across(.cols = c(value), .fns = list(min = min, max = max))) %>% ungroup() %>%
  group_by(param) %>% 
  mutate(percval_min = round(100*value_min/max(value_max),1)) %>% 
  mutate(percval_max = round(100*value_max/max(value_max),1))
alpha_summary_043

divdf3 <- tibble(env_sample = names(rich_3200), numasvs = unname(rich_3200), shan = unname(shan_3200), hill = unname(hill_3200))
divdf_3200_2 <- tibble(env_sample = rep(names(rich_3200), 3), value = c(unname(rich_3200), unname(shan_3200), unname(hill_3200)), 
                    param = c(rep("numasvs", 44), rep("shan", 44), rep("hill", 44)))
difdfenv_3200_2 <- left_join(divdf_3200_2, env_table, by = "env_sample") %>% dplyr::select(env_sample, value, param, month, depthbin) %>%
  left_join(., grp_df_3200, by = "env_sample") %>%
  mutate(cluster = factor(pcaclust, levels = cluster_order_3200, ordered = TRUE))

alpha_summary_3200 <- difdfenv_3200_2 %>% group_by(pcaclust, param) %>%
  #mutate(percval = value/max(value)) %>% 
  summarise(across(.cols = c(value), .fns = list(min = min, max = max))) %>% ungroup() %>%
  group_by(param) %>% 
  mutate(percval_min = round(100*value_min/max(value_max),1)) %>% 
  mutate(percval_max = round(100*value_max/max(value_max),1))



divdfenv04 <- left_join(divdf04, env_table, by = "env_sample") %>% dplyr::select(env_sample, numasvs, shan, hill, month, depthbin, station) %>%
  left_join(., grp_df_043, by = "env_sample") %>%
  mutate(cluster = factor(pcaclust, levels = cluster_order_043, ordered = TRUE)) %>% 
  mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T)) %>% 
  mutate(depthbin = gsub("epi", "Epipelagic", depthbin)) %>% 
  mutate(depthbin = gsub("meso", "Mesopelagic", depthbin))

divdfenv3200 <- left_join(divdf3, env_table, by = "env_sample") %>% dplyr::select(env_sample, numasvs, shan, hill, month, depthbin) %>%
  mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T)) %>% 
  left_join(., grp_df_3200, by = "env_sample") %>%
  mutate(cluster = factor(pcaclust, levels = cluster_order_3200, ordered = TRUE))%>% 
  mutate(depthbin = gsub("epi", "Epipelagic", depthbin)) %>% 
  mutate(depthbin = gsub("meso", "Mesopelagic", depthbin))


divdfenv04 %>% group_by(cluster) %>% summarise_if(is.numeric, min) %>% dplyr::select(cluster,numasvs, shan, hill) %>% mutate(percrich = numasvs/max(divdfenv04$numasvs))
divdfenv04 %>% group_by(cluster) %>% summarise_if(is.numeric, max) %>% dplyr::select(cluster,numasvs, shan, hill) %>% mutate(percrich = numasvs/max(divdfenv04$numasvs))
divdfenv3200 %>% group_by(cluster) %>% summarise_if(is.numeric, min) %>% dplyr::select(cluster,numasvs, shan, hill) %>% mutate(percrich = numasvs/max(divdfenv3200$numasvs))
divdfenv3200 %>% group_by(cluster) %>% summarise_if(is.numeric, max) %>% dplyr::select(cluster,numasvs, shan, hill) %>% mutate(percrich = numasvs/max(divdfenv3200$numasvs))

cols04 <- c("#999999", "#E69F00", "#56B4E9", "purple", "darkblue", "darkgreen")
cols3 <- c("#999999", "#E69F00", "#56B4E9", "purple", "darkblue", "darkgreen")




topm <- 15
pointsize <- 2
richplot043 <- ggplot(data = divdfenv04, aes(x = month, y = numasvs, shape = depthbin, color = cluster))+
  geom_jitter(width = 0.1, size = pointsize)+
  theme_cowplot(font_size = 8)+
  scale_color_manual(values = cols04)+
  theme(plot.margin = unit(c(topm,0,6,0), "pt"))+
  ylab("ASV richness")+
  xlab("Month")+
  labs(color = "Cluster", shape = "Depth zone")+
  ylim(c(0,420))
richplot043

hillplot043 <- ggplot(data = divdfenv04, aes(x = month, y = hill, shape = depthbin, color = cluster))+
  geom_jitter(width = 0.1, size = pointsize)+
  theme_cowplot(font_size = 8)+
  scale_color_manual(values = cols04)+
  theme(plot.margin = unit(c(topm,0,6,0), "pt"))+
  ylab("Hill's evenness")+
  xlab("Month")+
  labs(color = "Cluster", shape = "Depth zone")+
  ylim(c(0,0.5))

hillplot043

shanplot043 <- ggplot(data = divdfenv04, aes(x = month, y = shan, shape = depthbin, color = cluster))+
  geom_jitter(width = 0.1, size = pointsize)+
  theme_cowplot(font_size = 8)+
  scale_color_manual(values = cols04)+
  theme(plot.margin = unit(c(topm,0,6,0), "pt"))+
  ylab("Shannon's diversity")+
  xlab("Month")+
  labs(color = "Cluster", shape = "Depth zone")+
  ylim(c(0,150))

richplot3200 <- ggplot(data = divdfenv3200, aes(x = month, y = numasvs, shape = depthbin, color = cluster))+
  geom_jitter(width = 0.1, size = pointsize)+
  theme_cowplot(font_size = 8)+
  scale_color_manual(values = cols3)+
  theme(plot.margin = unit(c(topm,0,6,0), "pt"))+
  ylab("ASV richness")+
  xlab("Month")+
  labs(color = "Cluster", shape = "Depth zone")+
  ylim(c(0,420))
richplot3200

hillplot3200 <- ggplot(data = divdfenv3200, aes(x = month, y = hill, shape = depthbin, color = cluster))+
  geom_jitter(width = 0.1, size = pointsize)+
  theme_cowplot(font_size = 8)+
  scale_color_manual(values = cols3)+
  theme(plot.margin = unit(c(topm,0,6,0), "pt"))+
  ylab("Hill's evenness")+
  xlab("Month")+
  labs(color = "Cluster", shape = "Depth zone")+
  ylim(c(0,0.5))

shanplot3200 <- ggplot(data = divdfenv3200, aes(x = month, y = shan, shape = depthbin, color = cluster))+
  geom_jitter(width = 0.1, size = pointsize)+
  theme_cowplot(font_size = 8)+
  scale_color_manual(values = cols3)+
  theme(plot.margin = unit(c(topm,0,6,0), "pt"))+
  ylab("Shannon's diversity")+
  xlab("Month")+
  labs(color = "Cluster", shape = "Depth zone")+
  ylim(c(0,150))

divlegnd <- get_legend(hillplot3200+theme(legend.position = "right"))

text043 <- text_grob(expression(paste("0.4-3 ", mu, "m")), rot = 90)
text3200 <- expression(paste("3-200 ", mu, "m"))


ptab <- plot_grid(richplot043 + theme(legend.position = "none"), hillplot043 + theme(legend.position = "none"), shanplot043 + theme(legend.position = "none"), 
                  richplot3200 + theme(legend.position = "none"), hillplot3200 + theme(legend.position = "none"), shanplot3200 + theme(legend.position = "none"),
                  align = "vh", labels = LETTERS[c(1:6)], hjust = -1, nrow = 2, label_size = 12)
ptab

ptab_wleg <- plot_grid(ptab, divlegnd, ncol = 2, rel_widths = c(1,.15))
ptab_wleg

pdf("richplot4_mar24.pdf", width = 7, height = 4)
ptab_wleg
dev.off()





