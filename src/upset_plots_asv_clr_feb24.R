library(UpSetR)
library(readr)
library(here)
library(tidyverse)
library(patchwork)

#### Read ASV tables prepared in script "prep_syndi_asvtable.R" ####
synditab_043_num <- read_delim(here("data", "synditab_043_rn_num.txt"), delim = "\t")
synditab_043_tax <- read_delim(here("data", "synditab_043_rn_tax.txt"), delim = "\t")

synditab_3200_num <- read_delim(here("data", "synditab_3200_rn_num_fromsubsamp.txt"), delim = "\t")
synditab_3200_tax <- read_delim(here("data", "synditab_3200_rn_tax_fromsubsamp.txt"), delim = "\t")

#### Read cluster grouping df, from script "cluster_barplot_figures_clr_new.R" ####
grp_043_df <- read_delim(here("data", "grp_043_dfnew.txt"), delim = "\t")
grp_3200_df <- read_delim(here("data", "grp_3200_dfnew.txt"), delim = "\t")

grp_043_df <- read_delim(here("data", "grp_df_043new.txt"), delim = "\t")
grp_3200_df <- read_delim(here("data", "grp_df_3200new.txt"), delim = "\t")

grp_df_043 <- read_xlsx("data/grp_df_043_pcajan24.xlsx")
grp_df_3200 <- read_xlsx("data/grp_df_3200_pcajan24.xlsx")

cluster_order_043 <- c("JanMar_Atl", "JanMar_Arc", "May_epi", "Aug_epi_N03", "MayAuNo_meso", "1000m")
cluster_order_3200 <- c("JanMar_Atl", "JanMar_Arc",   "May_epi", "Aug_epi", "AuNo_meso", "1000m")

#### from help file ####
# text.scale	
# Numeric, value to scale the text sizes, applies to all axis labels, tick labels, and numbers above bar plot. 
# Can be a universal scale, or a vector containing individual scales in the 
# following format: c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)

text_scale <- c(1.2,1.2,1.2,1.2,1.2,1.2)

#### 0.4-3 ASV level ####
asvtab2_043 <- cbind.data.frame(synditab_043_num, "rsum" = rowSums(synditab_043_num), synditab_043_tax)
asvtab3_043 <- asvtab2_043 %>% arrange(desc(rsum)) %>% filter(rsum>0)
dim(asvtab3_043)
asvtab3_num_043 <- asvtab3_043 %>% purrr::keep(is.numeric) %>% dplyr::select(-rsum)

upset043asv <- data.frame(t(asvtab3_num_043))
upset043asv[upset043asv>0] <- 1
rowSums(upset043asv)
upset043asv$env_sample <- names(asvtab3_043)[c(1:44)]
upset043asv2 <- left_join(upset043asv, grp_df_043, by = "env_sample")
head(upset043asv2[,c(1744:1746)], n = 44)

upset043asv3 <- upset043asv2 %>% group_by(pcaclust) %>% summarise_if(is.numeric, sum)

upset043asv4 <- data.frame(t(data.frame(upset043asv3[,-1])))
upset043asv4[upset043asv4>0] <-1
sum(upset043asv4)
upset043asv4$Identifier <- asvtab3_043$asv_id

names(upset043asv4) <- c(upset043asv3$pcaclust, "Identifier")

upsetrplot_asv_043 <- upset(upset043asv4, sets = cluster_order_043, sets.bar.color = "#56B4E9",
                        order.by = "freq", empty.intersections = "off", text.scale = text_scale)


upsetrplot_asv_043

upsetrplot_data_043 <- upsetrplot_asv_043$New_data
colSums(upsetrplot_asv_043$New_data[,c(1:5)])

#### 3-200 ASV level ####
asvtab2_3200 <- cbind.data.frame(synditab_3200_num, "rsum" = rowSums(synditab_3200_num), synditab_3200_tax)
asvtab3_3200 <- asvtab2_3200 %>% arrange(desc(rsum)) %>% filter(rsum>0)
dim(asvtab3_3200)
asvtab3_num_3200 <- asvtab3_3200 %>% purrr::keep(is.numeric) %>% dplyr::select(-rsum)

upset3200asv <- data.frame(t(asvtab3_num_3200))
upset3200asv[upset3200asv>0] <- 1
rowSums(upset3200asv)
upset3200asv$env_sample <- names(asvtab3_3200)[c(1:44)]
upset3200asv2 <- left_join(upset3200asv, grp_df_3200, by = "env_sample")
head(upset3200asv2[,c(1744:1746)], n = 44)

upset3200asv3 <- upset3200asv2 %>% group_by(pcaclust) %>% summarise_if(is.numeric, sum)
upset3200asv3[,c(1:5)]

upset3200asv4 <- data.frame(t(data.frame(upset3200asv3[,-1])))
upset3200asv4[upset3200asv4>0] <-1
sum(upset3200asv4)
upset3200asv4$Identifier <- asvtab3_3200$asv_id

names(upset3200asv4) <- c(upset3200asv3$pcaclust, "Identifier")

upsetrplot_asv_3200 <- upset(upset3200asv4, sets = cluster_order_3200, sets.bar.color = "#56B4E9",
                        order.by = "freq", empty.intersections = "off", text.scale = text_scale)

upsetrplot_asv_3200
colSums(upsetrplot_asv_3200$New_data[,c(1:5)])

par(mfrow = c(2,1))

ragg::agg_png("upsetr043_mar24.png", width = 8.1, height = 5, units = "in", bg = "white", res = 300, scaling = 1.1)
upsetrplot_asv_043
dev.off()

ragg::agg_png("upsetr3200_mar24.png", width = 8.1, height = 5, units = "in", bg = "white", res = 300, scaling = 1.1)
upsetrplot_asv_3200
dev.off()

plot_grid(upsetrplot_asv_043, upsetrplot_asv_3200, nrow = 2)

par(mfrow = c(2,1))
upsetrplot_asv_043

#upsetrplot_asv_043+upsetrplot_asv_3200 + plot_layout(nrow = 2)
