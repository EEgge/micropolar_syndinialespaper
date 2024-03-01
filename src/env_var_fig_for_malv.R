require(readr)
require(here)
require(tidyverse)
library(oce)

env_table0 <- read_delim(here("data", "env_data_depths_new.txt"), delim = "\t")
env_profiles0 <- read_delim(here("data", "env_data_profiles.txt"), delim = "\t")

axte <- 0.8
axti <- 1
annsize <- 8
lts <- tibble("station" = levels(as.factor(env_table$station)), "lityp" = c(1,2,1,2,23,3,4,1,2,23,1,2,23,1,2,23))

env_profiles1 <- env_profiles0 %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_profiles2 <- env_profiles1 %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a, DOC_uM, POC_uM, HNF, Large_HNF)  %>%  
  filter(depth_m <= 1000) %>% 
  mutate(sigma_oce = swSigmaTheta(salinity = CTD.S,
                                  temperature = CTD.T,
                                  pressure = depth_m))

env_profiles3 <- left_join(env_profiles2, lts, by ="station") %>% mutate(lityp = as.character(lityp))
env_profiles_melt <- reshape2::melt(env_profiles3, measure.vars = c("CTD.S", "CTD.T", "sigma","Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a", "DOC_uM", "POC_uM", "HNF", "Large_HNF", "sigma_oce")) #"variable" is now the type of parameter


env_points1 <- env_table0 %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_points2 <- env_points1 %>% select(month, station, depth_m, depthbin, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a, DOC_uM, POC_uM, HNF, Large_HNF)  %>% 
  mutate(depth_m = as.numeric(as.character(depth_m))) %>%  
  filter(station %in% unique(lts$station) & depth_m<= 1000)%>% 
  mutate(sigma_oce = swSigmaTheta(salinity = CTD.S,
                                  temperature = CTD.T,
                                  pressure = depth_m))


env_points3 <- left_join(env_points2, lts, by ="station") %>% mutate(lityp = as.character(lityp))
env_points_melt <- reshape2::melt(env_points3, measure.vars = c("CTD.S", "CTD.T", "sigma", "Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a", "DOC_uM", "POC_uM", "HNF", "Large_HNF", "sigma_oce")) %>% 
  mutate(depthfact = cut(depth_m, breaks = c(0,100,200,500,1001), include.lowest = T, right = F))


#"variable" is now the type of parameter
#env_points_melt <- pivot_longer(env_points3, cols = c("CTD.S", "CTD.T", "sigma", "Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a", "DOC_uM", "POC_uM")) #"variable" is now the type of parameter


envplot <- function(envvar_profile, envvar_point, ylab) {
  ggplot(data = envvar_profile)+ 
    aes(x = depth_m, y = value, lty = lityp, col =lityp, group = lityp)+
    geom_point(data = envvar_point, aes(x = depth_m, y = value, shape = lityp, col =lityp))+
    geom_line()+
    facet_grid(.~ month)+
    theme_classic()+
    theme(strip.background = element_rect(colour="white", fill="white"))+
    #scale_x_reverse()+
    scale_x_sqrt(labels = c("0", "20", "500", "1000"), breaks = c(0, 20, 500, 1000))+
    scale_linetype_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
    scale_color_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
    scale_shape_manual(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"), values = c(19,17,3,4,7))+
    theme(axis.text = element_text(size = rel(axte)))+
    theme(axis.title = element_text(size = rel(axti)),
          strip.text = element_text(size = 11))+
    xlab("Depth")+
    ylab(ylab)+
    coord_flip()+
    scale_x_sqrt(labels = c("0", "20", "500", "1000"), breaks = c(0, 20, 500, 1000))+
    scale_x_reverse()+
        theme(legend.position = "right",
          legend.text = element_text(size = 11))+
    #guides(shape = guide_legend(override.aes = list(size = 3)))+
    theme(legend.key.size = unit(2,"line"))
}


sal <- envplot(env_profiles_melt %>% filter(variable == "CTD.S"), env_points_melt%>% filter(variable == "CTD.S"), "Salinity [PSU]")
sal
temp <- envplot(env_profiles_melt %>% filter(variable == "CTD.T"), env_points_melt%>% filter(variable == "CTD.T"), expression(paste("Temperature [", degree,"C]")))
temp
sigma_oce <- envplot(env_profiles_melt %>% filter(variable == "sigma_oce"), env_points_melt%>% filter(variable == "sigma_oce"), expression(paste("Density ", sigma[theta], " [kg ", m^-3,"]")))
sigma_oce

sigma <- envplot(env_profiles_melt %>% 
                   filter(variable == "sigma"), env_points_melt%>% filter(variable == "sigma"), expression(paste("Density ", sigma)))
sigma

ragg::agg_png(filename = "sigma_oce_malv.png", width = 1400, height = 550, res = 300, scaling = 0.6)
sigma_oce
dev.off()

ragg::agg_png(filename = "sal_malv.png", width = 1400, height = 550, res = 300, scaling = 0.6)
sal
dev.off()

inorgN <- envplot(env_profiles_melt %>% filter(variable == "Total_inorgN_uM"), env_points_melt%>% filter(variable == "Total_inorgN_uM"), expression(paste("Total inorganic N", " [", mu,"M]")))
phos <- envplot(env_profiles_melt %>% filter(variable == "PO4_uM"), env_points_melt%>% filter(variable == "PO4_uM"), expression(paste("PO"[4]^"3-", " [", mu,"M]")))
sil <- envplot(env_profiles_melt %>% filter(variable == "SiOH4_uM"), env_points_melt%>% filter(variable == "SiOH4_uM"), expression(paste("SiOH"[4], " [", mu,"M]")))
chl <- envplot(env_profiles_melt %>% filter(variable == "Chl_a"), env_points_melt%>% filter(variable == "Chl_a"), expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")))

envfig <- ggarrange(sal+theme(legend.title =element_text(size = 10),
                              legend.text = element_text(size =9), axis.title.y = element_blank()), temp+theme(axis.title.y = element_blank()), sigma+theme(axis.title.y = element_blank()), 
                    inorgN+theme(axis.title.y = element_blank()), phos+theme(axis.title.y = element_blank()), sil+theme(axis.title.y = element_blank()), chl+theme(axis.title.y = element_blank()), nrow = 7, common.legend = TRUE, legend = "right",
                    labels = LETTERS[1:7])
envfig


poc <- envplot(env_profiles_melt %>% filter(variable == "DOC_uM"), env_points_melt%>% filter(variable == "DOC_uM"), expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")))
poc

pchsize <- 2.5
chlplot <- env_points_melt%>% filter(variable == "Chl_a") %>% 
  ggplot(., aes(x = station, y = value, color = depthfact, shape = depthbin))+
  geom_jitter(width = 0.1, size =pchsize)+
  scale_y_sqrt()+
  scale_color_viridis_d()+
  facet_grid(.~month, space = "free", scales = "free")+
  labs(y = expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")), x = "Month", color = "Depth", shape = "Zone")

docplot <- env_points_melt%>% filter(variable == "DOC_uM") %>% 
  ggplot(., aes(x = station, y = value, color = depthfact, shape = depthbin))+
  geom_jitter(width = 0.1, size =pchsize)+
  scale_color_viridis_d()+
  facet_grid(.~month, space = "free", scales = "free")+
  labs(y = expression(paste("DOC "," [", mu,"M]")), x = "Month", color = "Depth", shape = "Zone")

docplot

pocplot <- env_points_melt%>% filter(variable == "POC_uM") %>% 
  ggplot(., aes(x = station, y = value, color = depthfact, shape = depthbin))+
  geom_jitter(width = 0.1, size =pchsize)+
  scale_y_sqrt()+
  scale_color_viridis_d()+
  facet_grid(.~month, space = "free", scales = "free")+
  labs(y = expression(paste("POC "," [", mu,"M]")), x = "Month", color = "Depth", shape = "Zone")

pocplot

hnfplot <- env_points_melt%>% filter(variable == "HNF") %>% 
  ggplot(., aes(x = station, y = value, color = depthfact, shape = depthbin))+
  geom_jitter(width = 0.1, size =pchsize)+
  scale_y_sqrt()+
  scale_color_viridis_d()+
  facet_grid(.~month, space = "free", scales = "free")+
  labs(y = expression(paste("HNF, cells ", L^-1)), x = "Month", color = "Depth", shape = "Zone")

hnfplot

largehnfplot <- env_points_melt%>% filter(variable == "Large_HNF") %>% 
  ggplot(., aes(x = station, y = value, color = depthfact, shape = depthbin))+
  geom_jitter(width = 0.1, size =pchsize)+
  scale_y_sqrt()+
  scale_color_viridis_d()+
  facet_grid(.~month, space = "free", scales = "free")+
  labs(y = expression(paste("Large_HNF, cells ", L^-1)), x = "Month", color = "Depth", shape = "Zone")

largehnfplot

require(patchwork)

envplotmalv <- chlplot+docplot+pocplot+hnfplot+largehnfplot +plot_layout(nrow = 5, guides = "collect")


ragg::agg_png(filename = "envplotmalv.png", width = 1400, height = 2000, res = 300, scaling = 0.6)
envplotmalv
dev.off()
