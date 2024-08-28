### LML Analysis -- Paper 1
## Community changepoint analysis - clean - with upload of simple dataframe

### Libraries --------------
library(lattice)
library(MASS)
library(dplyr)
require(pscl) # alternatively can use package ZIM for zero-inflated 
library(lmtest)
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(ggridges)
library(ecp)
library(gridExtra)
library(ggnewscale)

## Upload data -----------------------

CPUE.w.sec = read.csv("CPUE.w.sec.csv") %>% as.data.frame() %>%
  column_to_rownames(var = "X")
years_included = c(1998:2019)
# Removing rare + stocked taxa ------------ 

## Taxa get removed if they are rare or if they are stocked given 

`%nin%` = Negate(`%in%`) # sets up a way to exclude if in a string

## LML species names 
species_names = c("Brown Bulllhead", "Creek Chub", "Common Shiner","Lake Trout","Central Mudminnow", "Pumpkinseed", "Rainbow Smelt", "Round Whitefish", "Smallmouth Bass", "Slimy Sculpin", "Brook Trout", "White Sucker")

codes = data.frame(species_codes =  c("BB","CC","CS","LT","MM","PS",
                                      "RS","RWF","SMB","SS","ST","WS")) %>% 
  arrange(species_codes)

codes = data.frame(species_names = species_names, 
                   species = codes$species_codes)
species = codes$species
# Colorblind pallete
cbbPalette <- c("#000000",  "#56B4E9", "#D55E00","#009E73","#CEC6C6", "#0072B2","#E69F00","#F0E442",  "#CC79A7")

## FBL 
#species_names = c("Creek Chub","Lake Trout", "Central Mudminnow", "Smallmouth Bass", "Brook Trout", "White Sucker")
#species_names_fall=   c("Creek Chub", "Central Mudminnow", "Smallmouth Bass", "White Sucker")

# CPUE changepoint ------------------------------


v = CPUE.w.sec %>% 
  mutate(y_s = rownames(CPUE.w.sec)) %>%
  pivot_longer(1:length(codes$species),
               names_to = "Species") %>%
  separate(y_s, 
           into = c("Year", "site"), sep = "_") %>%
  unite("ID", 
        c(site:Species), 
        sep = "_", 
        remove = F) %>%
  
  dplyr::select(-site) %>%
  mutate(value = value * 60 * 60 )


#### Graphing Loop --------------------------
graph_list = list() # Create list of graphs for plotting
for(i in 1:length(codes$species)){
  
  # Set up data frame
  x = v %>% filter(Species == codes$species[i]) %>%
    mutate(value = as.numeric(value)) %>%
    dplyr::select(-Species) %>%
    pivot_wider(values_from = value,
                names_from = ID) %>%
    replace(is.na(.), 0) %>%
    dplyr::select(-Year) %>%
    as.matrix()
  
  rownames(x) = (years_included) # Sequence of years 
  
  # Run changepoint analysis 
  output = e.divisive(x, 
                      R = 499, 
                      alpha = 1, 
                      min.size = 2,
                      sig.lvl = .05)
  
  # Format data
  v_mod = left_join(v, (data.frame(Year = as.character(years_included), 
                                   color = output$cluster)))
  
  
  # If there are multiple breakpoints plot the means of the chunks
  
  if((2 %in% v_mod$color) == TRUE){ ## Graphs with changepoints 
    dat_graph = v_mod %>% filter(Species == species[i])%>%
      mutate(value = round(value))
    
    mean_cpue = dat_graph %>%
      group_by(color) %>% 
      dplyr::summarize(mean = mean(value))
    
    #cluster_means = left_join(dat_graph, mean_cpue) # combining the clusters with a mean cpue for that cluster
    
    cluster_line = left_join(dat_graph, mean_cpue) %>%
      select(Year, color, mean) %>% 
      unique()
    
    graph_list[[i]] = ggplot() + 
      geom_jitter(data = dat_graph, aes(x = as.numeric(Year), 
                                        y = value,
                                        col = as.character(color)), 
                  alpha = .2, 
                  width = .2, 
                  size = .9)+
      scale_colour_manual(values=cbbPalette) +
      ylab(paste(species_names[i], "Indv / Hour") )+ xlab("Year") + 
      theme(legend.position="none") + 
      geom_vline(xintercept = 2000, linetype = "dashed") +
      theme_minimal()+ 
      theme(legend.position = "none") +
      new_scale_color() + 
      scale_colour_manual(values = c("#707173",
                                     "#7088b8",
                                     "#E69F00", 
                                     "#6fa373",
                                     "#E69F00"))+
      geom_line(data = cluster_line, 
                aes(x = as.numeric(Year),
                    y = mean,
                    col = as.character(color),),
                lwd = 1) + 
      theme(text = element_text(size = 7)) + 
      theme(axis.text.x = element_text(angle =90, size = 9)) + 
      theme(axis.text.y = element_text(size = 9)) +
      xlim(1997, 2019) + 
      theme(plot.margin = margin(.25, .5, .25, .5, "cm")) 
    
  } else { ## Graphs with regressions
    dat_graph = v_mod %>% filter(Species == species[i])%>%
      mutate(value = round(value)) 
    
    dat_graph.yr = dat_graph %>%
      filter(Year > 2000)
    
    
    dat_graph.pred = dat_graph %>%
      filter(Year > 2000)   %>% 
      mutate(Year = as.numeric(Year)) %>%
      mutate(Year = scale(Year)[,1]) 
    
    
    mean_cpue = dat_graph %>%
      group_by(color) %>% 
      dplyr::summarize(mean = mean(value))
    
    fred = left_join(dat_graph, mean_cpue)
    
    try(M4 <- zeroinfl(value ~ (Year) | (Year),
                       dist = 'negbin',
                       data = dat_graph.pred))
    
    
    Pred<- predict(M4,newdata = dat_graph.pred, type = "response")
    
    pred_data = cbind(dat_graph.yr, Pred)
    
    graph_list[[i]] = ggplot() + 
      theme_minimal()+
      geom_jitter(data = dat_graph,
                  aes(x = as.numeric(Year),y = value),
                  color = "black", 
                  alpha = .2,
                  width = .2,
                  size = .9)+
      geom_line(data = pred_data,
                aes(x = as.numeric(Year),y = Pred),
                col = "#726F6F",
                lwd = 1) + 
      geom_vline(xintercept = 2000,
                 linetype = "dashed") +
      ylab(paste(species_names[i], "Indv / Hour") ) +
      xlab("Year") + 
      scale_colour_manual(values=cbbPalette) + 
      theme(legend.position="none") +
      theme(text = element_text(size = 7)) + 
      theme(axis.text.x = element_text(angle =90, size = 9)) + 
      theme(axis.text.y = element_text(size = 9)) +
      theme(plot.margin = margin(.25, .5, .25, .5, "cm")) +
      xlim(1997, 2019) 
    
  }
}
# Grid Arrange Graphic 
do.call("grid.arrange", c(graph_list, ncol=4))
