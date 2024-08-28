
### LML Analysis -- Paper 1


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
`%nin%` = Negate(`%in%`) # sets up a way to exclude if in a string

## Functions source -----------
# This is just setup at the project working directory. Use option in upper right corner of R to get into project directory. For example, on my computer,its stored in my family one-drive
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML_SMB_removal")

### Note - this source file will upload the data files. But you need to make sure to correct the source location for those data files for this to work...
source("Function_Source_Files/AFRP_Functions.R")

### Data -------------------

#sample = read.csv("MA2276_Code//Data/FISH_SAMPLE_edited.csv")
sample = read.csv("../AFRP/MA2276_Code/Data/FISH_SAMPLE_2022_editedsitenumbers.csv")



BEF_data_unfiltered =left_join(fish, sample, by = "YSAMP_N") %>% 
  left_join(sites, by = "SITE_N") %>% 
  left_join(shoreline_length, by = "SITE_N") %>%
  separate(SITE_N,  into = c("GEAR", "WATER","SITE")) %>%
  filter(WATER == "LML" & GEAR == "BEF"& GEAR_CODE == "NAF" & YEAR < 2020 & MONTH %in% c(5,6)) %>%
  mutate(SITE = as.character(SITE))

## Fixing site issues 


###bef_unfiltered_saved = BEF_data_unfiltered
#BEF_data_unfiltered = bef_unfiltered_saved

site_bin = c(1,4,5,9,13,16,18, 21, 23,27,31) ## for the simulations

post_2002 = BEF_data_unfiltered %>% filter(YEAR == 2005) %>% select(SITE) %>% unique() %>% mutate(SITE_num = parse_number(SITE)) %>% 
  mutate(site_bin = .bincode(SITE_num, site_bin)) ## Create matrices for each cluster of years that were sampled in one way or another


BEF_1998 = BEF_data_unfiltered %>%
  filter(YEAR == 1998) %>% 
  select(SITE) %>% unique() %>% 
  mutate(site_bin = c(1,1,1,1,2,3,3,3,4,4,5,7,8,11,10,10,4,5)) %>% 
  mutate(SITE_num = c(1:18)) 
  
  
bef_1998 = BEF_data_unfiltered %>% filter(YEAR == 1998) %>% 
  left_join(BEF_1998)


bef_1999 = BEF_data_unfiltered %>%
  filter(YEAR == 1999) %>% 
  dplyr::select(SITE) %>% unique() %>% 
  filter(SITE !="NA") %>% 
  mutate(site_bin = c(2,3,4,1,5,7,8, 10)) %>% 
  mutate(SITE_num = c(1:8))

BEF_1999 = BEF_data_unfiltered %>% filter(YEAR == 1999) %>% 
  left_join(bef_1999)

site_matrix = rbind(bef_1999, BEF_1998, post_2002) %>%
  as.data.frame()


BEF_data_2002 = BEF_data_unfiltered %>%
  filter(YEAR > 1999) %>%
  left_join(post_2002)

BEF_data_unfiltered = left_join(BEF_data_unfiltered, site_matrix) %>% ## Make sure to use this BEF_data_unfiltered for final graphs
  mutate(SITE_cat = case_when(YEAR %nin% c(1998, 1999) ~ parse_number(SITE),
                              YEAR %in% c(1998, 1999) ~ as.numeric(SITE_num))) %>%
  dplyr::select(-SITE) %>% 
  rename(SITE = SITE_cat) %>%
  filter(SITE != "NA") 






BEF_data_unfiltered = rbind(bef_1998, BEF_1999, BEF_data_2002) %>%
  mutate(SITE_cat = case_when(YEAR %nin% c(1998, 1999) ~ parse_number(SITE),
                              YEAR %in% c(1998, 1999) ~ as.numeric(SITE_num))) %>%
  dplyr::select(-SITE) %>% 
  rename(SITE = SITE_cat) %>%
  filter(SITE != "NA")



# Removing rare + stocked taxa ------------ 

## Taxa get removed if they are rare or if they are stocked given 



rare_threashold = 50 ## change this based on preference. Here it is filtering out rare fish (NRD and BND)

rare = BEF_data_unfiltered %>% ## This defines rare species 
  group_by(SPECIES) %>% 
  summarise(frequency = n()) %>% 
  filter(frequency < rare_threashold)

stocked = c("LLS", "RT", "ST") ## Stocked fish in LML to be excluded from analysis


BEF_data = BEF_data_unfiltered%>%
  filter(SPECIES %nin% c(stocked, rare$SPECIES)) %>% 
  filter(YEAR < 2020) %>% 
  filter(SPECIES != "SMB" | YEAR != 2000 | DAY_N < 160) %>% 
  filter(YEAR != 2002) #%>% 
  #filter(YEAR != 1999)## Filter out BEF SMB data from the year 2000 that's later than DAY_N 160. Change this around depending on how you want to filter 2000... 
  

BEF_data %>%
  filter(SPECIES == "SMB") %>%
  select(SPECIES, YEAR, DAY_N) %>% 
  unique()



## LML 
species_names = c("brown bullhead", "creek chub", "common shiner","lake trout","central mudminnow", "pumpkinseed", "rainbow smelt", "round whitefish", "smallmouth bass", "slimy sculpin", "white sucker")

codes = data.frame((species_codes = unique(BEF_data$SPECIES))) %>% 
  arrange(species_codes)

codes = data.frame(species_names = species_names, 
                   species = codes$X.species_codes)
species = codes$species
# Colorblind pallete
cbbPalette <- c("#000000",  "#56B4E9", "#D55E00","#009E73","#CEC6C6", "#0072B2","#E69F00","#F0E442",  "#CC79A7")

## FBL 
#species_names = c("Creek Chub","Lake Trout", "Central Mudminnow", "Smallmouth Bass", "Brook Trout", "White Sucker")
#species_names_fall=   c("Creek Chub", "Central Mudminnow", "Smallmouth Bass", "White Sucker")

# CPUE changepoint ------------------------------

#### Data setup ---------------------
#CPUE.w.sec = rbind(CPUE.w.sec_all, bef_1999)

CPUE.w.sec = ((CPUE_wide_seconds(BEF_data) %>%
                 unite("Group", c(YEAR, SITE)) %>% 
                 column_to_rownames(., var = "Group") %>% 
                 mutate(sumrow = rowSums(.)) %>%
                 filter(sumrow>0) %>%
                 dplyr::select(-sumrow)))

write.csv(CPUE.w.sec, file="CPUE.w.sec.csv")
#dev.off()

CPUE.w.sec.a = ((CPUE_wide_seconds_avg(BEF_data) %>% 
                   column_to_rownames(., var = "YEAR")))

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

# abundance for bass before vs. after


mean(CPUE.w.sec.a$SMB[4:21]) / mean(CPUE.w.sec.a$SMB[1:3])

#### Graphing Loop --------------------------
graph_list = list() # Create list of graphs for plotting
for(i in 1:length(codes$species)){
  
  # Set up data frame
  x = v %>% filter(Species == codes$species[i]) %>%
    mutate(value = as.numeric(value)) %>%
    dplyr::select(-Species) %>%
    pivot_wider(values_from = value,
                names_from = ID) %>% 
    arrange(Year) %>% 
    mutate(Year = as.numeric(Year)) %>%

    replace(is.na(.), 0) %>%
    dplyr::select(-Year) %>%
    as.matrix()
  
  x[1,19:32] = "NA" #add back in when not doing simulation
  x[2,9:32] = "NA"
  
  rownames(x) = c(1998:2019)[c(-5)] # Sequence of years 
  
  # Run changepoint analysis 
  output = e.divisive(x, 
                      R = 1000, 
                      alpha = 1, 
                      min.size = 2,
                      sig.lvl = .05)
  print(output$p.values)
  
  # Format data
  v_mod = left_join(v, (data.frame(Year = rownames(CPUE.w.sec.a), 
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
      theme(axis.text.x = element_text(angle =90, size = 9, vjust=.6)) + 
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
      
      
      
      
      pred_data = 0
      if( max(summary(M4)[1]$coefficients$count[1:2,4]) < .05){
        
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
          theme(axis.text.x = element_text(angle =90, size = 9, vjust = .60)) + 
          theme(axis.text.y = element_text(size = 9)) +
          theme(plot.margin = margin(.25, .5, .25, .5, "cm")) +
          xlim(1997, 2019) 
        
        
        
      } else {
        
        graph_list[[i]] = ggplot() + 
          theme_minimal()+
          geom_jitter(data = dat_graph,
                      aes(x = as.numeric(Year),y = value),
                      color = "black", 
                      alpha = .2,
                      width = .2,
                      size = .9)+
          
          geom_vline(xintercept = 2000,
                     linetype = "dashed") +
          ylab(paste(species_names[i], "Indv / Hour") ) +
          xlab("Year") + 
          scale_colour_manual(values=cbbPalette) + 
          theme(legend.position="none") +
          theme(text = element_text(size = 7)) + 
          theme(axis.text.x = element_text(angle =90, size = 9, hjust = .6)) + 
          theme(axis.text.y = element_text(size = 9)) +
          theme(plot.margin = margin(.25, .5, .25, .5, "cm")) +
          xlim(1997, 2019) 
      }
        
        
      
      

      
      print(paste(species[i]))
      print(summary(M4))
    }
}

# Grid Arrange Graphic 
do.call("grid.arrange", c(graph_list, ncol=4))


## CPUE Changepoint Summary

cp_data = read.csv("../AFRP/MA2276_Code/Data/LML_CP_data.csv") %>% left_join(codes)
species_data = read.csv("../AFRP/MA2276_Code/Data/LML_SPECIES_DATA.csv")  %>%
  arrange(MEAN_LML_LENGTH) %>% 
  rename(species = SPECIES) 

cp_data = left_join(cp_data, species_data) %>% arrange(MAX_LML_LENGTH)

cp_data %>% ggplot(aes(x = year, 
                       y = species_names,
                       label = direction_shape,
                       color = direction_shape)) +
  geom_text(size = 10, key_glyph = "rect") +
  xlim(2000, 2020) + 
  geom_segment(y = "lake trout", yend = "lake trout", x = 2001, xend = 2019,
             col = "#00BFC4", 
             lwd = 1.5) + 
  #geom_segment(y = "White Sucker", yend = "White Sucker", x = 2001, xend = 2020,
             #col = "#00BFC4", 
             #lwd = 1.5) +
  geom_segment(y= "creek chub", yend = "creek chub", x = 2001, xend = 2019,
             col = "#F8766D",
             lwd = 1.5) + 
  geom_segment(y= "round whitefish", yend = "round whitefish", x = 2001, xend = 2019,
               col = "#F8766D",
               lwd = 1.5) +
  geom_segment(y= "common shiner", yend = "common shiner", x = 2001, xend = 2019,
               col = "#F8766D",
               lwd = 1.5) +
  labs(color = "Relative Abundance") +
  labs(color = "Relative Abundance") + 
  scale_color_manual(labels = c("Decrease", "Increase"), 
                     values = c("#F8766D", "#00BFC4")) + 
  ylab("") + 
  xlab("Year") +
  scale_y_discrete(limits = unique(cp_data$species_names)) + 
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=9)) + 
  theme_minimal() + 
  geom_vline(xintercept = 2000, linetype = 2)
## Write a PDF 
pdf(file = paste("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/MA2276_Code/Graphics/LMLP1/max_length_CPUE_CP.pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4)
  
dev.off()
   
### Length Changepoint


BEF_data = BEF_data %>% filter(YEAR < 2020) %>% filter(YEAR != 2002)
# Splitting up the change point for length editing it 

graph_list_length = list()
for(i in 1:length(species)){
  
  #pdf(file = paste("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/MA2276_Code/Graphics/LMLP1/",species[i], "_length.pdf", sep = ""),   # The directory you want to save the file in
     # width = 4, # The width of the plot in inches
     # height = 4) # The height of the plot in inches
  year_min = 1997
  x = BEF_data %>% filter(SPECIES == species[i]) %>% 
    filter(YEAR >= year_min) %>%
    select(LENGTH, YEAR, SITE) %>% na.omit() %>% 
    group_by(YEAR, SITE) %>%
    summarize(median_L = mean(LENGTH,na.rm = T)) %>%
    ungroup() %>%
    filter(SITE != "") %>%
    pivot_wider(names_from = SITE, values_from = median_L) 
  
  output = e.divisive(x, R= 10000,
                      alpha = 1,
                      min.size = 2, 
                      sig.lvl  = .05)
  print(output$p.values)
  output_dat = data_frame(YEAR = unique(x$YEAR), 
                          cluster = output$cluster)
  
  dat = BEF_data %>%
    filter(YEAR >= year_min) %>%
    filter(SPECIES == species[i]) %>% 
    left_join(output_dat)
  
  lm_obj = summary(lm(dat$LENGTH ~ dat$YEAR, 
                      na.rm = T))
  
  
  if((2 %in% output_dat$cluster) == TRUE){
    
    
    mean_length = dat %>%
      group_by(cluster) %>% 
      summarize(mean = mean(LENGTH, na.rm = T))
    
    fred = left_join(dat, mean_length)
    
    

    

      
    cluster_means = left_join(fred, mean_length)
      

    
    
    graph_list_length[[i]] = ggplot() +
      theme_minimal() +
      geom_jitter(data = dat, aes(x = as.numeric(YEAR), 
                                  y = LENGTH, 
                                  col = as.character(cluster),
                                  alpha = .2,
                                  width = .2,
                                  ), size = .9) +
      labs(col = paste(species[i])) + 
      ylab("Length") +
      ylab(paste(species_names[i], "TL (mm)")) + 
      theme(text = element_text(size = 9)) + 
      theme(legend.position="none") +
      xlab("Year") +
      xlim(1997, 2021) + 
      scale_colour_manual(values=cbbPalette) + 
      geom_vline(xintercept = 2000, linetype = "dashed") + 
      new_scale_color() + 
      scale_colour_manual(values = c("#707173", 
                                              "#7088b8",
                                              "#E69F00", 
                                              "#6fa373"))+
      geom_line(data = cluster_means, aes(x = YEAR, 
                                          y = mean,
                                          col = as.character(cluster)),
                    size = 1) +
      theme(axis.text.x = element_text(angle = 90))
    
    
    
    
    
  } else {
    
    data_graph = dat %>%
      filter(SPECIES == species[i])
    
    summary_lm = summary(lm(data = data_graph, LENGTH ~ YEAR))[[4]] %>% as.matrix() 
    
    
    if(max(summary_lm[,4]) < .05){
      
      graph_list_length[[i]] = ggplot(data = data_graph, 
                                    aes(x = YEAR, y = LENGTH)) +
      theme_minimal() +
      geom_jitter(aes(col = 'red',
                      alpha = .2),
                  width = .2, 
                  size = .9) +
      labs(col = paste(species[i])) + 
      ylab("Length") +
      geom_smooth(method = "lm", color = "#726F6F", size = 1, se = F) + 
      ylab(paste(species_names[i], "TL (mm)")) + 
      theme(text = element_text(size = 9)) + 
      theme(legend.position="none",
            #axis.title.x = element_blank(),
      ) + 
      xlim(1997, 2021) + 
      xlab("Year") +
      scale_colour_manual(values=cbbPalette)+ 
      geom_vline(xintercept = 2000, linetype = 2) +
      theme(axis.text.x = element_text(angle = 90, hjust = -10)) +
      ylim(0,800)
      
      print(summary(lm(data = data_graph, LENGTH ~ YEAR)))
      
    } else {
      
      graph_list_length[[i]] = ggplot(data = data_graph, 
                                      aes(x = YEAR, y = LENGTH)) +
        theme_minimal() +
        geom_jitter(aes(col = 'red',
                        alpha = .2),
                    width = .2, 
                    size = .9) +
        labs(col = paste(species[i])) + 
        ylab("Length") +
        ylab(paste(species_names[i], "TL (mm)")) + 
        theme(text = element_text(size = 9)) + 
        theme(legend.position="none",
              #axis.title.x = element_blank(),
        ) + 
        xlim(1997, 2021) + 
        xlab("Year") +
        scale_colour_manual(values=cbbPalette)+ 
        geom_vline(xintercept = 2000, linetype = 2) +
        theme(axis.text.x = element_text(angle = 90))
      
    }
  
  
  
  
  }
  
  #print(h)
  
  #dev.off()
  
}
#install.packages("ggnewscale")
do.call("grid.arrange", c(graph_list_length, ncol=4))



species_length_changepoint = read.csv("../AFRP/MA2276_Code/Data/changepoint_length.csv") %>% left_join(codes)

cp_data = left_join(species_length_changepoint, species_data) %>% arrange(MAX_LML_LENGTH)


# In that data file i added in - signs next to species that dont actually have a trend to get rid of an extra color on the graphs... seemms to have worked but is not elegant 

species_length_changepoint %>% 
  
  ggplot(aes(x = year, 
                                          y = species_names,
                                          label = direction,
                                          color = direction)) + 
  geom_text(size = 10,key_glyph = "rect") +
  geom_vline(xintercept = 2000, linetype = 2) +
  xlim(2000, 2020) + 
  theme_minimal() +
  #theme(legend.position = 'none') + 
  geom_segment(y = "common shiner", yend = "common shiner", x = 2001, xend = 2019,
             col = "#00BFC4", 
             lwd = 1.5) + 
  
  geom_segment(x = 2001, xend = 2019, 
               y = "brown bullhead", yend = "brown bullhead",
               lwd = 1.5, col = "#F8766D") +


  geom_segment(y = "white sucker", yend = "white sucker", x = 2001, xend = 2019,
             col = "#00BFC4", 
             lwd = 1.5) +
  labs(color = "Total Length (mm)", size = 10) + 
  scale_color_manual(labels = c("Decrease", "Increase"), values = c("#F8766D", "#00BFC4")) + 
  ylab("") + 
  xlab("Year") +
  scale_y_discrete(limits = unique(cp_data$species_names))+
  theme(legend.title=element_text(size=10),
                            legend.text=element_text(size=9))



pdf(file = paste("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/MA2276_Code/Graphics/LMLP1/max_length_LENGTH_CP.pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4)

dev.off()

##------------------------------


cat = BEF_data %>% filter(SPECIES == "WS") %>% select(LENGTH, YEAR) 
summary(lm(cat$LENGTH ~ cat$YEAR))



### ------------box plots final ----------------------------

CPUE.w.sec %>% rownames_to_column(var = "year_site") %>%
  separate(year_site, into = c("year", "site")) %>% 
  group_by(year) %>%
  filter(year %in% c(2000,2001, 2019)) %>% 
  pivot_longer(3:13, values_to = "CPUE", names_to = "species") %>% 
  left_join(codes) %>%
  mutate(CPUE = CPUE * 60 * 60) %>%
  ggplot(aes(x = year, y = CPUE)) +
  theme_bw() +
  geom_boxplot()+ 
  #scale_y_continuous(
   # labels = scales::number_format(accuracy = 0.01)) +
  xlab("") + ylab("CPUE (Individuals / Hour)") +
  facet_wrap(~species_names, scales = "free_y")


BEF_data %>%
  rename(species = SPECIES) %>%
  left_join(codes) %>%
  group_by(species) %>% filter(YEAR %in% c(2000, 2001, 2019)) %>%
  
  ggplot(aes(x = as.character(YEAR), y = LENGTH)) +
  theme_bw() + geom_boxplot() + 
  facet_wrap(~species_names, scales = "free_y") + xlab("") + ylab("Total Length (mm)")



### ------------ Wilcoxon Tests ----------------------------

c.data = v %>% select(Year, Species,ID,  value) %>% rename(year = Year) %>% separate(ID, into = c("site", "sp.del")) %>% select(-sp.del) %>% rename(species = Species) %>% rename(CPUE = value)



Build the means and standard deviations into this table... ## Text so i don't foorget to do this
wilcox_list = list()

for(i in 1:length(species)){
  
    dat = c.data %>% 
        filter(species == species[i])
    
    # data frames
    data.2000 = dat %>% filter(year == 2000)
    data.2001 = dat %>% filter(year == 2001)
    data.2019 = dat %>% filter(year == 2019)
  
    # Wilcox tests 
    a = wilcox.test((data.2000 %>% filter(species==species[i]))$CPUE, (data.2001 %>% filter(species==species[i]))$CPUE, exact = F)
    b = wilcox.test((data.2000 %>% filter(species==species[i]))$CPUE, (data.2019 %>% filter(species==species[i]))$CPUE, exact = F)
    c = wilcox.test((data.2001 %>% filter(species==species[i]))$CPUE, (data.2019 %>% filter(species==species[i]))$CPUE, exact = F)
    
    
    # Build dataframe 
    
    wilcox_list[[i]] = data.frame(name = rep(species[i], 3),
                                  year = c(2000, 2001, 2019),
                                  mean = c(mean(data.2000$CPUE), mean(data.2001$CPUE),mean(data.2019$CPUE)),
                                  p.value = c(a$p.value, b$p.value, c$p.value),
                                  W_statistic = c(a$statistic, b$statistic, c$statistic))
    
    
    
    
}


wilcox_list
wilcox_table = rbind(wilcox_list[[1]], 
      wilcox_list[[2]],
      wilcox_list[[3]],
      wilcox_list[[4]],
      wilcox_list[[5]],
      wilcox_list[[6]],
      wilcox_list[[7]],
      wilcox_list[[8]],
      wilcox_list[[9]],
      wilcox_list[[10]],
      wilcox_list[[11]]
      ) %>% as.data.frame() %>% mutate(p.value = round(p.value, digits = 4)) %>% 
  mutate(p.value = replace(p.value, p.value < 0.0001, "<0.0001")) 
wilcox_table
write.csv(wilcox_table, file = "wilcox_table_CPUE.csv")

## ------------------------- Wilcox Length ---------------------



wilcox_list_length = list()
species_no.mm = species[c(-1)]
species_no.mm = species
for(i in 1:length(species_no.mm)){
  dat = BEF_data %>% filter(SPECIES == species_no.mm[i]) 
  
  data.2000 = dat %>% filter(YEAR == 2000)
  
  data.2001 = dat %>% filter(YEAR == 2001)
  data.2019 = dat %>% filter(YEAR == 2019)
  
  if(length(data.2000[,1]) > 4){
    a = wilcox.test(na.omit(data.2000$LENGTH), na.omit(data.2001$LENGTH), exact = F)
  } else {
    a = data.frame(p.value = "NA", statistic = "NA")
  }
  
  if(length(data.2000[,1]) > 4 & length(data.2019[,1]) > 4){
    b = wilcox.test(data.2000$LENGTH, data.2019$LENGTH, exact = F)
  } else {
    b = data.frame(p.value = "NA", statistic = "NA")
  }
  
  if(length(data.2019[,1]) > 4 & length(data.2001[,1]) > 4 ){
    c = wilcox.test(data.2001$LENGTH, data.2019$LENGTH, exact = F)
  } else {
    c = data.frame(p.value = "NA", statistic = "NA")
  }
  
  wilcox_list_length[[i]] = data.frame(name = rep(species_no.mm[i], 3),
                                       year = c(2000, 2001, 2019), 
                                       length = c(mean(data.2000$LENGTH),mean(data.2001$LENGTH),mean(data.2019$LENGTH) ),
                                       p.value = c(a$p.value, b$p.value, c$p.value), 
                                       W_statistic = c(a$statistic, b$statistic, c$statistic))

}

wilcox_list_length

wilcox_table = rbind(wilcox_list_length[[1]], 
                     wilcox_list_length[[2]],
                     wilcox_list_length[[3]],
                     wilcox_list_length[[4]],
                     wilcox_list_length[[5]],
                     wilcox_list_length[[6]],
                     wilcox_list_length[[7]],
                     wilcox_list_length[[8]],
                     wilcox_list_length[[9]],
                     wilcox_list_length[[10]],
                     wilcox_list_length[[11]]) 
  


write.csv(wilcox_table, file = "wilcox_table_length.csv")


### end ------------------------------------------------





BEF_data %>% group_by(SPECIES, YEAR) %>% summarize(max = max(LENGTH, na.rm=T)) %>%
  ggplot(aes(x = YEAR, y = max, col = SPECIES)) + geom_point() + geom_smooth(method= "lm", se = F)
  

BEF_data %>% group_by(SPECIES, YEAR) %>% summarize(min = min(LENGTH, na.rm=T)) %>%
  ggplot(aes(x = YEAR, y = min, col = SPECIES)) + geom_point() + geom_smooth(method= "lm", se = F)


BEF_data %>% group_by(SPECIES) %>% summarize(max = max(LENGTH, na.rm= T))
BEF_data %>% group_by(SPECIES) %>% summarize(min = min(LENGTH, na.rm= T))











v %>% separate(ID, into = c("site", "Species")) %>% filter(Year == 2021) %>% mutate(site= as.numeric(site)) %>% filter(site %in% c(14,15,16,18,19,20,21,22)) %>%
  mutate(site = as.factor(site)) %>%
  ggplot(aes(x = site, y = value, fill = Species)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +ylab("CPUE ind/hour")



