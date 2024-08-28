
## Please load in data associated with LML_P1_community-response.R to get the BEF_data and v data frame

## Sites sampled in May vs. June in 2012
month_bin = c(0,5,6)

may_sites = BEF_data %>% 
  filter(YEAR == 2012) %>% 
  select(SITE, MONTH) %>% 
  unique() %>% 
  mutate(value = 1) %>% 
  pivot_wider(names_from = MONTH, values_from = value) %>% 
  filter(`5` ==1 & `6` ==1)

species = unique(BEF_data$SPECIES) %>%
  as.data.frame() %>% 
  rename(species = ".") %>% 
  arrange(species) 
species = species$species
## All sites and months
 for(i in species){
    graph = BEF_data %>% 
    filter(SPECIES == i) %>%
    filter(YEAR %in% c(2009:2014)) %>%
    mutate(MONTH_binned = .bincode(MONTH, month_bin)) %>% 
    group_by(YEAR, SPECIES, SITE, MONTH,MONTH_binned, DSAMP_N, EFFORT) %>% 
    summarize(CPUE_min = n()) %>% 
    mutate(CPUE_min = (CPUE_min / EFFORT)*60) %>%
    ggplot(aes(x = (YEAR),
               y = CPUE_min,
              col = as.factor(MONTH_binned))) + 
    theme_minimal() +
    geom_point() + 
    facet_wrap(~SITE) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(col = paste(i, "Month Sampled")) 
  print(graph)
}
 
## Only sites that were sampled in both may and june
## Have to remove bullhead or the for loop won't run
for(i in species[-1]){
  graph = BEF_data %>% 
    filter(SPECIES == i) %>%
    filter(YEAR %in% c(2009:2014)) %>%
    filter(SITE %in% may_sites$SITE ) %>%
    mutate(MONTH_binned = .bincode(MONTH, month_bin)) %>% 
    group_by(YEAR, SPECIES, SITE, MONTH,MONTH_binned, DSAMP_N, EFFORT) %>% 
    summarize(CPUE_min = n()) %>% 
    mutate(CPUE_min = (CPUE_min / EFFORT)*60) %>%
    ungroup() %>%
    complete(MONTH_binned, YEAR, SPECIES, SITE) %>% 
    mutate(CPUE_min = replace_na(CPUE_min, 0)) %>%
    ggplot(aes(x = (YEAR),
               y = CPUE_min,
               col = as.factor(MONTH_binned))) + 
    geom_jitter(width = .3) + 
    theme_minimal() +
    facet_wrap(~SITE) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(col = paste(i, "Month Sampled")) + 
    geom_vline(aes(xintercept = 2012))
  print(graph)
}

## Looking at differences between 2012 may/june sampling

graph_list = list()
for(i in species[-1]){
  graph = BEF_data %>% 
    filter(SPECIES == i) %>%
    filter(YEAR == 2012) %>%
    filter(SITE %in% may_sites$SITE ) %>%
    mutate(MONTH_binned = .bincode(MONTH, month_bin)) %>% 
    group_by(YEAR, SPECIES, SITE, MONTH,MONTH_binned, DSAMP_N, EFFORT) %>% 
    summarize(CPUE_min = n()) %>% 
    mutate(CPUE_min = (CPUE_min / EFFORT)*60) %>%
    ungroup() %>%
    complete(MONTH_binned, YEAR, SPECIES, SITE) %>% 
    mutate(CPUE_min = replace_na(CPUE_min, 0)) %>%
    ggplot(aes(x =as.factor(MONTH_binned),
               y = CPUE_min)) + 
    
    theme_minimal() +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(col = paste(i, "Month Sampled"))  + 
    xlab("") + 
    scale_x_discrete(labels= c("May","June")) +
    ylab(paste(i, "CPUE (min)"))
  
  graph_list[[i]] = graph
}
do.call("grid.arrange", c(graph_list[-3], ncol=3))



BEF_data %>% 
  filter(YEAR == 2012) %>%
  filter(SITE %in% may_sites$SITE ) %>%
  mutate(MONTH_binned = .bincode(MONTH, month_bin)) %>% 
  group_by(YEAR, SPECIES, SITE, MONTH,MONTH_binned, DSAMP_N, EFFORT) %>% 
  summarize(CPUE_min = n()) %>% 
  mutate(CPUE_min = (CPUE_min / EFFORT)*60 * 60) %>%
  ungroup() %>%
  complete(MONTH_binned, YEAR, SPECIES, SITE) %>% 
  mutate(CPUE_min = replace_na(CPUE_min, 0)) %>%
  rename("species" = SPECIES) %>%
  left_join(codes) %>%
  ggplot(aes(x =as.factor(MONTH_binned),
             y = CPUE_min)) + 
  
  theme_minimal() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  
  xlab("") + 
  scale_x_discrete(labels= c("May","June")) +
  ylab(paste( "CPUE (indv / hour)")) + 
  facet_wrap(~species_names, scales = "free_y")


## t.test for differences in the year 2012
for(i in species[c(-1,-4,-8)]){
  ws = BEF_data %>% 
    filter(SPECIES == i) %>%
    filter(YEAR == 2012) %>%
    filter(SITE %in% may_sites$SITE ) %>%
    mutate(MONTH_binned = .bincode(MONTH, month_bin)) %>% 
    group_by(YEAR, SPECIES, SITE, MONTH,MONTH_binned, DSAMP_N, EFFORT) %>% 
    summarize(CPUE_min = n()) %>% 
    mutate(CPUE_min = (CPUE_min / EFFORT)*60 *60) %>%
    ungroup() %>%
    complete(MONTH_binned, YEAR, SPECIES, SITE) %>% 
    mutate(CPUE_min = replace_na(CPUE_min, 0)) 
  
  a = ws %>% filter(MONTH_binned == 1)
  b = ws %>% filter(MONTH_binned == 2)
  print(i)
  print(t.test(a$CPUE_min,b$CPUE_min))

}


## t.test for differences in the year 2012 all sites
for(i in species[c(-1,-4,-8)]){
  ws = BEF_data %>% 
    filter(SPECIES == i) %>%
    filter(YEAR == 2012) %>%
    #filter(SITE %in% may_sites$SITE ) %>%
    mutate(MONTH_binned = .bincode(MONTH, month_bin)) %>% 
    group_by(YEAR, SPECIES, SITE, MONTH,MONTH_binned, DSAMP_N, EFFORT) %>% 
    summarize(CPUE_min = n()) %>% 
    mutate(CPUE_min = (CPUE_min / EFFORT)*60) %>%
    ungroup() %>%
    complete(MONTH_binned, YEAR, SPECIES, SITE) %>% 
    mutate(CPUE_min = replace_na(CPUE_min, 0)) 
  
  a = ws %>% filter(MONTH_binned == 1)
  b = ws %>% filter(MONTH_binned == 2)
  print(i)
  print(t.test(a$CPUE_min,b$CPUE_min))
  
}

v %>%
  filter(Year %in% c(2010:2014)) %>% ggplot(aes(x = Year, y = value)) +
  geom_boxplot() +
  facet_wrap(~Species, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90)) 


## After comparing the sites that were sampled during both the May and June period, it looks like there are true zeros for cold-water species in May where these species usually occur. Justifies the incorporation of the unusual June data from electrofisher malfunction.