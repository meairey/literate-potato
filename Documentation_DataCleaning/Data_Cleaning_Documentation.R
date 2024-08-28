## Assessing impacts of year 2000 depletion through time on SMB populations 

## Please load in data associated with LML_P1_community-response.R to get the BEF_data data frame
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
#install.packages("lme4")
#install.packages("Matrix")
library(Matrix)
library(lme4)
## Overall Depletion 
oo <- options(repos = "https://cran.r-project.org/")
#install.packages("Matrix")
#install.packages("lme4")
options(oo)
library(Matrix)
library(lme4)
#install.packages("lme4", type = "source")
library(lme4)
library(nlme)
library(gridExtra)

BEF_data = BEF_data_unfiltered %>%
  filter(SPECIES %nin% c(stocked, rare$SPECIES)) %>% 
  filter(YEAR < 2020) %>% 
  filter(SPECIES != "SMB" | YEAR != 2000 | DAY_N < 153) %>% 
  filter(YEAR != 2002) #%>% 
#filter(YEAR != 1999)## Filter out BEF SMB data from the year 2000 that's later than DAY_N 160. Change this around depending on how you want to filter 2000... 

#1:length(species))[c(-7,-8)] put this bck in for loop

graph_list_lme2000 = list()
graph_list_lme2001 = list()
for(i in (1:length(species))[c(-7,-8)]){

  species_interest = species[i]
  
  lme.data_2000 = BEF_data %>% filter(YEAR == 2000) %>%
    filter(SPECIES == species_interest) %>% 
    select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
    group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
    summarize(total_count = n()) %>% 
    mutate(CPUE = (total_count / EFFORT)*60*60) %>%
    ungroup() %>% 
    mutate(SITE = (as.factor(SITE))) %>% 
    mutate(CPUE = as.numeric(CPUE)) %>%
    na.omit()

  j_2000 = lme.data_2000 %>% 
    na.omit() %>%
    ggplot(aes(x = DAY_N, y = CPUE)) + 
    geom_point() +
    geom_smooth(method = lm, col = "black") + 
    theme_minimal() + 
    ylab(paste(species_interest, "CPUE")) +
    xlab("") +
    xlim(134, 175)
    

  lme_2000 = nlme::lme(CPUE~DAY_N, random = ~1 | SITE, data = lme.data_2000)
  print(species[i])
  
  print(summary(lme_2000))
  graph_list_lme2000[[i]] = j_2000

}

do.call("grid.arrange", c(graph_list_lme2000[c(-7,-8,-11)], ncol=3))

for(i in (1:length(species))[-11]){
  
  species_interest = species[i]
  lme.data_2001 = BEF_data %>% filter(YEAR == 2001) %>%
    filter(SPECIES == species_interest) %>% 
    select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
    group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
    summarize(total_count = n()) %>% 
    mutate(CPUE = (total_count / EFFORT)*60*60) %>%
    ungroup() %>% 
    mutate(SITE = (as.factor(SITE))) %>% 
    mutate(CPUE = as.numeric(CPUE))
  
  j_2001 = lme.data_2001 %>% ggplot(aes(x = DAY_N, y = CPUE)) + 
    geom_point() +
    geom_smooth(method = lm, color = "black") + 
    theme_minimal() + 
    ylab(paste(species_interest, "CPUE")) +
    xlab("") +
    xlim(128, 170)

  
  lme_2001 = nlme::lme(CPUE~DAY_N, random = ~1 | SITE, data = lme.data_2001)
  print(species[i])
  
  print(summary(lme_2001))
  graph_list_lme2001[[i]] = j_2001
}

BEF_data %>% filter(YEAR == 2000) %>%
  select(SPECIES, SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(SPECIES, DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60*60) %>%
  ungroup() %>% 
  mutate(SITE = (as.factor(SITE))) %>% 
  mutate(CPUE = as.numeric(CPUE)) %>%
  rename("species" = SPECIES) %>%
  left_join(codes) %>%
  ggplot(aes(x = DAY_N, y = CPUE)) + 
  geom_point() +
  geom_smooth(method = lm, color = "black", se = F) + 
  theme_minimal() + 
  ylab("CPUE (indv / hour)") +
  xlab("") +
  #xlim(128, 170) + 
  facet_wrap(~species_names, scales = "free_y")



do.call("grid.arrange", c(graph_list_lme2001[-11], ncol=4))

## Temperature gradient over 2000 and 2001

BEF_data %>% filter(YEAR %in% c(2000, 2001)) %>% 
  select(YEAR, DAY_N, TEMP_SAMP) %>% unique()

BEF_data %>% 
  select(YEAR, TEMP_SAMP) %>% unique()

##It looks like there is quite a bit of depletion, but I think that the graph is misleading...
## Don't forget to double check the filter on BEF_data to get the range of dates you want included
species_interest = "PS"

BEF_data %>% filter(YEAR == 2001) %>%
  filter(SPECIES == species_interest) %>% 
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) %>% 
  ggplot(aes(x = DAY_N, y = log10(CPUE))) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  ylab("CPUE Ind/Min") 

BEF_data %>% filter(YEAR == 2000) %>%
  filter(SPECIES == species_interest) %>% 
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) 
## Observation - some sites show depeletion, other do not show depeletion. 
#### Many of those sites look like there were not many fish to depelete anyways? Compared to site 001?
#### Note some sites (like 022 and 023 increase)

BEF_data %>% filter(YEAR == 2001 & SPECIES == species_interest) %>%
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) %>% 
  ggplot(aes(x = DAY_N, y = CPUE)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~SITE) + 
  ylab("CPUE Ind/Min") +
  theme(axis.text.x = element_text(angle = 90))

## It looks like the first day of sampling for sites varied by as much as 20 days
BEF_data %>% filter(YEAR == 2000 & SPECIES == species_interest) %>%
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) %>% 
  ungroup() %>% 
  group_by(SITE) %>% 
  summarize(first_day = first(DAY_N)) %>%
  ggplot(aes(x = order(SITE, first_day), y = first_day)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Site") + ylab("First Day of Sampling")

## Looking at depletion by site            
depletion = BEF_data %>% filter(YEAR == 2000 & SPECIES == species_interest) %>%
  arrange(SITE) %>%
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) %>% 
  group_by(SITE) %>%
  do(day_regression = lm(CPUE ~ DAY_N , data = .))
  

## Not many of these are significant... not even site 1. These aren't normal data so this isn't perfect. But I'm not sure there is statistical proof of depletion through time at these sites

## Check below - these should make data frames that list sites with significant depletion


p.value = lapply(depletion$day_regression, function(x) summary(x)) %>%
  lapply(., function(x) p.value = max(x$coefficients[,4])) %>% unlist() %>%
  as.data.frame() %>% 
  rename("p.value" = ".") %>% mutate(depletion$SITE) %>% filter(p.value < .05)

slope = lapply(depletion$day_regression, function(x) summary(x)) %>%
  lapply(., function(x) slope = x$coefficients[,1][2])%>% unlist() %>%
  as.data.frame() %>% 
  rename("slope" = ".") %>% mutate(depletion$SITE) 

left_join(p.value, slope)

## Without filter there are 7/28 (25%) sites with significant depletion, with filter there is only 1 significant depletion
## Explore through the above regressions using....

site_depletion = "005" ## enter the site here that you are interested in. It has to be a character. Then you get summary of regression
summary(depletion[[2]][[which(depletion$SITE == site_depletion)]])


## Conclusion - I'm not sure there is any correct answer here... SMB are the treatment, not a response variable. I'd either pick averaging across all data points for 2000 or filtering within a date range like 130 - 160. The first value from each site seems like a bad representation of the data given that a site first sampled on Day 135 != CPUE from site first sampled day 145+



## First vs. last 

first = BEF_data %>% filter(YEAR == 2000 & SPECIES == species_interest) %>%
 # filter(SITE %nin% c("001", "005")) %>%
  filter(DAY_N < 160) %>%
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) %>%
  ungroup() %>% 
  group_by(SITE) %>%
  summarize(first = first(CPUE)) %>% na.omit()
last = BEF_data %>% filter(YEAR == 2000 & SPECIES == species_interest) %>%
  filter(DAY_N < 160) %>%
  #filter(SITE %nin% c("001", "005")) %>%
  select(SITE, EFFORT, DSAMP_N, DAY_N) %>% 
  group_by(DAY_N, DSAMP_N, SITE, EFFORT) %>% 
  summarize(total_count = n()) %>% 
  mutate(CPUE = (total_count / EFFORT)*60) %>%
  ungroup() %>% 
  group_by(SITE) %>%
  summarize(first = last(CPUE)) %>% na.omit()

wilcox.test(first$first,last$first, alternative = "less")
t.test(last$first, mu = mean(first$first), alternative = "less")

wilcox.test(first$first,last$first, alternative = "two.sided")
t.test(last$first, mu = mean(first$first), alternative = "less")


, mu = 4.34, "less")

cbind(first$first, last$first) %>%
  as.data.frame() %>% 
  mutate(site = c(1:length(last$first))) %>% 
  pivot_longer(c(V1:V2), names_to = "first_last", values_to = "CPUE") %>%
  ggplot(aes(x = first_last, y = CPUE)) + 
  geom_boxplot()
