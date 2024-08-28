### Libraries --------------

library(tidyverse)
library(vegan)

`%nin%` = Negate(`%in%`) # sets up a way to exclude if in a string


## Functions source -----------
# This is just setup at the project working directory. Use option in upper right corner of R to get into project directory. For example, on my computer,its stored in my family one-drive
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/literate-potato/")

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




post_2002 = BEF_data_unfiltered %>% filter(YEAR == 2005) %>% 
  select(SITE) %>% unique() %>%
  mutate(SITE_num = parse_number(SITE)) 


BEF_1998 = BEF_data_unfiltered %>%
  filter(YEAR == 1998) %>% 
  select(SITE) %>% unique() %>% 

  mutate(SITE_num = c(1:18)) 


bef_1998 = BEF_data_unfiltered %>% filter(YEAR == 1998) %>% 
  left_join(BEF_1998)


bef_1999 = BEF_data_unfiltered %>%
  filter(YEAR == 1999) %>% 
  dplyr::select(SITE) %>% unique() %>% 
  filter(SITE !="NA") %>% 

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






#### Data setup ---------------------
#CPUE.w.sec = rbind(CPUE.w.sec_all, bef_1999)

CPUE.w.sec = ((CPUE_wide_seconds(BEF_data) %>%
                 unite("Group", c(YEAR, SITE)) %>% 
                 column_to_rownames(., var = "Group") %>% 
                 mutate(sumrow = rowSums(.)) %>%
                 filter(sumrow>0) %>%
                 dplyr::select(-sumrow)))


## Write this to the Data folder
# This folder is included in the .gitignore so it will not show up on github
write.csv(CPUE.w.sec, file="Data/CPUE.w.sec.csv")


## Write BEF_Data.csv for LML_P1_community-response.R
write.csv(BEF_data, file = "Data/BEF_Data.csv")
