## Supplamental Statistics for Supp. Section ----------------------------


pre = BEF_data_unfiltered %>% filter(YEAR < 2001) %>% filter(SPECIES =="SMB") 
post = BEF_data_unfiltered %>% filter(YEAR > 2001) %>% filter(SPECIES == "SMB")
before_after = c(1996,2000, 2025)

pre$WEIGHT %>% na.omit() %>% median()
post$WEIGHT %>% na.omit() %>% median()
spring_fall = c(0,8,11)
## Summary stats for supplemental 
BEF_data %>% select(MONTH, YEAR, SITE,EFFORT, DAY_N) %>% 
  mutate(MONTH_color = .bincode(MONTH, spring_fall)) %>%
  unique() %>% group_by(YEAR) %>% 
  summarize(m = mean(DAY_N, na.rm = T)) %>% 
  mutate(Year_bin = .bincode(YEAR, before_after)) %>%
  ungroup() %>% group_by(Year_bin) %>%
  summarize(m = sd(m))
## Parsing site names 
BEF_data %>% select(MONTH, YEAR, SITE,EFFORT, DATE_COL) %>% 
  unique() %>% 
  group_by(SITE,YEAR)


## Pre vs post effort mean of total effort per year  

BEF_data %>% select(MONTH, YEAR, SITE,EFFORT, DATE_COL) %>% 
  unique() %>% 
  group_by(SITE,YEAR) %>% 
  mutate(YEAR_prepost = .bincode(YEAR, before_after)) %>% 
  
  ungroup() %>% 
  group_by(YEAR, YEAR_prepost) %>% 
  summarize(sum = sum(EFFORT, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(YEAR_prepost) %>% 
  summarize(mean = sd(sum))

BEF_data %>% select(YEAR, SITE) %>% unique() %>% 
  group_by(YEAR) %>% 
  summarize(count = n())


# Log10 of sum of effort --- Supplemental -----------
BEF_data %>% select(MONTH, YEAR, SITE,EFFORT, DATE_COL) %>% 
  unique() %>% 
  group_by(SITE,YEAR) %>%
  summarize(effort = sum(EFFORT)) %>% 
  mutate(effort = log10(effort)) %>%
  ggplot(aes(x =(YEAR), 
             y = effort)) + 
  geom_point() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("log10 Effort (s)") + 
  xlab("Year")

## Trying to break it up by year instead of by year,site 
### color scale the day of year and effort total per day per year



BEF_data %>% select(MONTH, YEAR, SITE,EFFORT, DATE_COL, DSAMP_N) %>% 
  unique() %>% 
  group_by(YEAR,DATE_COL, DSAMP_N) %>%
  summarize(effort = mean(EFFORT)) %>% 
  mutate(EFFORT = log10(effort)) %>%
  ggplot(aes(x =(YEAR), 
             y = EFFORT)) + 
  geom_point() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("log10(Effort)") + geom_smooth()
## I need to fix sites here - its a little weird with the partial sites 

## Trying to make a table of effort/site/year information 

BEF_data %>% dplyr::select(YEAR, MONTH, SITE, EFFORT, DATE_COL, DSAMP_N)

## Julian day --- Supplemental -----
BEF_data %>% 
  dplyr::select(MONTH, YEAR, SITE,EFFORT, DAY_N, GEAR_CODE) %>% 
  group_by(YEAR, DAY_N) %>%
  summarize(EFFORT = sum(EFFORT)) %>%
  ggplot(aes(x = YEAR,
             y = DAY_N,
             col = round(EFFORT,digits = 2)#, 
             #col = as.character(as.factor(GEAR_CODE)),
             #shape = MONTH_color
             
  )) +
  theme_minimal() +
  geom_vline(xintercept = 2000, linetype = 2) +
  #geom_jitter(width = .2, size = 1.2, height = 0) + 
  geom_point() + 
  ylab("Day of Year") + 
  labs(color = "Total Effort (sec)") +  
  #scale_shape_discrete(labels = c("Spring", "Fall")) + 
  xlab("") + 
  theme(axis.text = element_text(size = 13), 
        #text = element_text(size = 13), 
        axis.text.x = element_text(angle = 90, vjust = .5)) +
  
  #ylim(123,175) + 
  xlim(1998, 2019) + 
  scale_color_viridis_c(trans = "log", breaks = c(2980, 22026, 162754, 1202602))
    
    #+
  #guides(colour = guide_legend(override.aes = list(size=3))) + 
  #guides(shape = guide_legend(override.aes = list(size=3)))  
  #scale_color_manual(values = cbbPalette, labels = c("Day Bass Only", "Night All Fish", "Night Bass Only")) + xlim(1998, 2019) + ylim(100, 325)
  
##


### Julian Day -- FBL and LML

BEF_data_unfiltered %>%
  select(WATER, YEAR, MONTH, SITE_N, EFFORT, DAY_N, GEAR_CODE) %>%
  group_by(WATER, YEAR, MONTH, DAY_N, GEAR_CODE) %>%
  summarize(EFFORT = sum(EFFORT)) %>%
  filter(GEAR_CODE == "NAF", MONTH < 7) %>%
  ggplot(aes(x = YEAR,
             y = DAY_N,
             col = round(EFFORT,digits = 2))) +
  theme_minimal() +
  geom_vline(xintercept = 2000, linetype = 2, WATER == "LML") +
  #geom_jitter(width = .2, size = 1.2, height = 0) + 
  geom_point() + 
  ylab("Day of Year") + 
  labs(color = "Total Effort (sec)") +  
  #scale_shape_discrete(labels = c("Spring", "Fall")) + 
  xlab("") + 
  theme(axis.text = element_text(size = 13), 
        #text = element_text(size = 13), 
        axis.text.x = element_text(angle = 90, vjust = .5)) +
  
  #ylim(123,175) + 
  xlim(1998, 2019) + 
  scale_color_viridis_c(trans = "log", breaks = c(2980, 22026, 162754, 1202602)) + 
  facet_wrap(~WATER)


julian_dat_cvs = BEF_data %>% dplyr::select(MONTH, YEAR, SITE,EFFORT, DAY_N, GEAR_CODE) %>% 
  mutate(MONTH_color = as.character(as.factor(.bincode(MONTH, spring_fall)))) %>%
  unique() %>% 
  select(YEAR, MONTH, DAY_N, EFFORT) %>% 
  group_by(YEAR, MONTH, DAY_N) %>% 
  summarize(total_effort = sum(EFFORT))

#write.csv(julian_dat_cvs, file = "julian_dat.csv")
#dev.off()


## means 
BEF_data %>% select(MONTH, YEAR, SITE,EFFORT, DAY_N) %>% 
  mutate(MONTH_color = .bincode(MONTH, spring_fall)) %>%
  unique() %>% group_by(YEAR) %>% 
  summarize(m = mean(DAY_N)) %>% 
  mutate(Year_bin = .bincode(YEAR, before_after)) %>%
  ungroup() %>% group_by(Year_bin) %>%
  summarize(m = sd(m))



BEF_data %>% mutate(MONTH_color = .bincode(MONTH, spring_fall))
## Looking to see if CPPUE of SMB through time 
BEF_data %>% 
  filter(SPECIES == "SMB") %>% 
  filter(YEAR == 2000) %>% 
  group_by(DAY_N, SITE, YEAR, EFFORT) %>%
  
  summarize(count = n()) %>%
  mutate(cpue = count/EFFORT) %>%
  ggplot(aes(x = DAY_N, y = cpue )) + geom_point()

## Notes from meeting with Kurt 
# Go through and add back in the 2012 data
# Go through and make sure to weed out the later 2000 dates for just the SMB

## Checking # of sites through time
