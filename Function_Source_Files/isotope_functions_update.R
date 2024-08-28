





### CV of different fish groups 
# ---------------------- Library Loads ---------------------------------------
library(RColorBrewer)
library(abind)

#----------------------------- Params ---------------------------------------- 

x.p = seq(-38,-20,length.out=100); 
y.p = seq(0,9, length.out = 100) 
parms <- list();
parms$n.iter=2*10^4; 
parms$n.burnin=1*10^3; 
parms$n.thin=10
parms$n.chains=3     
priors=list();
priors$R=1*diag(2);
priors$k=2; 
priors$tau.mu=1.0E-3 # Vague priors 
Nsamples=1000
n.posts <- 1000;
p.ell <- 0.90 # How much data to include? Standard ellipses --> p.ell = .9
n.points = 1000




overlap = function(data_input, comm, dr, posterior){
  ## Remove this later
  data_overlap = data_input %>% 
    filter(community == comm)
  
  #data_overlap = combo %>% filter(community == 1)
  spp=length(unique(data_overlap$group))
  print(data_overlap)
  
  name_list = names(posterior) 
  
  
  name_matrix = expand.grid(name_list, name_list) %>% 
    separate(Var1, into = c("C1","G1"), remove = F) %>%
    separate(Var2, into = c("C2", "G2"), remove = F) %>%
    filter(C1 == comm & C2 == comm) %>%
    select(Var1, Var2) %>%
    filter(Var1 != Var2)
  
  
  names_modified = expand.grid(name_list, name_list) %>% 
    separate(Var1, into = c("C1","group"), remove = F) %>%
    separate(Var2, into = c("C2", "G2"), remove = F)%>%
    filter(C1 == comm & C2 == comm) %>%
    mutate(group = as.numeric(group))%>%
    left_join(legend) %>%
    select(-color, -common) %>%
    rename("G1" = group) %>%
    rename(group = "G2") %>%
    rename(Code1 = CODE) %>%
    mutate(group = as.numeric(group)) %>%
    left_join(legend) %>%
    select(-color, -common) %>%
    rename(G2 = group) %>%
    unite(col = Name, c(Code1, CODE), sep = " v " ) %>%
    filter(Var1 != Var2) %>%
    select(Name)
  
  
  ma = matrix(NA, nrow = dr, ncol = length(name_matrix[,1]))
  
  for(i in 1:length(name_matrix[,1])){ 
    
    overlap = bayesianOverlap(ellipse1 = as.character(name_matrix[i,1]), 
                              ellipse2 = as.character(name_matrix[i,2]), 
                              posterior, 
                              draws = dr,
                              p.interval = 0.95,
                              n = dr)
    
    bayes.prop <- (overlap[,3] / (overlap[,2] + overlap[,1] - overlap[,3]))
    
    ma[,i] = bayes.prop
    
  }
  
  
  
  full_olap_mat = ma %>% as.data.frame() %>%
    rename_with(~ setNames((unique(names_modified[,1])), .), everything())
  
  full_olap_mat[dr+1,] = unique(names_modified[,1])
  
  rownames(full_olap_mat) = c((data.frame(n = "com",c = comm, post = c(1:dr)) %>% unite("rowname", c(n, c, post)))$rowname, "SppPair")
  
  
  #return(as.data.frame(t(olap_mat))) ## This is a reminder to pull 90% CI intervals out of this at a point when i can play with code and don't need to be writing
  return(t(as.data.frame(full_olap_mat)))
  
}



## -------------- Data function ----------------------------

## Removed PD because not sure identification 
## Must put in Siber formatted data_input

data_setup = function(data_input, com_num){
  # Pre combo_lake = combo_2023 %>% filter(Water == lake$Water[[com_num]])
  combo_lake = data_codes %>% 
    filter(community == com_num)
  #filter(community == lake$Community[[com_num]])
  data = data_input %>% 
    filter(community == com_num)
  data = data[order(data$group),] %>% as.data.frame() 
  siber.example <- createSiberObject(data)
  posterior <- siberMVN(siber.example, parms, priors)
  
  both = list(siber.example, posterior, data,combo_lake)
  return(both)
}


## ------------------------ Ellipse data function ---------------------
ellip_data = function(numb_species, numb_posts, posterior){
  for(i in 1:numb_species){
    dat = vector()
    for(j in 1:numb_posts){
      ellipse_data =  ellipse::ellipse(x = matrix(c(posterior[[2]][[i]][j,1],
                                                    median(posterior[[2]][[i]][j,2]),
                                                    median(posterior[[2]][[i]][j,3]),
                                                    median(posterior[[2]][[i]][j,4])),
                                                  2,2), 
                                       centre = c(median(posterior[[2]][[i]][j,5]),
                                                  median(posterior[[2]][[i]][j,6])), level = .95, 
                                       npoints = n.points) %>% as.data.frame() %>%
        summarise_all(list("min"=min, "max"=max)) %>% as.matrix()
      
      dat = rbind(dat, ellipse_data)
    }
    
    ellip[,,i] = dat
    
  }
  
  return(ellip)
}

## Legends

## Setting up data frames --------------------------------
#### Includes lipid corrections

mwt_20 = read.csv(file = "Data/2020_mwt.csv") %>% 
  filter(FISH_N != "") %>%
  unique()  %>%
  filter(MA_ID != 50) %>%
  select(FISH_N, LENGTH, WEIGHT, MA_ID)
mwt_21 = read.csv(file = "Data/AIREY_MWT_2021.csv")%>% 
  filter(FISH_N != "") %>%
  unique() %>%
  select(FISH_N, LENGTH, WEIGHT)

m = read.csv(file = "Data/isotope_measurement.csv") %>% 
  left_join(mwt_20, by = "FISH_N") %>%
  left_join(mwt_21, by = "FISH_N") %>%
  filter(ISO_YSAMP_N != "")

bef = read.csv("Data/FISH_MEASUREMENT_LML.csv") %>% group_by(FISH_N) %>%
  slice(1)
gln = read.csv("Data/FISH_MEASUREMENT_PRL.csv")

s = read.csv(file = "Data/isotope_sample.csv")%>%
  filter(ISO_YSAMP_N != "")

full_data1 = left_join(m, s, by = "ISO_YSAMP_N") %>%
  mutate(TAXA = str_replace(TAXA,"BM","BNM")) %>%
  mutate(TL = case_when(is.na(LENGTH.x)== T ~ as.numeric(LENGTH.y),
                        is.na(LENGTH.x) == F ~ as.numeric(LENGTH.x))) %>% 
  left_join(bef %>% filter(FISH_N != ""), by = "FISH_N") %>%
  mutate(TL = case_when(is.na(TL)== T ~ as.numeric(LENGTH),
                        is.na(TL) == F ~ as.numeric(TL))) %>% 
  select(-LENGTH.y, -LENGTH.x, -WEIGHT.y, -WEIGHT.x, -LENGTH, -WEIGHT) %>%
  left_join(gln %>% filter(FISH_N != ""), by = "FISH_N") %>% 
  mutate(TL = case_when(is.na(TL)== T ~ as.numeric(LENGTH),
                        is.na(TL) == F ~ as.numeric(TL))) %>%
  mutate(C.N = PER_C / PER_N) %>%
  mutate(D13C = D13C - 3.32 + .99 * C.N)


## Incorporating  new data
dd_communities = c(2,3,22)
dumb_fix = c(1:23)[-dd_communities]

## Filter out just fis data and replace names that need to be changed
full_data = full_data1 %>% filter(CATEGORY == "FISH") %>%
  filter(MATERIAL == "TISSUE") %>%
  select(WATER, YEAR, MONTH, everything()) %>%
  mutate(SEASON = case_when(MONTH %in% c(5:7) ~ "spring", MONTH %in% c(8:10) ~ "fall")) %>%
  unite(ID, c(WATER, YEAR, SEASON), remove=F) %>%
  mutate(community = (as.factor(ID))) %>%
  mutate(community= as.numeric(community)) %>%
  mutate(TAXA = str_replace(TAXA, "FSD-NRD", "NRD")) %>%
  mutate(TAXA = str_replace(TAXA, "FSD", "NRD")) 


## Filter out just fish data and replace names that need to be changed





legend = full_data %>% select(TAXA) %>% arrange(TAXA) %>% unique() %>%
  filter(TAXA != "") %>%
  rename(CODE = TAXA) %>%
  mutate(common =  c("Brown Bullhead", "Banded Killifish", "Blacknose Dace","Bluntnose Minnow", "Creek Chub", "Cutlips Minnow", "Common Shiner", "Fathead Minnow", "Golden Shiner", "Landlocked Salmon","Central Mudminnow","N. Pearl Dace", "Chrosomus spp.", "Pumpkinseed Sunfish","Redbreasted Sunfish", "Rainbow Smelt","Smallmouth Bass", "Brook Trout", "White Sucker")) %>% 
  mutate(group = as.numeric(as.factor(CODE)),
         color = as.numeric(as.factor(CODE)))


data = full_data %>%
  ungroup() %>%
  select(community,TAXA, D13C, D15N) %>%
  rename(CODE = TAXA) %>%
  left_join(legend) %>%
  select(D13C, D15N, group, community) %>%
  rename(iso1 = D13C, iso2 = D15N) %>% ## Format data to go into siber object 

  mutate(group = as.numeric(group)) %>%
  arrange(community, group)  %>% ## Arrange data in order
  filter(group != "NA") %>%
  group_by(community, group) %>%
  mutate(count = n()) %>%
  filter(count > 3) %>%
  select(-count) %>%
  filter(community %in% dumb_fix)

# set up lake matrix for look up tables
lake_2023_new = full_data %>%
  select(community, ID) %>%
  unique() %>% 
  arrange(community)

