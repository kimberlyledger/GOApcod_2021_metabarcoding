
```{r}
sites <- metadata %>%
  dplyr::select(Sample_ID, location1)

rep_table <- asv_table_filter3 %>%
  left_join(sites, by = "Sample_ID") %>%
  arrange(location1)
```

```{r}
#library(stringi)
tmp.df <- rep_table %>%
  filter(location1 == 3) # comment this to test the ASV-centric df
```

```{r}
# format the dataframe appropriately
wide.loc.frame <- tmp.df %>%
  ungroup() %>%  
  mutate(reads = ifelse(reads > 0, 1, 0)) %>% # change counts to presence/absence
  dplyr::select(ASV, Sample_ID, reads) %>%
  group_by(ASV, Sample_ID, reads) %>%
  unique() %>% # eliminate duplicate entries for the same ASV and sample (these seem to be a thing for aquaF2)
  pivot_wider(names_from = Sample_ID, values_from = reads) %>%
  ungroup()
```


```{r}
#library(rstan)

# format must be dataframe, not tibble for converting rownames properly
wide.loc.frame <- data.frame(wide.loc.frame)

# helper function for maintaining rownames in matrix format
matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

# convert df to matrix
wide.loc.matrix <- matrix.please(wide.loc.frame)

ASV <- row.names(wide.loc.matrix)

n_replicates = ncol(wide.loc.matrix)

# then call the Stan code to work on the dataset we created above. 
# reformatting our data from above to match the Stan code's inputs
testData <<- 
  data.frame(K = n_replicates,   #trials per species (row)
             N = rowSums(wide.loc.matrix),  #detections per species
             z = ifelse(rowSums(wide.loc.matrix) > 0, 1, 0)  #was it ever detected at this site? (integer that helps estimate psi)
  )

mySOMmodel <- stan(file = "Stan_SOM_demo.stan",
                   data = list(
                     S = nrow(testData),
                     K = testData$K,
                     N = testData$N,
                     z = ifelse(testData$N > 0, 1, 0)
                   ),
                   chains = 4,   #number of chains
                   iter = 5000   #number of iterations per chain
)

```

```{r}
library(broom)
library(tibble)

mySOMmodel_results2 <- tidy(tibble(as.data.frame(mySOMmodel)))
```

```{r}
library(bayesplot)

p1 <- mcmc_areas(mySOMmodel, 
                 pars = c("psi", "p11", "p10"))
p2 <- mcmc_intervals(mySOMmodel, 
                     regex_pars = "Occupancy_prob")

plot(p1)
plot(p2)
```

```{r}
#See histogram of occupancy probabilities for unique patterns of presence
mySOMmodel_results %>% 
  dplyr::filter(grepl("Occupancy_prob", column)) %>% 
  ggplot() + geom_histogram( aes(mean)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "Unique Patterns of Presence") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_4

#plot_3
plot_4
#ggsave(plot=plot_3,"Output_plots/Histogram_of_Mean_Occupancy_Probabilities_unique_patterns_of_presence.png", device = "png", width = 12, height = 8, units = "in")
```



```{r}
#Now Re-attach SOM to unique patterns
mySOMmodel_results %>% 
  dplyr::filter(grepl("Occupancy_prob", column)) %>%
  separate(column, into=c("column","N"), sep="([\\[\\]])") %>%    ### i'm guessing this is the column i want to join later on so i'm naming it "N"
  tidyr::pivot_wider(., names_from=column, values_from = c(mean, sd))-> mySOMmodel_results_wide

mySOMmodel_results_wide$N <- as.numeric(mySOMmodel_results_wide$N)

# this chunk is from gruinard_decon  - think i'll need to use this once I have data from multiple sites... 
#testData %>% 
#  dplyr::left_join(mySOMmodel_results_wide) %>% 
#  group_by(seq_number,Site) %>% 
#  summarise(., max_Occupancy_prob = max(mean_Occupancy_prob)) -> unique_data_SOM

joinedData <- testData %>% 
  dplyr::left_join(mySOMmodel_results_wide) %>%
  cbind(ASV) # i loss the row names here so add them back in... (they are the ASVs)

hist(joinedData$mean_Occupancy_prob)
```

now i will need to save the mean occupancy prob for each ASV... and then figure out how to do this over and over again for every site.  




```{r}
######################    
#create test data
######################

Nspecies <- 2
Ndates <- 3
Nsites <- 3
Ntechreplicates <- 3

# mock data, given true parameters and using those to generate observations            
psi_given <- seq(0.3, 0.9, length.out = Nspecies*Nsites) %>% 
  sample()
p11_given <- seq(0.7, 0.9, length.out = Nspecies) %>% 
  sample()
p10_given <- seq(0.01, 0.05, length.out = Nspecies) %>% 
  sample()

#set up data with different hierarchical levels, in case you want to play with these later
testData <- expand.grid(
  "Species" = c(1:Nspecies), # model N species
  "Date" = c(1:Ndates), # samples collected on different dates
  "Site" = c(1:Nsites) # at different sites
) %>%
  mutate(
    K = Ntechreplicates,
    N = NA
  ) %>%
  arrange(Species)

#create unique identifier for combinations of site-date-species; for use in hierarchical modeling
SDS <- unite(data = testData,
             col = SDS,
             c("Site", "Date", "Species")
) %>% pull(SDS)
testData$SiteDateSpecies <- match(SDS, unique(SDS)) #index for unique site-date-species combinations

#create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = testData,
            col = SS,
            c("Site", "Species")
) %>% pull(SS)
testData$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations

testData <- testData %>% 
  mutate(psi_given = psi_given[testData$SiteSpecies],
         p11_given = p11_given[testData$Species],
         p10_given = p10_given[testData$Species])

#Given parameters (psi, p11, p10) assigned above, generate pattern of detections in mock data
for (i in 1:nrow(testData)){
  if(rbinom(1,1,testData$psi_given[i]) == 0) {testData$N[i] <- 0} else {
    testData$N[i] <- rbinom(1, testData$K[i], testData$p11_given[i])
  }
}

```


#the stan model used here: https://github.com/zjgold/gruinard_decon/blob/master/gruinard_decontam_script.R
```{r}
##Stan Model
sink("Stan_SOM_hierarchical.stan")
cat(
  "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // number of detections among these replicates
    int z[S];   // integer flag to help estimate psi parameter
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0,upper=1> psi[Nloc];  //commonness parameter
    real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
    real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
}

transformed parameters{/////////////////////////////////////////////////////////////////////
}

model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
    for (i in 1:S){
			z[i] ~ bernoulli(psi[L[i]]);
			p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
			N[i] ~ binomial(K[i], p[i]);
	}; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
}

generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015
  
  for (i in 1:S){
  Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  / ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  + (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
  );
  }
 }
  
",
fill=TRUE)
sink()
```


the stan model
```{r}
##Stan Model
sink("Stan_SOM_hierarchical.stan")
cat(
  "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // number of detections among these replicates
    int z[S];   // integer flag to help estimate psi parameter
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0,upper=1> psi[Nloc];  //commonness parameter
    real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
    real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
}

transformed parameters{/////////////////////////////////////////////////////////////////////
}

model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
    for (i in 1:S){
			z[i] ~ bernoulli(psi[L[i]]);
			p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
			N[i] ~ binomial(K[i], p[i]);
	}; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
}

generated quantities{
}

",
fill=TRUE)
sink()
```

run the model 
```{r}
#####################
#run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++
#####################      
myHierarchicalModel <- stan(file = "Stan_SOM_hierarchical.stan", 
                            data = list(
                              S = nrow(occu_df_TEST_ASV1),
                              Species = occu_df_TEST_ASV1$Species,
                              Nspecies = length(unique(occu_df_TEST_ASV1$Species)),
                              L = occu_df_TEST_ASV1$SiteSpecies,
                              Nloc = length(unique(occu_df_TEST_ASV1$SiteSpecies)),
                              K = occu_df_TEST_ASV1$K,
                              N = occu_df_TEST_ASV1$N,
                              z = ifelse(occu_df_TEST_ASV1$N > 0, 1, 0)
                            ), 
                            chains = 4,   #number of chains
                            iter = 4000   #number of iterations per chain
)

myHierarchicalStanResults <- tidy(tibble(as.data.frame(myHierarchicalModel)))   
plot(myHierarchicalModel)
```

okay, so now i have 83 psi values - the occupancy of ASV1 at the site - and a p11 (true positive) and p10 (false positive rate)

plot a histogram of occupancy probabilities for ASV1
```{r}
myHierarchicalStanResults%>% 
  filter(grepl("psi", column)) %>% 
  ggplot() + geom_histogram( aes(mean)) +
  ggtitle(label = "Histogram of Mean Psi", subtitle = "Unique Patterns of Presence") +
  ylab("Count") + xlab("Mean Psi")
```

probably want to use a 0.7 probability cut off (this is just according to ASV1) . 

now i need to save this output to reconnect to my filtered ASV dataframe
```{r}
myHierarchicalStanResults <- myHierarchicalStanResults %>%
  filter(grepl("psi", column)) %>%
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])")

myHierarchicalStanResults$SiteSpecies <- as.numeric(myHierarchicalStanResults$SiteSpecies)
```


```{r}
asv1_psi <- occu_df_TEST_ASV1 %>% 
  select(Species, Site, SiteSpecies) %>%
  unique() %>%
  left_join(myHierarchicalStanResults, by = "SiteSpecies") %>% 
  group_by(Species, Site)
asv1_psi
```

######### OR DO I DO THIS???????? ##########








############# need to write "sodm.by.AVS" function #############

for testing, filter to just ASV1 
```{r}
#occu_df_ASV1 <- occu_df %>%
#  filter(Species == 1)

occu_df_TEST <- occu_df %>%
  filter(Species == 1)
```

create unique identifier for combinations of site-biologicalrep-ASV; for use in hierarchical modeling
```{r}
SDS <- unite(data = occu_df_TEST, col = SDS, c("Site", "BiologicalRep", "Species")) %>% pull(SDS)

occu_df_TEST$SiteRepSpecies <- match(SDS, unique(SDS)) #index for unique site-biologicalrep-species combinations
```

create unique identifier for combinations of site-ASV; for use in hierarchical modeling
```{r}
SS <- unite(data = occu_df_TEST, col = SS, c("Site", "Species")) %>% pull(SS)

occu_df_TEST$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations
```


the stan code does not want to run all together, but seems to work fine if i split the dataframe into ASVs before running it... 

run the model - this only works with Species == 1... maybe something gets messed up with the order/numbers/etc... 
STAN model requires species names to be numbers.  Maybe there also can't be missing numbers...
```{r}
       #####################
       #run Stan model
       #note this will take a while the first time you run a particular model, because it needs to compile from C++
       #####################      
       myHierarchicalModel <- stan(file = "Stan_SOM_hierarchical_with_occuprob.stan", 
                             data = list(
                               S = nrow(occu_df_TEST),
                               Species = occu_df_TEST$Species,
                               Nspecies = length(unique(occu_df_TEST$Species)),
                               L = occu_df_TEST$SiteSpecies,
                               Nloc = length(unique(occu_df_TEST$SiteSpecies)),
                               K = occu_df_TEST$K,
                               N = occu_df_TEST$N,
                               z = ifelse(occu_df_TEST$N > 0, 1, 0)
                             ), 
                             chains = 4,   #number of chains
                             iter = 4000   #number of iterations per chain
       )
       
       myHierarchicalStanResults <- tidy(tibble(as.data.frame(myHierarchicalModel)))   
       plot(myHierarchicalModel)
```

calculate mean occupancy probabilty for each SiteSpecies using the OccuProb values generated for SiteRepSpecies 

plot - these values are all a bit larger than psi... maybe a 0.8 cutoff for ASV1 would be okay here 
```{r}
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>% 
  ggplot() + geom_histogram( aes(mean)) +
  ggtitle(label = "Histogram of Occupancy Probabilities") +
  ylab("Count") + xlab("Occupancy Probability")
```

the output has 248 occupancy probabilities that i need to summarize 
probably would be good to save my p11 and p10 for each ASV

now i need to save this output to reconnect to my filtered ASV dataframe
```{r}
myHierarchicalStanResults <- myHierarchicalStanResults %>%
  filter(grepl("Occupancy_prob", column)) %>%
  separate(column, into=c("column","SiteRepSpecies"), sep="([\\[\\]])")

myHierarchicalStanResults$SiteRepSpecies <- as.numeric(myHierarchicalStanResults$SiteRepSpecies)
```

```{r}
asv1_occuprob <- occu_df_ASV1 %>% 
  select(Species, Site, SiteSpecies, SiteRepSpecies) %>%
  left_join(myHierarchicalStanResults, by = "SiteRepSpecies") %>% 
  group_by(Species, Site, SiteSpecies) %>%
  summarise(max_Occupancy_prob = max(mean))
asv1_occuprob
```



  
  spread(key = "pcr_replicate", value = "reads", fill = 0, drop = F)

%>%
  group_by(ASV, Sample_ID, reads) %>%
  unique() %>% # eliminate duplicate entries for the same ASV and sample (these seem to be a thing for aquaF2)
  pivot_wider(names_from = Sample_ID, values_from = reads) %>%
  ungroup()













occu_df_temp <- occu_df %>%
  #filter(Species == 1) %>%
  select(Site, K, N) 

occu_df_temp <- occu_df_temp[, -1]   ## get rid of the ASV column that the 'select' function doesn't remove... 

occu_df_temp <- occu_df_temp %>%
  unite(K.N, K, N, sep = ".") 

myrep_oneASV <- c(rep(1:3, length.out = 45), c(1,2), rep(1:3, length.out = 201))   ## site 30 only had 2 biological reps
myrep_full <- rep(myrep_oneASV, times = 1828)

occu_df_temp$myrep <- myrep_full

occu_df_by_site <- occu_df_temp[,-1] %>%
  
  occu_df_by_site <- occu_df_by_site[,-1] %>%
  pivot_wider(names_from = "myrep", values_from = K.N)

occu_df_by_site %>%
  dplyr::group_by(myrep) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)


occu_df_unique <- occu_df_by_site[,-1] %>%  ##remove site column
  unique()







pivot table to have samples by ASV 
```{r}
normalized_pivot <- normalized %>%
  dplyr::select(Sample_ID, sample_type, ASV, Normalized_reads) %>%
  tidyr::pivot_wider(names_from = "ASV", values_from = "Normalized_reads")
```

how many samples, controls are in this data frame? 
  ```{r}
normalized_pivot %>%
  group_by(sample_type) %>%
  summarise(count = n())
```

okay, so there are only environmental samples here. 

i need replicate info so split Sample_ID
```{r}
normalized_pivot <- normalized_pivot %>%
  tidyr::separate(col = "Sample_ID", into = c("extraction_ID", "rep_id", "rep2"), sep = "-", remove = F)
```
the warning is for the samples with both extraction and PCR replicates... for now i'll write some code to treat the extraction reps just as different pcr reps of the same sample.. 

```{r}
normalized_pivot_no_ext_reps <- normalized_pivot[1:636,]
normalized_pivot_ext_reps <- normalized_pivot[637:840,]

normalized_pivot_ext_reps <- normalized_pivot_ext_reps %>%
  unite(rep_id, rep2, sep = "_", col = "rep_id") %>%
  mutate(rep_id = ifelse(rep_id == "1_A", "A", 
                       ifelse(rep_id == "1_B", "B",
                              ifelse(rep_id == "1_C", "C",
                                     ifelse(rep_id == "2_A", "D",
                                            ifelse(rep_id == "2_B", "E",
                                                   ifelse(rep_id == "2_C", "F", rep_id)))))))
  
normalized_pivot_no_ext_reps <- normalized_pivot_no_ext_reps[,-4]

combined <- bind_rows(normalized_pivot_no_ext_reps, normalized_pivot_ext_reps)
```


```{r}
replicates <- combined$rep_id
```

remove the metadata assosiated with the normalized reads 
```{r}
norm_reads <- combined[,-c(1:4)]

#put in zeros for NA is the data frame... double check this appropriate for calculating distance matrix later on. 
norm_reads <- replace(norm_reads, is.na(norm_reads), 0)

row.names(norm_reads) <- combined$Sample_ID
```

calculate distances between samples 
```{r}
distmat <- vegan::vegdist(norm_reads)
```


this code is copied from 20180220_Tides_and_eDNA_RPK.Rmd
```{r}
distList=list(NA); distList.tri=list(NA); index=1
          for (i in unique(replicates)){
          	rowMatch<-which(replicates%in%i)
          	distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
          			distList.tri[[index]]<-	distList[[index]][upper.tri(distList[[index]])]
          	index=index+1
          }

normparams=MASS::fitdistr(unlist(distList.tri), "normal")$estimate  #fit normal distribution to bray-curtis dissimilarities. lognormal, beta, etc, had less-good fits.
probs=pnorm(unlist(distList.tri), normparams[1], normparams[2])
outliers =which(probs>0.80)
minOutlier<-min(unlist(distList.tri)[outliers]) #minimum outlier value
	
 #remove outliers
      #distList<-distList[-which(lapply(distList, length)==1)] 	   ### this line produces an error - come back and make sure this still works as it should. 
      namesOutliers=list(NA)
      for (i in 1:length(distList)){
      	namesOutliers[[i]]<-intersect(
                                names(which(colSums(distList[[i]]>=minOutlier)>0)),
                                names(which.max(rowMeans(distList[[i]])))
                              )
      }
```





```{r}
replicate_outliers <- unlist(namesOutliers)
replicate_outliers
```

there are NO outliers detected with a 0.95 probability threshold... 

okay, using the 80% threshold gives me some samples to look at... 

filter data (need to update some of the object names below to get this code to work...)
```{r}
#keepers_step5 <- combined %>%
#  dplyr::filter(!Sample_ID %in% replicate_outliers)

#ids_keep <- keepers_step5$extraction_ID
```


here i will just filter using the pairwise distances among PCR replicates
```{r}
pcr.distances <- all.distances.to.plot %>%
  filter(Distance.type == "PCR.replicates")

pcr.distances.value <- pcr.distances$value

#normparams <- fitdistr(all_pairwise_distances, "normal")$estimate
normparams <- MASS::fitdistr(pcr.distances.value, "normal")$estimate                                      
#  probs <- pnorm(all_pairwise_distances, normparams[1], normparams[2])
probs <- pnorm(pcr.distances.value, normparams[1], normparams[2])
outliers <- which(probs>0.95)

discard <- pcr.distances[outliers,]

samples_to_discard <- data.frame(unique(c(discard$Sample1, discard$Sample2)))
colnames(samples_to_discard) <- "new_ID"
```

after my first round of filtering, i wonder if the 95% threshold is too stringent... let's plot what pairwise distances were removed
```{r}
pcr.distances$probs <- probs 
pcr.distances

#what is the minimum pairwise distance that is removed by the 95% threshold? 
min_value <- min(pcr.distances$value[pcr.distances$probs > 0.95])

ggplot (pcr.distances) +
  geom_histogram (aes ( x = value, after_stat(ndensity)), position = "dodge",  alpha = 0.9, bins = 50) +
  geom_vline(xintercept = min_value, linetype = "dashed", color = "red") +
  labs (x = "Pairwise dissimilarity", y = "density" ,
        Distance.type = "Distance") +
  guides (fill = "none")

```

okay to convince me if this threshold needs to be moved, let's visualize the samples that are just at this limit...
```{r}
x <- pcr.distances %>%
  filter(value >0.75, value < 0.80) %>%
  arrange(value)
x
```

```{r}
asv_table_filter4 %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 100) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

e02091-1-A and e02091-1-B don't pass the 95% filter... 

```{r}
asv_table_filter4 %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 103) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

e02106-B,C,F; e02107-A,F; e02108-A,D
hmmm.. these don't seem too bad... 

```{r}
asv_table_filter4 %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 37) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

e00479-B,C - again doesn't seem too bad.  

how about samples from site 82 - it's pairwise is just below 0.8 
```{r}
asv_table_filter4 %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 82) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

e02046-A,C  - those look a bit different... 

MAYBE I DON"T WANT TO FILTER BY PAIRWISE DISTANCES... this seems to drop reps like e02046-A that are similar to other reps(e02046-B) just because they are different from another (e02026-C). when just e02026-C should be the one removed..  

okay so the samples from above that are removed are: e00385-B, e00388-B, e00396-C, e00409-B, e00421-C and e00426-A.  This all makes sense visually for these few extraction replicates. cool. 


how similar are these to the other extractions from their locations??? 

```{r}
asv_table_filter4 %>%
  filter(extraction_ID == 'e00388')

asv_table_filter4 %>%
  #filter(!extraction_ID %in% ids_removed) %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 7) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

alright none of these look all that similar.. site7 is bit of a mess... 


```{r}
asv_table_filter4 %>%
  filter(extraction_ID == 'e00426')

asv_table_filter4 %>%
  #filter(!extraction_ID %in% ids_removed) %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 27) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

pcr reps in e00424 and e00425 are more similar to one another than e00426.  i think we will just drop all pcr reps from e00426 (?).


```{r}
asv_table_filter4 %>%
  filter(extraction_ID == 'e00421')

asv_table_filter4 %>%
  #filter(!extraction_ID %in% ids_removed) %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 29) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

the pairwise distance between e00421-A and e00421-C exceeded the threshold...


```{r}
asv_table_filter4 %>%
  filter(extraction_ID == 'e00474')

asv_table_filter4 %>%
  #filter(!extraction_ID %in% ids_removed) %>%
  group_by(Sample_ID) %>%
  mutate(sum=sum(reads)) %>%
  mutate(prop = reads/sum) %>%
  filter(location1 == 36) %>%
  ggplot(aes(x=Sample_ID, y=prop, fill=ASV)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~extraction_ID, scales = 'free', ncol = 3) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95),
    legend.position = "none",
    legend.title = element_blank()
  )  
```

the pairwise distance between e00474-A and e00474-C exceeded the threshold...

okay, i could go on like this forever... i need to make a decision now on how exactly to filter the pcr reps. 




PLOTTING THE NEWLY FILTERED DATA! - remember that i ignored the samples with 2 extraction reps for now. so none of those sites have been filtered for pcr rep dissimilarity. 


