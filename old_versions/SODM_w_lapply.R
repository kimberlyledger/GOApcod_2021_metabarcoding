

### my SODM with a lapply loop 


temporary parameters used during testing
```{r}
#feature_table <- occu_df
#asv_filter_list <- as.list(5)
```

my hierarchical SODM
```{r}
sodm.by.ASV <- function(feature_table, asv_filter_list){
  
  #filter by the species in spp_list
  #temp.df <- feature_table[feature_table$Species == asv_filter, ]
  temp.df <- feature_table[feature_table$Species %in% asv_filter_list, ]
  temp.df$Species_1 <- 1                  ## the stan model doesn't like it when the starting values of species is not 1 so just doing this as work around
  
  #create unique identifier for combinations of site-biologicalrep-ASV; for use in hierarchical modeling
  SDS <- unite(data = temp.df, col = SDS, c("Site", "BiologicalRep", "Species")) %>% pull(SDS)
  temp.df$SiteRepSpecies <- match(SDS, unique(SDS)) #index for unique site-biologicalrep-species combinations
  
  #create unique identifier for combinations of site-ASV; for use in hierarchical modeling
  SS <- unite(data = temp.df, col = SS, c("Site", "Species")) %>% pull(SS)
  temp.df$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations
  
  #####################
  #run Stan model
  #note this will take a while the first time you run a particular model, because it needs to compile from C++
  #####################      
  myHierarchicalModel <- stan(file = "Stan_SOM_hierarchical_with_occuprob.stan", 
                              data = list(
                                S = nrow(temp.df),
                                Species = temp.df$Species_1,
                                Nspecies = length(unique(temp.df$Species_1)),
                                L = temp.df$SiteSpecies,
                                Nloc = length(unique(temp.df$SiteSpecies)),
                                K = temp.df$K,
                                N = temp.df$N,
                                z = ifelse(temp.df$N > 0, 1, 0)
                              ), 
                              chains = 4,   #number of chains
                              iter = 4000   #number of iterations per chain
  )
  
  myHierarchicalStanResults <- tidy(tibble(as.data.frame(myHierarchicalModel)))   
  #plot(myHierarchicalModel) 
  
  #extract the information I want from the stan results and put into a table 
  
  ## occupancy probabilities 
  myHierarchicalStanResults_occu <- myHierarchicalStanResults %>%
    filter(grepl("Occupancy_prob", column)) %>%
    separate(column, into=c("column","SiteRepSpecies"), sep="([\\[\\]])")
  
  myHierarchicalStanResults_occu$SiteRepSpecies <- as.numeric(myHierarchicalStanResults_occu$SiteRepSpecies)
  
  occupancy_prob <- temp.df %>% 
    select(Species, Site, SiteSpecies, SiteRepSpecies) %>%
    left_join(myHierarchicalStanResults_occu, by = "SiteRepSpecies") %>% 
    group_by(Species, Site, SiteSpecies) %>%
    summarise(max_Occupancy_prob = max(mean))
  
  ## true positive 
  myHierarchicalStanResults_p11 <- myHierarchicalStanResults %>%
    filter(grepl("p11", column)) %>%
    separate(column, into=c("column","Species_1"), sep="([\\[\\]])") %>%
    select(column, mean, sd)
  
  myHierarchicalStanResults_p11$Species <- asv_filter_list
  
  ## false positive 
  myHierarchicalStanResults_p10 <- myHierarchicalStanResults %>%
    filter(grepl("p10", column)) %>%
    separate(column, into=c("column","Species_1"), sep="([\\[\\]])") %>%
    select(column, mean, sd)
  
  myHierarchicalStanResults_p10$Species <- asv_filter_list
  
  return(list(occupancy_prob, myHierarchicalStanResults_p11, myHierarchicalStanResults_p10))
}
```


my short species list used for testing out the function
```{r}
#asv_filter_list <- list(1,2,3)
```

make list of the ASVs (now in my occu_df as "Species") that can be cycled over 
```{r}
asv_filter_list <- occu_df$Species %>%
  unique() %>%
  as.list()
```

cycle over the ASVs using lapply and generate csv outputs that I can then filter based on occupancy probabilties
```{r}
output_list <- lapply(asv_filter_list, sodm.by.ASV, feature_table = occu_df)
```

combine outputs from each ASV (remember here it's called 'Species' into single object, but keep occupancy prob, p11, and p10 outputs separate
```{r}
my_occupancy_list <- lapply(output_list, '[[', 1)
my_p11_list <- lapply(output_list, '[[', 2)
my_p10_list <- lapply(output_list, '[[', 3)
```

now take them each out of the list format
```{r}
my_occupancy <- do.call(rbind, my_occupancy_list)
my_p11 <- do.call(rbind, my_p11_list)
my_p10 <- do.call(rbind, my_p10_list)
```


plot occupancy probabilities
```{r}
asv1_occuprob %>% 
  ggplot() + geom_histogram(aes(max_Occupancy_prob)) +
  ggtitle(label = "Histogram of Max Occupancy Probabilities Per Biological Replicate") +
  ylab("Count") + xlab("Occupancy Probability")
```


my feeling is that i should be using the second stan model and filter by the max occupancy probability calculated for each SiteSpecies (aka a biological replicate)

```{r}
mythreshold <- 0.8 

asvs_to_remove_bySite <- asv1_occuprob %>%
  filter(max_Occupancy_prob < mythreshold) %>%
  mutate(ASV = paste0("ASV", Species)) %>%
  select(ASV, Site)
```



