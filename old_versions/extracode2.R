my make-shift attempt at occupancy modeling..... 

**need to come back and figure out how to do this for real...**
  
  right now this code retains only ASVs present in > or = 0.3 of the site replicates (most sites have 3 bio reps with 3 pcr reps)
  
  ```{r}
  sites <- metadata %>%
    dplyr::select(Sample_ID, location1)
  
  rep_table <- asv_table_filter3 %>%
    #filter(ASV == "ASV1") %>%
    left_join(sites, by = "Sample_ID") %>%
    arrange(location1)
  
  # okay, by doing the left_join i reintroduced the NA sites.... 
  rep_table <- rep_table %>%
    filter(location1 != "NA")
  
  ### i used this next section of code with just ASV1 to determine the replicate pattern 
  rep_counts <- data.frame(table(rep_table$location1))
  sum(rep_counts$Freq) #840
  
  #there has to be a better way to code this but for now.... 
  myreps <- c(rep(1:9, 4), rep(1:7, 1), rep(1:9, 10), rep(1:6, 1), rep(1:9, 3), rep(1:8, 1), rep(1:9, 28), rep(1:8, 3), rep(1:9, 20), rep(1:12, 1), rep(1:18, 11))
  
  #####
  full_reps <- rep(myreps, 1837)
  
  #now i need a "visit" variable... the data frame must be arranged by AVS first 
  rep_table_byASV <- rep_table %>%
    arrange(ASV)
  
  #rep_table_byASV$replicate <- myreps
  rep_table_byASV$replicate <- full_reps
  
  rep_table_byASV <- rep_table_byASV[,-c(1:2)] %>%
    tidyr::pivot_wider(names_from = "replicate", values_from = "reads")
  
  rep_table_byASV_1 <- rep_table_byASV[,c(1:2)]
  rep_table_byASV_2 <- rep_table_byASV[,c(3:20)]
  
  rep_table_byASV_2[rep_table_byASV_2 > 0] <- 1
  
  occupancy <- rep_table_byASV_2 %>%
    mutate(numb_occupied = rowSums(., na.rm = T),
           n_replicates = rowSums(!is.na(.)),
           prop_occupied = numb_occupied/n_replicates)
  
  occupancy_sites <- cbind(rep_table_byASV_1, occupancy[,21])
  
  occu_df <- occupancy_sites %>%
    group_by(ASV) %>%
    summarise(max_prop = max(prop_occupied))
  
  occu_asv_filter <- occu_df %>%
    filter(max_prop >= 0.3)       #### this values is important! 
  
  asvs_to_keep <- occu_asv_filter$ASV
  ```
  
  filter the data frame 
  ```{r}
  asv_table_filter4 <- asv_table_filter3 %>%
    filter(ASV %in% asvs_to_keep)
  ```
  
  what ASVs did this toss out??? 
    ```{r}
  asv_removed2 <- asv_table_filter3 %>%
    filter(!ASV %in% asvs_to_keep)
  ```
  
  ASVs removed from samples
  ```{r}
  asv_removed2 %>%
    filter(sample_type == "sample") %>%
    ggplot(aes(x=ASV, y=reads, fill=Sample_ID)) +
    geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(
      axis.text.x = element_blank(),
      legend.position = "none",
      legend.title = element_blank()
    )  
  ```
  
  hmmm... doesn't seem like i should get rid of all these. real occupancy modeling might help here.or maybe run occupancy modeling at just the sample level (not the site).  or maybe save this even until after taxonomic assignment?? 

```{r}
asv_removed2 %>%
  filter(sample_type == "sample") %>% 
  group_by(ASV) %>%
  summarize(total_reads = sum(reads)) %>%
  arrange(desc(total_reads)) %>%
  head()
```

ASV85 = Gadidae
ASV88 = Gadidae
ASV91 = Stichaeus punctatus
ASV95 = Gadus chalcogrammus
etc....  

**this step is when lots of ASVs are dropped from the data set** - not necessarily a bad thing but need to check this
