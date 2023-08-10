asv cleaning - GOApcod2021
================
Kimberly Ledger
2023-05-23

this script “decontaminates” the ASVs and sample reads from the mifish
libraries (May 9 2023 and May 12 2023) that inclded samples from the GOA
pcod survey in 2021

i will try to follow the eDNA decontamination pipeline suggested here:
<https://github.com/zjgold/gruinard_decon>. But make some of my own
adaptations along the way.

load libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
    ## ✔ tibble  3.2.1      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(fitdistrplus) #for fitdist()
```

    ## Loading required package: MASS
    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## Loading required package: survival

``` r
#library(microDecon) #for decontamination based on negatives 
```

the data i am starting with is the ASV table output from DADA2 pipeline
in my sequence_filtering.Rmd in the eDNA_metabarcoding folder

``` r
asv_table <- read.csv("/genetics/edna/workdir/GOApcod_2021/combined/trimmed/filtered/outputs/ASVtable.csv") %>%
  rename(Sample_ID = X)

asv_table$Sample_ID <- as.factor(asv_table$Sample_ID)
```

# 1. Remove ASVs present less than 1% of samples (1% of 1014 samples is about 10)

``` r
rownames(asv_table) <- asv_table$Sample_ID

asv_present <- asv_table[-1]
asv_present[asv_present > 0] <- 1
asv_count <- data.frame(colSums(asv_present)) %>%
  filter(colSums.asv_present. > 10)
asv_count_filter <- rownames(asv_count)

asv_table_filter1 <- asv_table[, asv_count_filter]
```

# 2. Estimate index hopping

subtract the proportion of reads that jumped into the control samples
from each environmental sample

we need sample metadata to do this…

``` r
metadata <- read.csv("/genetics/edna/workdir/GOApcod_2021/GOA2021_metadata_20230515.csv")

#illumina output changed "_" to "-"
metadata$Sample_ID <- gsub("_", "-", metadata$Sample_ID) 
```

add column to asv table that labels the sample type

``` r
asv_temp <- asv_table_filter1
asv_temp$Sample_ID <- asv_table$Sample_ID

samp_type <- metadata %>%
  dplyr::select(Sample_ID, sample_type) %>%
  left_join(asv_temp, by = "Sample_ID")
```

identify maximum proportion of reads for each ASV found in the positive
controls

``` r
pos_asvs <- samp_type %>%
  filter(sample_type == "positive_control") %>%
  pivot_longer(cols = c(3:143), names_to = "ASV", values_to = "reads") %>%
  group_by(Sample_ID) %>%
  mutate(TotalReadsPerSample = sum(reads)) %>%
  mutate(Prop = reads/TotalReadsPerSample) %>%
  group_by(ASV) %>%
  summarise(max_prop = max(Prop))
```

now subtract this max proportion for each ASV from environmental samples
and all negative controls

``` r
indexhop_table <- samp_type %>%
  filter(sample_type != "positive_control") %>% ## working all samples except the positive controls 
  pivot_longer(cols = c(3:143), names_to = "ASV", values_to = "reads") %>%
  group_by(Sample_ID) %>%
  mutate(TotalReadsPerSample = sum(reads)) %>%
  left_join(pos_asvs, by = "ASV") %>%
  mutate(IndexHoppingReads = TotalReadsPerSample*max_prop) %>%
  mutate(reads_IndexHop_removed = reads - IndexHoppingReads) %>%
  mutate(reads_IndexHop_removed = if_else(reads_IndexHop_removed < 0, 0, reads_IndexHop_removed))
```

clean up the table by removing columns no longer needed

``` r
asv_table_filter2 <- indexhop_table %>%
  dplyr::select(Sample_ID, sample_type, ASV, reads_IndexHop_removed) %>%
  rename(reads = reads_IndexHop_removed)
```

# 3. Discard PCR replicates with low numbers of reads

calculate reads per sample - just consider the field blanks and
environmental samples here

``` r
all_reads <- asv_table_filter2 %>%
  filter(sample_type == "field_blank" | sample_type == "sample") %>%
  group_by(Sample_ID) %>%
  summarize(ReadsPerSample = sum(reads))
```

fit a normal distribution

``` r
fit <- fitdist(all_reads$ReadsPerSample, "gamma", lower=c(0,0), start=list(scale=1,shape=1))

all_reads %>%  
  mutate(prob = pgamma(ReadsPerSample, shape = fit$estimate[[2]], scale = fit$estimate[[1]], lower.tail = TRUE,
       log.p = FALSE)) -> all_reads
```

identify and remove the outliers

``` r
low_dist_probability_cutoff <- 0.05
minimum_read_cutoff <- 1000

outliers <- all_reads %>% 
  filter(prob < low_dist_probability_cutoff  | ReadsPerSample < minimum_read_cutoff) # changed to 0.05 to save the two samples
outlierIDs <- outliers$Sample_ID
```

note: the minimum read cutoff was not actually necessary because all
samples above 0.05 also had at least 1000 reads

``` r
asv_table_filter3 <- asv_table_filter2 %>%
  filter(!Sample_ID %in% outlierIDs)
```

# 4. Account for contaminants in extraction and pcr negative controls

will ignore field negatives for now…

try out the microDecon package. microDecon compares the prevalence of
ASVs in blanks to the prevalence in samples and removes contaminant ASVs
and subtracts contaminant sequences.

the data frame must be a column of OTU IDs, at least one column from a
blank sample, at least one column from actual samples

``` r
# asv_table_filter3$sample_type <- as.factor(asv_table_filter3$sample_type)
# 
# for_decon <- asv_table_filter3 %>%
#   mutate(sample_type = factor(sample_type, levels = c("extraction_blank", "PCR_blank", "field_blank", "sample"))) %>%
#   arrange(sample_type) %>%
#   dplyr::select(!sample_type) %>%
#   pivot_wider(names_from = "Sample_ID", values_from = "reads")
```

``` r
#decontaminated <- decon(data = for_decon, numb.blanks = 39, numb.ind = 918, taxa = FALSE, regression = 2)
```

well, i only get error messages for different variations of this….

since the decon function is not running for me, i’ll try out a similar
approach by hand…

``` r
temp <- asv_table_filter3 %>%
  filter(sample_type == "extraction_blank" | sample_type == "PCR_blank") %>%
  filter(Sample_ID != "e00562-A") %>%  # remove because something funky happened with this one 
  group_by(Sample_ID) %>%
  mutate(TotalReadsPerSample = sum(reads)) %>%
  mutate(Prop = reads/TotalReadsPerSample) %>%
  group_by(ASV) %>%
  summarise(max_prop = max(Prop, na.rm = T),
            quant3_prop = quantile(Prop, probs = 0.75, na.rm = T))
```

biggest contamination issue is with ASV1 (which is pcod) ASV73 = Homo
sapien so we’ll remove all of this one any ways… ASV3 = Oncorhynchus
gorbuscha ASV2 = Clupea pallasii

clearly there needs to be some decontamination that takes place here but
i’m not exactly sure what that should be…

maybe skipping to the next step of using site occupancy modeling will
help solve this issue? or just setting a read count threshold for every
ASV to subtract (i.e. not using proportions)???

``` r
nc_reads <- asv_table_filter3 %>%
  filter(sample_type == "extraction_blank" | sample_type == "PCR_blank") %>%
  filter(Sample_ID != "e00562-A") %>%  # remove because something funky happened with this one 
  group_by(ASV) %>%
  summarise(max_reads = max(reads),
            quant3_reads = quantile(reads, probs = 0.75, na.rm = T))
```

subtract read numbers from samples

``` r
negcontcontam_table <- asv_table_filter3 %>%
  left_join(nc_reads, by = "ASV") %>%
  mutate(reads_negcontcontam_removed = reads - max_reads) %>%
  mutate(reads_negcontcontam_removed = if_else(reads_negcontcontam_removed < 0, 0, reads_negcontcontam_removed))
```

``` r
asv_table_filter4 <- negcontcontam_table %>%
  dplyr::select(Sample_ID, sample_type, ASV, reads_negcontcontam_removed) %>%
  rename(reads = reads_negcontcontam_removed)
```

# 5. Site Occupancy Modeling…

well i haven’t been able to crack how to do this yet… but for now i will
calculate the % occurrence for each ASV for each site. if an ASV does
not occur in at least 50% of replicates from at least one site, it is
removed.

for each ASV, i will create a pres-abs matrix. each matrix will be 83
sites by 9 visits (aka replicates). well actually there are more than 9
“visits” for some sites.

the data frame can have missing observation records (NAs). some sites
have more replicates than others and that’s okay.

there are 141 ASVs currently in the data set… which is too many to try
and check this by hand…

``` r
sites <- metadata %>%
  dplyr::select(Sample_ID, location1)

rep_table <- asv_table_filter4 %>%
  #filter(ASV == "ASV1") %>%
  left_join(sites, by = "Sample_ID") %>%
  arrange(location1)

# okay, by doing the left_join i reintroduced the NA sites.... 
rep_table <- rep_table[1:115197,]

### i used this next section of code with just ASV1 to determine the replicate pattern 
rep_counts <- data.frame(table(rep_table$location1))
#sum(rep_counts$Freq) #817 

#there has to be a better way to code this but for now.... 
myreps <- c(rep(1:9, 4), rep(1:7, 1), rep(1:9, 8), rep(1:6, 1), rep(1:9, 1), rep(1:6, 1), rep(1:9, 3), rep(1:8, 1), rep(1:9, 22), rep(1:5, 1), rep(1:9, 2), rep(1:3, 1), rep(1:9, 2), rep(1:7, 1), rep(1:8, 2), rep(1:9, 5), rep(1:7, 1), rep(1:9, 1), rep(1:8, 1), rep(1:9, 7), rep(1:8, 1), rep(1:9, 4), rep(1:11, 1), rep(1:18, 1), rep(1:14, 1), rep(1:18, 9))

#####
full_reps <- rep(myreps, 141)

#now i need a "visit" variable... the data frame must be arranged by AVS first 
rep_table_byASV <- rep_table %>%
  arrange(ASV)

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
  filter(max_prop >= 0.5)

asvs_to_keep <- occu_asv_filter$ASV
```

``` r
asv_table_filter5 <- asv_table_filter4 %>%
  filter(ASV %in% asvs_to_keep)
```

# 6. Dissimilarity between PCR (biological) replicates

This step removes samples for which the dissimilarity between PCR
replicates exceeds the normal distribution of dissimilarities observed
in samples. The objective of this step is to remove any technical
replicates that look like they do not belong.

following the gruinard_decon again, i’ll calculate an eDNA index

``` r
normalized <- asv_table_filter5 %>%
  group_by(Sample_ID) %>%
  mutate(Tot = sum(reads),
         Prop_reads = reads/Tot) %>%
  dplyr::group_by(ASV) %>%
  mutate(Colmax = max(Prop_reads, na.rm = TRUE),
         Normalized_reads = Prop_reads/Colmax)
```

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

pivot table to have samples by ASV

``` r
normalized_pivot <- normalized %>%
  dplyr::select(Sample_ID, sample_type, ASV, Normalized_reads) %>%
  tidyr::pivot_wider(names_from = "ASV", values_from = "Normalized_reads")
```

i need replicate info so split Sample_ID

``` r
temp <- normalized_pivot %>%
  dplyr::filter(sample_type == "sample") %>%
  tidyr::separate(col = "Sample_ID", into = c("ID", "rep"), sep = "-", remove = F)
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 224 rows [690, 691,
    ## 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707,
    ## 708, 709, ...].

the warning is for the samples with both extraction and PCR replicates…
for now i’ll just treat all replicates of those as PCR reps. but should
probably deal with this later on.

``` r
replicates <- temp$ID
```

``` r
norm_reads <- temp[,-c(1:4)]

row.names(norm_reads) <- temp$Sample_ID
```

    ## Warning: Setting row names on a tibble is deprecated.

calculate distances between samples

``` r
distmat <- vegan::vegdist(norm_reads)
```

this code is copied from 20180220_Tides_and_eDNA_RPK.Rmd

``` r
distList=list(NA); distList.tri=list(NA); index=1
          for (i in unique(replicates)){
            rowMatch<-which(replicates%in%i)
            distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
                    distList.tri[[index]]<- distList[[index]][upper.tri(distList[[index]])]
            index=index+1
          }

normparams=MASS::fitdistr(unlist(distList.tri), "normal")$estimate  #fit normal distribution to bray-curtis dissimilarities. lognormal, beta, etc, had less-good fits.
probs=pnorm(unlist(distList.tri), normparams[1], normparams[2])
outliers =which(probs>0.95)
minOutlier<-min(unlist(distList.tri)[outliers]) #minimum outlier value
    
 #remove outliers
      distList<-distList[-which(lapply(distList, length)==1)]   
      namesOutliers=list(NA)
      for (i in 1:length(distList)){
        namesOutliers[[i]]<-intersect(
                                names(which(colSums(distList[[i]]>=minOutlier)>0)),
                                names(which.max(rowMeans(distList[[i]])))
                              )
      }
```

``` r
replicate_outliers <- unlist(namesOutliers)

decontam <- temp %>%
  dplyr::filter(!Sample_ID %in% replicate_outliers)
```

okay so the ASVs and the samples in ‘decontam’ have made it through the
decontamination methods…

and this table will have the read numbers (not normalized)

``` r
asv_table_filter6 <- asv_table_filter5 %>%
    dplyr::filter(!Sample_ID %in% replicate_outliers)
```

``` r
write.csv(asv_table_filter6, "decontaminated_reads.csv")
```
