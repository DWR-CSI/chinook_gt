#!/usr/bin/env Rscript
#### Genotyping with the MykroHap panel ####

## Borrowed some code from Neil for this ##

library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
haps_file <- args[1]
project_name <- args[2]
locus_indices <- args[3] # use latest version
na_string_var = "ND"

# Read in the unfiltered raw haplotype data
hap_raw <- read_csv(haps_file)
# Drop the first column if the first column is an index
if (colnames(hap_raw)[1] %in% c("X1","")) {
  hap_raw <- hap_raw %>% select(-1)
}

# We'll define a run name for subsequent use
run_name <- project_name

# Let's make sure we don't have a bunch of haplotypes with one or multiple
# Ns at the beginning and end. We will trim those off first and see what we're left with
hap_raw$haplo <- str_replace(hap_raw$haplo, "^N+", "")
hap_raw$haplo <- str_replace(hap_raw$haplo, "N+$", "")


# And now any Ns that remain are in the middle of haplotypes and we don't want those
# although most have low read depth and will be tossed anyway
# We will leave the haplotypes with Xs as most of them also have very low read depth and will
# be filtered later, the real ones are quality variants (presumably with Mendelian inheritance)
hap2remove <- hap_raw %>% filter(str_detect(haplo, "N")) %>% pull(haplo)

# Now compile problematic 3 hap individuals and loci
three_hap_fish <- hap_raw %>% 
  filter(!haplo %in% hap2remove, depth > 4) %>% 
  arrange(locus, indiv.ID, rank) %>%
  group_by(indiv.ID,locus) %>%
  mutate(rank = row_number()) %>%
  mutate(allele.balance = depth/depth[1]) %>%
  filter(allele.balance > 0.2999) %>% 
  filter(rank > 2) %>%
  mutate(tmp_id = paste0(indiv.ID,locus)) %>% distinct(tmp_id) %>% pull(tmp_id)

three_hap2chuck <- hap_raw %>% 
  filter(!haplo %in% hap2remove, depth > 4) %>% 
  arrange(locus, indiv.ID, rank) %>%
  group_by(indiv.ID,locus) %>%
  mutate(rank = row_number()) %>%
  mutate(allele.balance = depth/depth[1]) %>%
  filter(allele.balance > 0.2999) %>%
  mutate(tmp_id = paste0(indiv.ID,locus)) %>%
  filter(tmp_id %in% three_hap_fish) %>%
  mutate(ab_test = depth/depth[2]) %>% 
  filter(ab_test > 0.5) %>%
  filter(any(rank > 2)) %>%
  ungroup() %>%
  distinct(tmp_id) %>% pull(tmp_id)

three_hap2flag <- hap_raw %>% 
  filter(!haplo %in% hap2remove, depth > 4) %>% 
  arrange(locus, indiv.ID, rank) %>%
  group_by(indiv.ID,locus) %>%
  mutate(rank = row_number()) %>%
  mutate(allele.balance = depth/depth[1]) %>%
  filter(allele.balance > 0.2999) %>%
  mutate(tmp_id = paste0(indiv.ID,locus)) %>%
  filter(tmp_id %in% three_hap_fish) %>%
  mutate(ab_test = depth/depth[2]) %>% 
  filter(ab_test > 0.5) %>%
  filter(!any(rank > 2)) %>% distinct(tmp_id, .keep_all = TRUE) %>%
  select(tmp_id, group) %>%
  mutate(three_hap = "true", gtseq_run = run_name) %>%
  select(gtseq_run, group, indiv.ID, locus, three_hap)

three_hap2track <- hap_raw %>% 
  filter(!haplo %in% hap2remove, depth > 4) %>% 
  arrange(locus, indiv.ID, rank) %>%
  group_by(indiv.ID,locus) %>%
  mutate(rank = row_number()) %>%
  mutate(allele.balance = depth/depth[1]) %>%
  filter(allele.balance > 0.2999) %>% 
  filter(rank > 2) %>% ungroup() %>% 
  distinct(locus) %>%
  mutate(three_hap = "true", gtseq_run = run_name) %>% select(gtseq_run, locus, three_hap)


# Next filter step is to remove all hap 3 and higher, then remove the three_hap2chuck
# Also enact the total read depth of 20 Diana filter at this stage.

hap_int <- hap_raw %>% 
  filter(!haplo %in% hap2remove, depth > 4) %>% 
  arrange(locus, indiv.ID, rank) %>%
  group_by(indiv.ID,locus) %>%
  mutate(rank = row_number()) %>%
  mutate(allele.balance = depth/depth[1]) %>%
  filter(allele.balance > 0.2999, rank < 3) %>% #removes all 3rd ranked  hap and higher
  mutate(tmp_id = paste0(indiv.ID,locus)) %>%
  filter(!tmp_id %in% three_hap2chuck) %>% #removes unresolvable three hap occurences
  mutate(tot_depth = sum(depth)) %>%
  filter(tot_depth > 19) %>% # retains total depth of 20 or higher.
  ungroup() 

# Need to identify the homozygotes so we can duplicate
homoz2dup <- hap_int %>% 
  group_by(indiv.ID,locus) %>%
  count() %>% filter(n < 2) %>%
  mutate(tmp_id = paste0(indiv.ID,locus)) %>% pull(tmp_id)

# for the homozygotes we need to add the 2nd haplo to the dataframe before spreading it
homoz2add <- hap_int %>% 
  mutate(tmp_id = paste0(indiv.ID,locus)) %>%
  filter(tmp_id %in% homoz2dup) %>%
  mutate(rank = 2) %>%
  mutate(allele.balance = NA)

# bind rows so that each fish has 2 haplo per locus
hap_final <- bind_rows(hap_int,homoz2add) %>%
  mutate(run = run_name) %>%
  select(-tmp_id) %>%
  select(run,everything())

# that should be a long haplotype-frame
# discount double-check
hap_final %>% 
  group_by(locus) %>%
  count(indiv.ID) %>% 
  ungroup() %>%
  distinct(n)

# Keep and export haplotypes
# now go wide to a 2 column format to calculate missing data
hap2col <- hap_final %>%
  select(group,indiv.ID,locus,rank,haplo) %>%
  unite(tmp, locus:rank,sep = ".") %>%
  spread(tmp,haplo)
hap2col[is.na(hap2col)]  <- na_string_var
hap2col_output <- stringr::str_c(project_name, "_", format(Sys.Date(), "%Y%m%d"),"_genotypes", ".txt")
write_tsv(hap2col, file = hap2col_output)
# calculate missing data per fish
hap_missdata <- 
  hap2col %>%
  select(-"group") %>%
  mutate(n_miss = rowSums(. == na_string_var)/2) %>%
  select(indiv.ID, n_miss) %>%
  arrange(desc(n_miss))

# how many missing data at more than 5 loci
hap_missdata %>% 
  filter(n_miss > 20) %>%
  nrow(.)

# QED. Let's print it out and do some analyses.

missing_data_output <- stringr::str_c(project_name, "_", format(Sys.Date(), "%Y%m%d"), "_missingdata", ".txt")
write_tsv(hap_missdata, file = missing_data_output)

#############################################
###  Begin Arthur's rescoring code here   ###
#############################################

baseline_data <- hap2col
locus_index <- read_csv(locus_indices)

#look for new alleles for locus
pivot_data<-baseline_data%>%
  pivot_longer(cols=3:ncol(baseline_data),
               names_to = 'locus',
               values_to = 'allele')
new_data<-pivot_data %>%
  mutate(locus = substr(locus, 1, nchar(locus)-2)) %>%
  select(locus, allele) %>%
  distinct() %>%
  drop_na(allele) %>%
  filter(allele != '')



new_alleles<-new_data %>% 
  anti_join(locus_index) %>%
  mutate(index = NA)

#append new alleles
new_indices<-locus_index %>%
  rbind(new_alleles)


locus_list<-unique(new_alleles$locus)
for_data<-data.frame()
for(i in 1:length(locus_list)){
  l<-locus_list[i]
  j<-locus_index%>%filter(locus==l)
  m<-max(j$index)
  d<-new_indices%>%
    filter(locus==l)
  for(n in 1:nrow(d)){
    d$index[n]<-ifelse(is.na(d$index[n]),max(d$index,na.rm=T)+1,d$index[n])
  }
  for_data<-for_data%>%rbind(d)
}
new_indices<-new_indices%>%
  filter(!locus%in%locus_list)
new_indices_long<-new_indices%>%
  rbind(for_data)

#get wide locus indices
new_indices_wide<-new_indices_long%>%
  pivot_wider(names_from =index,values_from = allele)

#get fish locus index
fish_data<-pivot_data %>%
  mutate(locus = substr(locus, 1, nchar(locus)-2)) %>%
  inner_join(new_indices)


########## Jeff's part for spreading the recoded genos

### need to change "haplo" to "index" in hap_raw's equivalent (fish_data)

#get wide fish indices
fish_data_wide<-fish_data %>%
  pivot_wider(names_from =locus,values_from = index)

hap_int <- fish_data %>% 
  arrange(locus, indiv.ID) %>%
  group_by(indiv.ID,locus) %>%
  mutate(rank = row_number()) %>%
  mutate(tmp_id = paste0(indiv.ID,locus)) %>%
  ungroup() 

# Need to identify the homozygotes so we can duplicate
homoz2dup <- hap_int %>% 
  group_by(indiv.ID,locus) %>%
  count() %>% filter(n < 2) %>%
  mutate(tmp_id = paste0(indiv.ID,locus)) %>% pull(tmp_id)

# for the homozygotes we need to add the 2nd haplo to the dataframe before spreading it
homoz2add <- hap_int %>% 
  mutate(tmp_id = paste0(indiv.ID,locus)) %>%
  filter(tmp_id %in% homoz2dup) %>%
  mutate(rank = 2) %>%
  mutate(allele.balance = NA)

# bind rows so that each fish has 2 haplo per locus
hap_final <- bind_rows(hap_int,homoz2add) %>%
  mutate(run = run_name) %>%
  select(-tmp_id) %>%
  select(run,everything())

# that should be a long haplotype-frame
# discount double-check
hap_final %>% 
  group_by(locus) %>%
  count(indiv.ID) %>% 
  ungroup() %>%
  distinct(n)

### to keep the base calls use this instead:
# now go wide to a 2 column format
hap2col_num <- hap_final %>%
  select(group,indiv.ID,locus,rank,index) %>%
  unite(tmp, locus:rank,sep = ".") %>%
  spread(tmp,index)

hap2col_num[is.na(hap2col_num)]  <- as.numeric("-9") #a last tidbit of cleanup for the missing data code

# results files - note that the code adds blank spaces to the filenames at times, still needs to be fixed
# makes a new index file with any new alleles found in this run
# makes a recoded fish file in two column format

hap2col_num_output <- stringr::str_c(project_name, "_", format(Sys.Date(), "%Y%m%d"),"_numgenotypes.txt")
write_csv(hap2col_num, file = hap2col_num_output, na = "#N/A")
new_indices_filename <- stringr::str_c(format(Sys.Date(), "%Y%m%d"), "_", project_name, "_locus_indices.csv")
write_csv(new_indices_long, file = new_indices_filename, quote = "all")