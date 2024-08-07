#!/usr/bin/env Rscript
#### Genotyping with the MykroHap panel ####

## Borrowed some code from Neil for this ##

library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
haps_file <- args[1]
project_name <- args[2]

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

# before we go to wide two-column format we need to convert our haplotypes into numbers
haplo_key <- hap_final %>% 
  select(haplo) %>% 
  distinct(haplo) %>% 
  mutate(n_haplo = 1:nrow(.))

hap_recode <- left_join(hap_final, haplo_key, by = "haplo")

### to convert haplotypes to numeric use this:
# now go wide to a 2 column format to calculate missing data
hap2col_num <- hap_recode %>%
  select(group,indiv.ID,locus,rank,n_haplo) %>%
  unite(tmp, locus:rank,sep = ".") %>%
  spread(tmp,n_haplo)

### to keep the base calls use this instead:
# now go wide to a 2 column format to calculate missing data
hap2col <- hap_final %>%
  select(group,indiv.ID,locus,rank,haplo) %>%
  unite(tmp, locus:rank,sep = ".") %>%
  spread(tmp,haplo)

# calculate missing data per fish
hap_missdata <- 
  hap2col %>%
  mutate(n_miss = rowSums(is.na(.))/2) %>%
  select(indiv.ID, group, n_miss) %>%
  arrange(desc(n_miss))

# how many missing data at more than 5 loci
hap_missdata %>% 
  filter(n_miss > 20) %>%
  nrow(.)

# QED. Let's print it out and do some analyses.
hap2col_output <- stringr::str_c(project_name, "_", format(Sys.Date(), "%Y%m%d"),"_genotypes", ".txt")
write.table(hap2col, file = hap2col_output, sep = "\t", row.names = F, quote = F)

hap2colnum_output <- stringr::str_c(project_name, "_", format(Sys.Date(), "%Y%m%d"),"_numgenotypes", ".txt")
write.table(hap2col_num, file = hap2colnum_output, sep = "\t", row.names = F, quote = F)

missing_data_output <- stringr::str_c(project_name, "_", format(Sys.Date(), "%Y%m%d"), "_missingdata", ".txt")
write.table(hap_missdata, file = hap2colnum_output, sep = "\t", row.names = F, quote = F)