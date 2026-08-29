#!/usr/bin/env Rscript



###############################################################################
#                                                                             #
#                           Required Libraries                                #
#                                                                             #
###############################################################################


library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(dplyr)
library(argparser)

###############################################################################
#                                                                             #
#                               Functions                                     #
#                                                                             #
###############################################################################



#'Locate the last instance of a pattern in a string.
#' 
#'Finds the starting position of the last instance of a pattern within a string.
#'This function is a wrapper for `stringr::str_locate_all` so pattern matching behavior follows that of `stringr::str_locate_all`
#'
#'@param string the string to be searched.
#'@param pattern the pattern to search for.
#'@returns An integer indicating the start position of the last instance of the pattern.
str_locate_last <- function(string, pattern){
  x <- str_locate_all(string, pattern)[[1]]
  last <- dim(x)[1]
  return(as.numeric(x[last,1]))

  }


#' Read data from an allele index file.
#'
#' Reads in the data from an allele index file. The file format is CSV with columns `locus`, `index` and `allele`.
#' The index may contain alleles with a string representing uncalled alleles.
#' 
#' @param index_file the name of the index file to load.
#' @param missing_string string used to represent uncalled alleles in the index.
#' @returns An object of class `index`.
#' An object of class `index` is a list containing the following components:
#' `original` The index as represented in the input file.
#' `clean` The index with representations of uncalled alleles removed
read_index <- function(index_file, missing_string = "ND"){
  out <- list()
  out$original <- read_csv(index_file)
  out$clean <- subset(out$original, allele != missing_string)
  class(out) <- "index"
  return(out)
}

#' Read data from a wide-format table of haplotypes.
#' 
#' It is assumed that in wide format, loci will be represented as pairs of columns with a separator and an allele identifier at the end of the string.
#' 
#' @param haplotypes_file the file to read.
#' @param allele_sep the separator between locus names and allele specifiers.
#' @param na_string the string used to identify missing data. Will be converted to `NA`.
#' @param sample_col the index of the colum n in the data that contains individual identifiers .
#' @param drop_cols indices of any columns that do not contain either the individual names or haplotype data. Dropped columns are saved in the `dropped_cols` element of the returned object.
#' @return An object of class `haplotype` contianing the following components:
#' `loci` Names of the loci.
#' `allele_slots` The string used to identify allele slots within loci.
#' `allele_sep` The string used to separate locus names from allele slots
#' `dropped_cols` Any columns dropped from the input.
#' `wide` A table of haplotypes in wide format.
#' `long` A table of haplotypes in long format.
read_haplotypes_wide <- function(haplotypes_file,
                                 allele_sep = ".",
                                 na_string = NULL,
                                 sample_col = 2,
                                 drop_cols = NULL){
  wide <- read_csv(haplotypes_file)
  # Set missing to NA
  if (!is.null(na_string)){
    wide[wide == na_string] <- NA}
  sample_col_name <- names(wide)[sample_col]
  # remove any columns to drop
  if (!is.null(drop_cols)){
    dropped_cols <- wide[,drop_cols]
    wide <- wide[, !(names(wide) %in% names(dropped_cols))]} 
  # Ensure sample col is the first col, rest are assumed to be data columns
  data_cols <- names(wide)[names(wide) != sample_col_name]
  wide <- wide[c(sample_col_name, data_cols)]
  #convert to long form
  long <- pivot_longer(wide, all_of(data_cols))
  names(long) <- c(sample_col_name,
                   "locus",
                   "allele")
  #strip allele specifiers from locus names and save as `allele_slot`
  allele_sep_pos <- unlist(sapply(long$locus, str_locate_last, fixed(allele_sep)))
  long$allele_slot <- map2_chr(long$locus, allele_sep_pos + 1, str_sub)
  long$locus <- map2_chr(long$locus, allele_sep_pos - 1, function(x, y){str_sub(x, end = y)})
  #assemble the returned haplotype object
  hap <- list(loci = unique(long$locus),
              allele_slots = unique(long$allele_slot),
              allele_sep = allele_sep,
              dropped_cols = dropped_cols,
              wide = wide,
              long = long)
  class(hap) <- "haplotype"
  return(hap)
}


#' Search for new haplotypes and assign indices
#' 
#' New haplotypes are given new index values starting from the next available index for the locus.
#' @param index the index to update
#' @param haplotypes the haplotype object to search for new haplotypes
#' @return an object of class `index`
update_index <- function(index, haplotypes){
  stopifnot(class(index) == "index", class(haplotypes) == "haplotype")
  #check if any haplotypes in the haplotypes file are not represented in the index.
  distinct_calls <- arrange(unique(select(filter(haplotypes$long, 
                                                 !is.na(allele)), 
                                          c("locus", "allele"))), 
                            locus, allele)
  new_calls <- distinct_calls[!(paste(distinct_calls$locus, 
                                      distinct_calls$allele, 
                                      sep = "_") %in% 
                                  paste(index$clean$locus, 
                                        index$clean$allele, 
                                        sep = "_")),]
  #
  # If no new haplotypes, nothing to do
  #
  if(dim(new_calls)[1] == 0){
    print("No new haplotypes found")
    return(index)}
  #
  # Found at least one new haplotype, update the index
  #
  loci_to_update <- idx$clean[idx$clean$locus %in% new_calls$locus,]
  new_index_starts <- summarize(idx$clean[idx$clean$locus %in% new_calls$locus,], .by = locus, start = max(index) + 1)
  new_calls$index <- NA
  for ( l in unique(new_calls$locus)){
    locus_new_haplotype_count <- dim(subset(new_calls, locus == l))[1]
    first_new_index <- subset(new_index_starts, locus == l)$start
    last_new_index <- first_new_index + (locus_new_haplotype_count - 1)
    new_indices <- seq(first_new_index, last_new_index)
    new_calls[new_calls$locus == l,]$index <- new_indices
  }
  print("New haplotypes found and indexed")
  print(new_calls)
  updated_idx <- list()
  updated_idx$original <- arrange(rbind(idx$original, 
                                        new_calls), 
                                  locus, 
                                  index)
  
  updated_idx$clean <- arrange(rbind(idx$clean, 
                                     new_calls), 
                               locus, 
                               index)
  class(updated_idx) <- "index"
  return(updated_idx)
}


#' Recode haplotypes to numerically-coded genotypes.
#' 
#' Haplotypes in the `haplotypes` argument are converted to numerical alleles supplied by the `index` argument. Components of the returned `genotype` object other than `wide` and `long`  are inherited from the `haplotypes` argument.
#' @param hapltypes an object of class `haplotype`
#' @param index an object of class `index`
#' @return an object of class `genotype` containing the following components:
#' `loci` Names of the loci.
#' `allele_slots` The string used to identify allele slots within loci.
#' `allele_sep` The string used to separate locus names from allele slots
#' `wide` A table of genotypes in wide format.
#' `long` A table of genotypes in long format.
index_haplotypes <- function(haplotypes, index){
  stopifnot(class(index) == "index", class(haplotypes) == "haplotype")
  allele_sep <- haplotypes$allele_sep
  loci <- haplotypes$loci
  genos_long <- left_join(haplotypes$long, index$clean)
  #
  # Move the numerically-indexed allele into the `allele` column
  # Drop the haplotype version of the allele.
  #
  genos_long$allele <- genos_long$index
  genos_long <- genos_long[,grep("index", names(genos_long), invert = T)]
  #
  # Create wide version of genotypes
  # Ensure loci are in same order as in the haplotypes
  #
  pre_wide <- genos_long
  pre_wide$loc_allele <- paste(pre_wide$locus, 
                                 allele_sep,
                               pre_wide$allele_slot, 
                                 sep = "")
  pre_wide <- pre_wide[names(pre_wide)[!(names(pre_wide) %in% c("locus", "allele_slot"))]]
  wide <- pivot_wider(pre_wide, names_from = "loc_allele", values_from = "allele")
  ordered_col_names <- rep(loci, each = length(haplotypes$allele_slots))
  ordered_col_names <- paste(ordered_col_names, 
                             allele_sep, 
                             rep(haplotypes$allele_slots, 
                                 length.out = length(ordered_col_names)), 
                             sep = "")
  wide <- wide[c(names(wide)[(!(names(wide) %in% ordered_col_names))], ordered_col_names)]
  # return a genotypes object with recoded alleles
  out <- list()
  out$loci <- loci
  out$allele_slots <- haplotypes$allele_slots
  out$allele_sep <- allele_sep
  out$dropped_cols <- haplotypes$dropped_cols
  out$wide <- wide
  out$long <- genos_long
  class(out) <- "genotype"
  return(out)
  }

#' Write an index to a CSV file.
#' 
#' The `original` component of the index is written to file.
#' @param  index the index object to write from
#' @param filename the base name of the file, note the .csv extension will be added automatically
#' @param add_date logical indicating if the data should be added to the base file name
write_indices <- function(index, filename = "new_indices", add_date = T){
  fn <- filename
  if(add_date){
    fn <- c(fn, as.character(Sys.Date()))
    fn <- paste(fn, collapse = "_")
  }
  fn <- paste(c(fn, "csv"), collapse = ".")
  write_csv(index$original, file = fn)
}

#' Write genotypes to a CSV file
#' 
#' @param genos the genotypes object to write from
#' @param filename the name of the file to write to. Note that the `.csv` extension is added automatically
#' @param add_date whether to add the date to the output file
#' @param locus_no_first_slot_label whether the first allele slot of each locus is labeles (loc.1, loc.2) or not (loc, loc.1)
#' @param allele_sep the string to use to separate locus names from allele slots.
#' @param na.string string to use for missing values
#' @param loci the loci to include in the output. If `NULL` all loci in the genos argument will be used
#' @param restore_cols columns from the `genos$dropped_cols` to restore. Value should be a list with one element for each column to restore consisting of two numbers, the index of the column within `genos$dropped_cols` and the position of the column in the output.
#' @param rename_cols columns to rename,argument should be a list of pairs of strings `old_name`, `new_name` for each column to rename
write_genotypes_wide <- function(genos,
                                 filename = "fish_recoded_wide",
                                 add_date = T,
                                 locus_no_first_slot_label = T,
                                 allele_sep = ".",
                                 na.string = "NA",
                                 loci = NULL,
                                 restore_cols = NULL,
                                 rename_cols = NULL){
  stopifnot(class(genos) == "genotype")
  
  fn <- filename
  if(add_date){
    fn <- paste(fn, as.character(Sys.Date()), sep = "_")}
  fn <- paste(fn, "csv", sep = ".")
  
  if(!is.null(loci)){
    loci_chk <- loci #mabe unneceassary, but play it safe
    for(l in loci_chk){
      if(!(l %in% genos$loci)){
        message("Warning locus ", l, " was requested for output file ", fn, " but is missing from data, skipping.")
        lidx <- which(loci == l)
        loci <- loci[-lidx]
        }
      }
    }else{loci <- genos$loci}
  if(locus_no_first_slot_label){
    allele_slots <- c("", head(genos$allele_slots, length(genos$allele_slots) - 1))}else{
      allele_slots <- genos$allele_slots}
  # Note that the requested loci may not be in the same order as they appear in `genos$wide`
  #
  wide_cols <- 1
  for(l in loci){
    wide_cols <- c(wide_cols, grep(paste(l, genos$allele_sep, sep = ""), names(genos$wide), fixed = T)) 
    }
  wide <- genos$wide[wide_cols]
  new_header <- paste(rep(loci, each = length(allele_slots)), allele_slots, sep = allele_sep)
  #clean up trailing allele_sep if needed
  if(locus_no_first_slot_label){
    new_header <- gsub("\\.$", "", new_header)}
  # add in the individual ID column name
  new_header <- c(names(wide)[1], new_header)
  names(wide) <- new_header
  # restore dropped columns if needed
  if(!is.null(restore_cols)){
    to_restore_names <- names(genos$dropped_cols)[sapply(restore_cols, "[", 1)]
    to_restore <- genos$dropped_cols[to_restore_names]
    wide_restore <- cbind(to_restore, wide)
    #that's the easy part, now we need to put the restored columns in their proper locations
    restored_cols_ordered <- names(wide_restore)
    # Note that whatever the original order in genos$dropped_cols, the cols to restore
    # are ordered sequentially in the order they were requested in the restore_cols argument
    dest_positions <- sapply(restore_cols, "[", 2)
    for(p in dest_positions){
      #pop <- dest_positions[p]
      right <- restored_cols_ordered[(p + 1):length(restored_cols_ordered)]
      right <- c(restored_cols_ordered[1], right)
      left <- restored_cols_ordered[1:p]
      left <- left[-1]
      restored_cols_ordered <- c(left, right)
      wide <- wide_restore[restored_cols_ordered]
            }
  }
  # rename cols if needed
  if(!is.null(rename_cols)){
    for(i in 1:length(rename_cols)){
      n <- rename_cols[[i]]
      old_name <- n[1]
      new_name <- n[2]
      names(wide)[,names(wide == old_name)] <- new_name
      }
  }
  #change missing values to `na.string`
  wide <- as.data.frame(wide)
  wide[is.na(wide)] <- na.string
  write_csv(wide, fn)
}


###############################################################################
#                                                                             #
#                           Script body                                       #
#                                                                             #
###############################################################################


argo <- arg_parser(description = "Recode haplotypes to numerically coded alleles",
                   name = "cdfw_recoder.R",
                   hide.opts = T)
argo <- add_argument(argo, 
                     arg = "haplotypes", 
                     help = "A file of haplotypes, as output by the pipeline.")
argo <- add_argument(argo, 
                     arg = "index", 
                     help = "index of haplotype to allele number mappings.")
argo <- add_argument(argo, 
                     arg = "--RoSA", 
                     help = "A file containing a list of RoSA loci, used to write an output file with just the RoSA loci.")
argo <- add_argument(argo, 
                     arg = "--rubias", 
                     help = "File of loci to include for Rubias output.")
argo <- add_argument(argo, 
                     arg = "--colony", 
                     help = "File of loci to include for Colony output.")
argo <- add_argument(argo, "--no_date", 
                     flag = T, 
                     help = "Disable date in output file names")


args <- parse_args(argo)


##Debug

print(args)


# Read in the various inputs

haps <- read_haplotypes_wide(args$haplotypes, drop_cols = 1)
idx <- read_index(args$index)

if (!is.na(args$RoSA)){
  RoSA_loci <- read_lines(args$RoSA)
}

if(!is.na(args$rubias)){
  Rubias_loci <- read_lines(args$rubias)
}

if(!is.na(args$colony)){
  Colony_loci <- read_lines(args$colony)
  }
  

# Indexing

idx_new <- update_index(idx, haps)
genos <- index_haplotypes(haps, idx_new)
write_indices(idx)

# Outputs

# generic

write_genotypes_wide(genos,
                     add_date = !args$no_date,
                     restore_cols = list(c(1,1)),
                     na.string = "-9")

# RoSA only

if (!is.na(args$RoSA)){
  write_genotypes_wide(genos,
                       filename = "fish_recoded_wide_RoSA_only",
                       loci = RoSA_loci,
                       add_date = !args$no_date,
                       na.string = "-9")
}

# Rubias

if(!is.na(args$rubias)){
  write_genotypes_wide(genos,
                       filename = "fish_recoded_wide_Rubias",
                       loci = Rubias_loci,
                       add_date = !args$no_date,
                       na.string = "0",
                       restore_cols = list(c(1,1)))
}

if(!is.na(args$colony)){
  write_genotypes_wide(genos,
                       filename = "fish_recoded_wide_Colony",
                       loci = Colony_loci,
                       add_date = !args$no_date,
                       na.string = "0")
  }
