#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

# Check if correct number of arguments are passed
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript script_name.R <vsearch_input> <sklearn_input> <output_path> <primer>")
}


# Assign input and output file paths from arguments
vsearch_input <- args[1]
sklearn_input <- args[2]
output_path <- args[3]
primer<-args[4]

################################################################################
## step 1 - Import all VSEARCH classified clustered seqs and filter candidate features
################################################################################

vs_raw <- read_delim(file = vsearch_input,
                     delim = "\t", col_names = TRUE)

## remove 'Ambiguous taxa' from strings
vs_raw$Taxon <- gsub("Ambiguous_taxa", "", vs_raw$Taxon)

## remove sequences without any assignment
vs_filt <- vs_raw %>% 
  filter(Taxon != "Unassigned")
  ## 556 remaining seqs

## split 'Taxon' field into kingdom-->species levels
vs_filt <- vs_filt %>% 
  separate(Taxon, sep=';',
           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))


vs_filt <- vs_filt %>%
  mutate(across(where(is.character), ~ if_else(
    . == "" | str_detect(str_trim(.), "^\\w+__$"), 
    NA_character_, 
    .
  )))

## retain only sequences with at least Family-rank information
vs_filt <- vs_filt %>% 
  filter(!is.na(Order)) %>% 
  filter(!is.na(Family))

## add classifier label
vs_filt$Classifier <- "VSEARCH"

## create list of FeatureIDs to filter from qiime .qza taxonomy artifact
vslist <- as.data.frame(vs_filt$`Feature ID`)
colnames(vslist) <- 'Feature ID'

write.csv(vslist, file=file.path(output_path,paste0(primer,"_filtd_taxlist_vsearch.txt")),
          row.names = FALSE, quote = FALSE)

## cleanup
rm(vs_raw, vslist)

################################################################################
## step 2 - Import all sklearn classified clustered, ...
## ...then gather only those seqs not retained by VSEARCH above ...
## ...and filter using the same principles applied above
################################################################################

## import sklearn classified seq features
sk_raw <- read_delim(file = sklearn_input,
                     delim = "\t", col_names = TRUE)

## exclude vsearch features retained above from sklearn
sk_filt <- sk_raw %>% 
  filter(!`Feature ID` %in% vs_filt$`Feature ID`)

## filter using same approach as in VSEARCH above
sk_filt$Taxon <- gsub("Ambiguous_taxa", "", sk_filt$Taxon)

sk_filt <- sk_filt %>% 
  filter(Taxon != "Unassigned")

sk_filt <- sk_filt %>% 
  separate(Taxon, sep=';',
           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))


sk_filt <- sk_filt %>%
  mutate(across(where(is.character), ~ if_else(
    . == "" | str_detect(str_trim(.), "^\\w+__$"), 
    NA_character_, 
    .
  )))

sk_filt <- sk_filt %>% 
  filter(!is.na(Order)) %>% 
  filter(!is.na(Family))

## add classifier label
sk_filt$Classifier <- "sklearn"

## create list of FeatureIDs to filter from qiime .qza taxonomy artifact
sklist <- as.data.frame(sk_filt$`Feature ID`)

colnames(sklist) <- 'Feature ID'
write.csv(sklist, file=file.path(output_path,paste0(primer,"_filtd_taxlist_sklearn.txt")),
          row.names = FALSE, quote = FALSE)

## cleanup
rm(sk_raw)

################################################################################
## step 3 - combine taxonomy results to generate data table
## also, perform sanity check to ensure featureIDs are distinct in both sklearn and vsearch
################################################################################

## combine vsearch and sklearn into single data frame
all_dat <- bind_rows(sk_filt, vs_filt)

## sanity check we don't have duplicate Feature IDs
finalrows <- nrow(all_dat)
uniquerows <- length(unique(all_dat$`Feature ID`))
if (finalrows != uniquerows) {
  warning("Duplicate Feature IDs found.")
}

## write to disk
write.csv(all_dat, file = file.path(output_path,paste0(primer,"_filtd_tax_dataframe_ALL.csv")),
          row.names = FALSE, quote = FALSE)
