#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(magrittr)
library(glue)
library(GenomicFeatures)

# Inputs
args = commandArgs(trailingOnly=TRUE)
file_bed <- args[1]
gff      <- args[2]

# Load modkit output
file_bed %>%
  read_tsv() %>%
  # bedMethyl file header starts with #, so remove from column names
  rename_with(~ sub("^#", "", .)) %>%
  # Clean up chrom (transcript ID)
  mutate(chrom = str_split_i(chrom, "[|]", 1)) ->
ds

# Build the reference database
txdb  <- makeTxDbFromGFF(gff, format = "gff3")

# Get all exons by transcript
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)

# Get the transcriptomic sites
tx_coords <- GRanges(seqnames = ds$chrom,
                     ranges = IRanges(start = ds$chromStart, end = ds$chromEnd),
                     valid_coverage = ds$valid_coverage,
                     percent_modified = ds$percent_modified)

# Map transcriptomic sites to genomic sites
genomic_coords <- mapFromTranscripts(tx_coords, exons)

# Determine if any transcriptomic sites were unmapped
mapped_indices   <- unique(mcols(genomic_coords)$xHits)
all_indices      <- seq_along(tx_coords)
unmapped_indices <- setdiff(all_indices, mapped_indices)
print(glue("{file_bed} had {length(unmapped_indices)} transcript coordinates unmapped"))

# Add metadata to genomic sites
ds_tx_coords <- as_tibble(as.data.frame(tx_coords)) %>% rowid_to_column("gr_number")
ds_gen_coords <- as_tibble(as.data.frame(genomic_coords))

ds_joined <- left_join(ds_gen_coords,
                       ds_tx_coords,
                       by = join_by("xHits" == "gr_number"),
                       suffix = c("_gen", "_txt"))
write_tsv(ds_joined, glue("{file_bed}.genomic"))
