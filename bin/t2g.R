#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(glue)
library(GenomicFeatures)

# Inputs
args = commandArgs(trailingOnly = TRUE)
file_bed   <- args[1]
annotation <- args[2]

# Output file
outfile <- glue("{file_bed}.genomic")
if (file.exists(outfile)) file.remove(outfile)

# Build the reference database
txdb  <- makeTxDbFromGFF(annotation)

# Get all exons by transcript
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)

# Chunk size for processing bed file memory-efficiently
chunk_size <- 1e6

# Chunk processing
process_chunk <- function(df, pos) {
  message(glue("Processing chunk starting line {pos} with {nrow(df)} rows"))

  df <- df %>%
    # bedMethyl file header starts with #, so remove from column names
    rename_with(~ sub("^#", "", .)) %>%
    # Clean up chrom (transcript ID)
    mutate(chrom = str_split_i(chrom, "[|]", 1))

  # Get the transcriptomic sites
  tx_coords <- GRanges(
    seqnames = df$chrom,
    # bedMethyl is 0-based, gtf/gff3 is 1-based, so mod position is chromEnd in bedMethyl
    ranges   = IRanges(start = df$chromEnd, end = df$chromEnd),
    valid_coverage   = df$valid_coverage,
    percent_modified = df$percent_modified
  )

  # Map transcriptomic sites to genomic sites
  genomic_coords <- mapFromTranscripts(tx_coords, exons)

  # Record unmapped sites for this chunk
  mapped_indices   <- unique(mcols(genomic_coords)$xHits)
  all_indices      <- seq_along(tx_coords)
  unmapped_indices <- setdiff(all_indices, mapped_indices)
  message(glue("Chunk had {length(unmapped_indices)} unmapped coordinates"))

  # Prepare metadata to add to genomic sites
  df_tx_coords  <- as_tibble(as.data.frame(tx_coords)) %>% rowid_to_column("gr_number")
  df_gen_coords <- as_tibble(as.data.frame(genomic_coords))

  # Add metadata to genomic sites
  df_joined <- df_gen_coords %>%
    left_join(df_tx_coords, by = join_by("xHits" == "gr_number"), suffix = c("_gen", "_txt"))

  # Append to output
  write_tsv(df_joined, outfile, append = TRUE)

  # Return null since callback has to return something
  invisible()
}

# Iterate over chunks
message("Starting chunked processing of BED data…")
read_tsv_chunked(
  file_bed,
  callback   = SideEffectChunkCallback$new(process_chunk),
  chunk_size = chunk_size
)

message("Finished. Output written to: ", outfile)
