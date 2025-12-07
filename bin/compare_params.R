#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)
library(purrr)
library(ggplot2)
library(glue)
library(rlist)
library(ggrepel)

################################################################################
# Inputs
################################################################################

args     <- commandArgs(trailingOnly = TRUE)
files    <- head(args, -2) # Results files, one per modkit parameter set
out_tsv  <- tail(args, 2)[1] # 2nd last arg is the filepath to write F1 scores
out_plot <- tail(args, 1) # Last arg is the plot filepath

################################################################################
# Load precision recall results for each modkit parameter set
################################################################################

# dir <- "/home/alex/Documents/repos/modkitopt/modkitopt/results"
# files <- list.files(path=dir, pattern='1{2}.*.tsv', full.names = TRUE)

ds_list <- list()
for (file in files)
{
  filter_threshold = str_split_i(file, "[_]", -2)
  mod_threshold = str_split_i(str_split_i(file, "[_]", -1), ".tsv", 1)

  # filter_threshold = str_split_i(str_split_i(file, "pileup_filter", 2), "[_]", 1)
  # mod_threshold = str_split_i(str_split_i(file, "_mod", 2), ".bed", 1)

  print(file)
  print(filter_threshold)
  print(mod_threshold)

  file %>%
    read_tsv() %>%
    mutate(params = glue("f: {filter_threshold}, m: {mod_threshold}"),
           f1 = 2 * precision * recall / (precision + recall)) ->
  ds

  ds_list <- list.append(ds_list, ds)
}

ds_to_plot <- map_dfr(ds_list, bind_rows)

################################################################################
# Plot the precision-recall curves for each modkit parameter set
################################################################################

best_f1 <- ds_to_plot %>% group_by(params) %>% slice_max(f1, n = 1) %>% ungroup()

# Write best F1 score per parameter set to .tsv file
write_tsv(best_f1, out_tsv)

theme_set(theme_classic())
ds_to_plot %>%
  ggplot(aes(x = recall, y = precision, colour = params)) +
    geom_line() +
    geom_point(data = best_f1, size = 3) +
    geom_label_repel(data = best_f1,
                     aes(label = sprintf("t=%.3f\nF1=%.3f", threshold, f1)),
                     family = "Helvetica",
                     label.size = NA,
                     show.legend = FALSE,
                     size = 2,
                     alpha = 0.8) +
    labs(x = "Recall",
         y = "Precision",
         fill = "Modkit parameters") +
    theme(text = element_text(family = "Helvetica"))

ggsave(out_plot, height = 5, width = 5)

################################################################################
# Return the modkit parameters and stoichiometry threshold that gives best F1
################################################################################

best_threshold <- best_f1 %>% slice_max(f1) %>% pull(threshold)
best_precision <- best_f1 %>% slice_max(f1) %>% pull(precision)
best_recall    <- best_f1 %>% slice_max(f1) %>% pull(recall)
best_f1_value  <- best_f1 %>% slice_max(f1) %>% pull(f1)
best_params    <- best_f1 %>% slice_max(f1) %>% pull(params)

# Extract best individual modkit parameters
best_filter_threshold <- str_split_i(str_split_i(best_params, ",", 1), " ", 2)
best_mod_threshold    <- str_split_i(str_split_i(best_params, ",", 2), " ", -1)

output <- "\n\nThe optimal modkit pileup parameters are:\n\n"
output <- paste(output, glue(">>> filter_threshold:\t{best_filter_threshold}\n\n"))
output <- paste(output, glue(">>> mod_threshold:\t{best_mod_threshold}\n\n\n"))
output <- paste(output, glue("With optimal stoichiometry threshold:\n\n\n"))
output <- paste(output, glue(">>> Threshold: {best_threshold}\n\n"))

cat(output)