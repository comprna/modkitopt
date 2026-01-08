#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)
library(purrr)
library(ggplot2)
library(stringr)
library(glue)
library(rlist)
library(ggrepel)
library(paletteer)
library(tidyr)

################################################################################
# Inputs
################################################################################

args            <- commandArgs(trailingOnly = TRUE)
files           <- head(args, -4) # Results files, one per modkit parameter set
out_tsv         <- tail(args, 4)[1] # 4th last arg is the filepath to write F1 scores
out_pr_plot     <- tail(args, 3)[1] # 3rd last arg is the precision-recall curve filepath
out_barplot     <- tail(args, 2)[1] # 2nd last arg is the barplot filepath
out_scatterplot <- tail(args, 1) # Last arg is the scatterplot filepath

################################################################################
# Load precision recall results for each modkit parameter set
################################################################################

ds_list <- list()
for (file in files)
{
  filter_threshold = str_split_i(file, "[_]", -2)
  mod_threshold = str_split_i(str_split_i(file, "[_]", -1), ".tsv", 1)

  file %>%
    read_tsv() %>%
    mutate(params = glue("{filter_threshold}, {mod_threshold}"),
           f1 = 2 * precision * recall / (precision + recall)) ->
  ds

  ds_list <- list.append(ds_list, ds)
}

ds_to_plot <- map_dfr(ds_list, bind_rows)

################################################################################
# Plot the precision-recall curves for each modkit parameter set
################################################################################

ds_to_plot %>%
  filter(params != "default, default") %>%
  group_by(params) %>%
  slice_max(f1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  relocate(params, threshold, precision, recall, f1) %>%
  arrange(desc(f1)) %T>%
  write_tsv(out_tsv) ->
best_f1

# We only want to label the top F1 score on the plot to avoid clutter
f1_to_plot <- best_f1 %>% filter(f1 == max(f1))

theme_set(theme_classic())
n_colours <- length(unique(ds_to_plot$params))
ds_to_plot %>%
  filter(params != "default, default") %>%
  ggplot(aes(x = recall, y = precision, colour = params)) +
    geom_line() +
    geom_point(data = best_f1, size = 3) +
    scale_color_paletteer_d("Polychrome::palette36") +
    geom_label_repel(data = f1_to_plot,
                     aes(label = sprintf("%s\nthreshold = %.3f\nF1 = %.3f",
                                         params,
                                         threshold,
                                         f1)),
                     family = "Helvetica",
                     label.size = NA,
                     show.legend = FALSE,
                     size = 4,
                     alpha = 0.8) +
    labs(x = "Recall",
         y = "Precision",
         colour = "Modkit parameters") +
    theme(text = element_text(family = "Helvetica"),
          legend.position = "bottom")

ggsave(out_pr_plot, height = 10, width = 8)

################################################################################
# Plot precision, recall, F1 for each parameter set (and corresponding best
# stoichiometry cut-off)
################################################################################

brewer_tan <- "#E5C494"
br_pink <- "#ff00a8"
br_green <- "#00af8d"
mk_violet <- "#8400CD"

best_f1 %>%
  select(params, f1, precision, recall, threshold) %>%
  arrange(desc(f1)) %>%
  mutate(params = factor(params, levels = params)) %>%
  pivot_longer(cols = c("recall", "precision", "f1", "threshold"), values_to = "metric") %>%
  mutate(name = factor(name, levels = c("f1", "precision", "recall", "threshold"))) %>%
  ggplot(aes(x = params, y = metric, fill = name)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    labs(x = "Modkit parameters (--filter_threshold, --mod_threshold)",
         y = "Value",
         fill = "Metric") +
    scale_fill_manual(labels = c("F1", "Precision", "Recall", "Stoich. cutoff"),
                      values = c(mk_violet, br_pink, br_green, brewer_tan)) +
    ylim(0, 1) +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
          legend.position = "bottom",
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24))

ggsave(out_barplot, height = 8, width = 14)

################################################################################
# Plot the F1 score across the threshold parameter space
################################################################################

mk_purple_tyrian <- "#68023F"
br_yellow <- "#fdb027"

best_f1 %>%
  filter(params != "default, default") %>%
  separate_wider_delim(params, delim = ", ", names = c("filter_threshold", "mod_threshold")) %>%
  mutate(filter_threshold = as.numeric(filter_threshold),
         mod_threshold    = as.numeric(mod_threshold)) ->
to_plot_f1

to_plot_f1_best <- to_plot_f1 %>% slice_max(f1)

to_plot_f1 %>%
  ggplot(aes(x = filter_threshold, y = mod_threshold, size = f1, fill = f1)) +
    geom_point(shape = 21, colour = "black", stroke = 0.3) +
    # Add a star on top of the best performing parameter combo
    geom_point(data = to_plot_f1_best,
               aes(x = filter_threshold, y = mod_threshold),
               shape = 18,
               colour = br_yellow,
               fill = br_yellow,
               size = 6,
               stroke = 1,
               show.legend = FALSE) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
    scale_fill_gradient(low = "white", high = mk_purple_tyrian) +
    scale_x_continuous(limits = c(0.4, 0.9959),
                       breaks = c(0.5, 0.75, 0.9, 0.95, 0.99, 0.995),
                       transform = "logit") +
    scale_y_continuous(limits = c(0.4, 0.9959),
                       breaks = c(0.5, 0.75, 0.9, 0.95, 0.99, 0.995),
                       transform = "logit") +
    guides(size = "legend", fill = "legend") +
    theme(text = element_text(family = "Helvetica"),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out_scatterplot, height = 4, width = 4)

################################################################################
# Return the modkit parameters and stoichiometry threshold that gives best F1
################################################################################

best_threshold <- best_f1 %>% slice_max(f1) %>% pull(threshold)
best_precision <- best_f1 %>% slice_max(f1) %>% pull(precision)
best_recall    <- best_f1 %>% slice_max(f1) %>% pull(recall)
best_f1_value  <- best_f1 %>% slice_max(f1) %>% pull(f1)
best_params    <- best_f1 %>% slice_max(f1) %>% pull(params)

# Extract best individual modkit parameters
best_filter_threshold <- str_split_i(best_params, ",", 1)
best_mod_threshold    <- str_split_i(best_params, " ", 2)

if (length(best_params) > 1)
{
  output <- "\n\nThe optimal modkit pileup parameters are as follows"
  output <- paste(output, "(in this case, more than one parameter set achieved the highest F1 score):\n\n")
  output <- paste(output, glue(">>> filter_threshold:\t{best_filter_threshold}\n\n"))
  output <- paste(output, glue(">>> mod_threshold:\t{best_mod_threshold}\n\n\n"))
  output <- paste(output, glue("With the optimal stoichiometry threshold to filter sites:\n\n\n"))
  output <- paste(output, glue(">>> Threshold: {best_threshold}\n\n\n"))
  output <- paste(output, glue("Achieving an F1 score of {round(best_f1_value, 3)}\n\n"))
} else
{
  output <- "\n\nThe optimal modkit pileup parameters are:\n\n"
  output <- paste(output, glue(">>> filter_threshold:\t{best_filter_threshold}\n\n"))
  output <- paste(output, glue(">>> mod_threshold:\t{best_mod_threshold}\n\n\n"))
  output <- paste(output, glue("With the optimal stoichiometry threshold to filter sites:\n\n\n"))
  output <- paste(output, glue(">>> Threshold: {best_threshold}\n\n\n"))
  output <- paste(output, glue("Achieving an F1 score of {round(best_f1_value, 3)}\n\n"))
}

cat(output)
