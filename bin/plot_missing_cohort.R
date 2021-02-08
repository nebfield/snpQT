#!/usr/bin/env Rscript

library("tidyverse")

# Plot missingness per case/control status
# Args
# 1: plink.missing file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> cs_missingness # plink is a dummy column

ggplot(cs_missingness, aes(x = P)) +
  geom_histogram() +
  theme_linedraw() +
  ylab("Variant count") +
  xlab("P-value") + 
  ggtitle("Missingness vs Case/Control status") -> cs_missingness_hist 


ggsave("cs_missingness.png", cs_missingness_hist, dpi = 300)