#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: before .het file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
    mutate(HET_RATE = (N.NM. - O.HOM.) / N.NM.) -> het

het %>%
  filter(HET_RATE < mean(HET_RATE) - 3 * sd(HET_RATE) |
         HET_RATE > mean(HET_RATE) + 3 * sd(HET_RATE)) %>%
  mutate(HET_DST = HET_RATE - mean(HET_RATE / sd(HET_RATE))) %>%
  write_tsv("het_failed_samples.txt") 
