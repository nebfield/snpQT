#!/usr/bin/env Rscript

library('tidyverse')
library('gridExtra')

# Args
# 1: plink.frq file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header=T) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> maf

ggplot(maf, aes(x = MAF)) +
  geom_histogram() +
  theme_linedraw() + 
  xlab("Minor allele frequency") + 
  ylab("Variant count") -> maf_hist

ggsave("maf.png", dpi = 300)
