#!/usr/bin/env Rscript

library('tidyverse')

# Args
# 1: plink.frq file path
# 2: maf threshold
# 3: before or after threshold


args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header=T) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> maf

ggplot(maf, aes(x = MAF)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[2]]), colour = "red")+
  theme_linedraw() + 
  xlab("Minor allele frequency") + 
  ylab("Variant count") -> maf_hist

ggsave(paste0("maf", "_", args[[3]], ".png"), dpi = 300)
