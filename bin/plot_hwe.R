#!/usr/bin/env Rscript

library('tidyverse')

# Args
# 1: plink.hwe file path
# 2: plot title

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = T) %>% 
  as_tibble(.) -> hwe

hwe %>%
  ggplot(., aes(x = P)) +
  geom_histogram() +
  theme_linedraw() +
  xlab("P-value") + 
  ylab("Variant count")+
  ggtitle(paste("Hardy-Weinberg Equilibrium (HWE) ", args[[2]]))
ggsave(paste0(args[[1]], ".png"))