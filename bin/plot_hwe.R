#!/usr/bin/env Rscript

library('tidyverse')

# Args
# 1: plink.hwe file path
# 2: plot title
# 3: hwe threshold
# 4: before or after threshold

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = T) %>% 
  as_tibble(.) -> hwe

hwe %>%
  ggplot(., aes(x = P)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[3]]), colour = "red")+
  theme_linedraw() +
  xlab("P-value") + 
  ylab("Variant count")+
  ggtitle(paste("Hardy-Weinberg Equilibrium (HWE) ", args[[2]]))
ggsave(paste0(args[[1]], "_", args[[4]], ".png"))
