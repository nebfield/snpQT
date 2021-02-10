#!/usr/bin/env Rscript

library('tidyverse')

# Args
# 1: plink.hwe file path
# 2: plot title
# 3: hwe threshold
# 4: before or after threshold
# 5: sub or zoom

args <- commandArgs(trailingOnly = TRUE)

hwe <- read_table(args[[1]])

fn <- ifelse(args[[5]] == "zoom", "zoom", "sub")

ggplot(hwe, aes(x = P)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[3]]), colour = "red") +
  theme_linedraw() +
  xlab("P-value") + 
  ylab("Variant count")+
  ggtitle(paste("Hardy-Weinberg Equilibrium (HWE) ", args[[2]]))
ggsave(paste0("hwe_", fn, "_", args[[4]], ".png"), device = "png", height = 7, width = 7)
