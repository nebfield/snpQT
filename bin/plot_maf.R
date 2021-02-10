#!/usr/bin/env Rscript

library('tidyverse')

# Args
# 1: plink.frq before file path
# 2: plink.frq after file path 
# 3: maf threshold

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
    mutate(type = "before") -> before

read_table(args[[2]]) %>%
    mutate(type = "after") %>%
    bind_rows(before) %>%
    mutate(type = fct_relevel(type, "before")) %>%
    mutate(plink = "plink") -> maf

ggplot(maf, aes(x = MAF)) +
    geom_histogram() +
    geom_vline(xintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) +
    theme_linedraw() +
    xlab("Minor allele frequency") +
    ylab("Variant count") 
ggsave("maf.png")
