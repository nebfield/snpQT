#!/usr/bin/env Rscript

library('tidyverse')

# Args
# 1: plink.hwe file path before
# 2: plink.hwe file path after
# 3: threshold
# 4: plot title

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
    mutate(type = "before") -> before

read_table(args[[2]]) %>%
    mutate(type = "after") %>%
    bind_rows(before) %>%
    mutate(type = fct_relevel(type, "before")) -> hwe

fn <- ifelse(args[[4]] != "", "hwe_zoom.png", "hwe_sub.png")

ggplot(hwe, aes(x = P)) +
    geom_histogram() +
    geom_vline(xintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) +
    theme_linedraw() +
    xlab("P-value") +
    ylab("Variant count")+
    ggtitle(paste("Hardy-Weinberg Equilibrium (HWE) ", args[[4]]))
ggsave(fn) 
