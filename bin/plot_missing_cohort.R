#!/usr/bin/env Rscript

library("tidyverse")

# Plot missingness per case/control status
# Args
# 1: plink.missing before file path
# 2: plink.missing after file path
# 3: missingness threshold

args <- commandArgs(trailingOnly = TRUE)


read_table(args[[1]]) %>%
    mutate(plink = "plink") %>% # dummy column
    mutate(type = "before") -> before

read_table(args[[2]]) %>%
    mutate(plink = "plink") %>%
    mutate(type = "after") %>%
    bind_rows(before) %>%
    mutate(type = fct_relevel(type, "before")) -> cs_missingness 

ggplot(cs_missingness, aes(x = P)) +
    geom_histogram() +
    geom_vline(xintercept = as.numeric(args[[3]]), colour = "red")+
    facet_grid(~ type) +
    theme_linedraw() +
    ylab("Variant count") +
    xlab("P-value") +
    ggtitle("Missingness vs Case/Control status")
ggsave("missingness_per_cohort.png")
