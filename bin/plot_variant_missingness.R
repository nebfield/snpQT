#!/usr/bin/env Rscript

library("tidyverse")
library(cowplot)

# Plot sample missingness 
# Args
# 1: plink.lmiss before file path
# 2: plink.lmiss after file path
# 3: missingness threshold

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
    mutate(type = "before") -> before

read_table(args[[2]]) %>%
    mutate(type = "after") %>%
    bind_rows(before) %>%
    mutate(type = fct_relevel(type, "before")) -> variant_missingness

variant_missingness %>%
  filter(type == "before") %>%
  count() %>%
  pull(n) -> n

ggplot(variant_missingness, aes(x = F_MISS)) +
    geom_histogram() +
    geom_vline(xintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) +
    theme_cowplot() +
	background_grid()+
	panel_border() +
	scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    ylab("Variant count") +
    xlab("Missing call rate") +
    ggtitle("Variant missingness rate")  
ggsave("variant_missingness_hist.png", height = 7, width = 12)

variant_missingness %>%
    ggplot(aes(x = SNP, y = F_MISS)) +
    geom_jitter(alpha=0.3) +
    geom_hline(yintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) + 
    theme_cowplot() +
	panel_border() +
    ggtitle("Variant missingness rate") +
    ylab("Missing call rate") +
    xlab(glue::glue("Variant (n = {n})")) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    coord_flip() 
ggsave("variant_missingness_scatter.png", height = 7, width = 12)
