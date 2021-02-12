#!/usr/bin/env Rscript

library("tidyverse")

# Plot sample missingness 
# Args
# 1: before plink imiss file path
# 2: after plink imiss file path
# 3: missingness threshold

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(type = "before") -> before 

read_table(args[[2]]) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(type = "after") %>%
  bind_rows(before) %>%
  mutate(type = fct_relevel(type, "before")) -> sample_missingness 

n <- nrow(sample_missingness)

ggplot(sample_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[3]]), colour = "red")+
  theme_linedraw() +
  facet_grid(~ type ) +
  xlab("Missing call rate") + 
  ylab("Sample count") +
  ggtitle("Sample missingness rate")
ggsave("sample_missingness_hist.png")

sample_missingness %>%
  ggplot(aes(x = IID, y = F_MISS)) +
    geom_jitter(alpha=0.2) +
    geom_hline(yintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) + 
    theme_bw() +  
    ylab("Missing call rate") +
    xlab(glue::glue("Samples (n = {n})")) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    coord_flip() +
    ggtitle("Sample missingness rate")

ggsave("sample_missingness_scatter.png")
