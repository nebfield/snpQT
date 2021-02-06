#!/usr/bin/env Rscript

library("tidyverse")

# Plot sample missingness 
# Args
# 1: plink.imiss file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(plink = "plink") -> sample_missingness # plink is a dummy column

n <- nrow(sample_missingness)

ggplot(sample_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  theme_linedraw() + 
  xlab("Missing call rate") + 
  ylab("Sample count") +
  ggtitle("Sample missingness rate")
ggsave("sample_missingness_hist.png")

ggplot(sample_missingness, aes(x = plink, y = F_MISS)) +
  geom_jitter() +
  theme_linedraw() + 
  ylab("Missing call rate") +
  xlab(glue::glue("Samples (n = {n})")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() +
  ggtitle("Sample missingness rate")

ggsave("sample_missingness_scatter.png")