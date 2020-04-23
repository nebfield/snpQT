#!/usr/bin/env Rscript

library("tidyverse")
library("gridExtra")

# Plot sample missingness 
# Args
# 1: plink.imiss file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(plink = "plink") -> sample_missingness # plink is a dummy column

ggplot(sample_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Missing call rate") + 
  ylab("Sample count") +
  ggtitle("Individual missingness") -> individual_missingness

n <- nrow(sample_missingness)
ggplot(sample_missingness, aes(x = plink, y = F_MISS)) +
  geom_jitter() +
  theme_classic() + 
  ylab("Frequency of missing variants") +
  xlab(glue::glue("Sample (n = {n})")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() +
  ggtitle("Individual missingness") -> missing_variants

ggsave("sample_missingness.png", arrangeGrob(individual_missingness, 
  missing_variants, ncol=2), dpi = 300)
