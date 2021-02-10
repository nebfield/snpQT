#!/usr/bin/env Rscript

library("tidyverse")

# Plot sample missingness 
# Args
# 1: plink.imiss file path
# 2: missingness threshold
# 3: before or after threshold

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(plink = "plink") -> sample_missingness # plink is a dummy column

n <- nrow(sample_missingness)

ggplot(sample_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[2]]), colour = "red")+
  theme_linedraw() + 
  xlab("Missing call rate") + 
  ylab("Sample count") +
  ggtitle("Sample missingness rate")
ggsave(paste0("sample_missingness_hist_", args[[3]], ".png"))

ggplot(sample_missingness, aes(x = plink, y = F_MISS)) +
  geom_jitter() +
  geom_vline(xintercept = as.numeric(args[[2]]), colour = "red")+
  theme_linedraw() + 
  ylab("Missing call rate") +
  xlab(glue::glue("Samples (n = {n})")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() +
  ggtitle("Sample missingness rate")

ggsave(paste0("sample_missingness_scatter_", args[[3]], ".png"))
