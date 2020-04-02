#!/usr/bin/env Rscript

library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)

x <- as_tibble(read.table(args[[1]], header = TRUE))
y <- as_tibble(read.table(args[[2]], header = TRUE))

ggplot(x, aes(x = F_MISS)) +
  geom_histogram(binwidth=0.05) +
  theme_classic() + 
  xlab("Missing call rate") + 
  ylab("Sample count") +
  ggtitle("Missing call rate for all samples") + 
  expand_limits(x = c(0, 1)) 

ggsave("sample_callrate.png", device = "png", dpi = 300)

ggplot(y, aes(x = F_MISS)) +
  geom_histogram(binwidth=0.05) +
  theme_classic() + 
  xlab("Missing call rate") + 
  ylab("Variant count") +
  ggtitle("Missing call rate for all variants") + 
  expand_limits(x = c(0, 1)) +
  scale_y_continuous(label=scales::comma)

ggsave("variant_callrate.png", device = "png", dpi = 300)
