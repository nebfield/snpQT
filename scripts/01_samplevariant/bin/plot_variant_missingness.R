#!/usr/bin/env Rscript

library("tidyverse")
library("gridExtra")

# Plot sample missingness 
# Args
# 1: plink.lmiss file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> variant_missingness # plink is a dummy column

ggplot(variant_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  theme_classic() +
  ggtitle("Proportion of sample missing per SNP") -> variant_hist 
  
n <- nrow(variant_missingness)
ggplot(variant_missingness, aes(x = plink, y = F_MISS)) +
  geom_jitter(alpha=0.2) +
  theme_classic() + 
  ylab("Proportion of sample missing per SNP") +
  xlab(glue::glue("Sample (n = {n})")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() -> variant_scatter

ggsave("variant_missingness.png", arrangeGrob(variant_hist, 
  variant_scatter, ncol=2), dpi = 300)
