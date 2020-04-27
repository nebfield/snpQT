#!/usr/bin/env Rscript

library('tidyverse')
library('gridExtra')

# Args
# 1: plink.frq file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header=T) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> maf

ggplot(maf, aes(x = MAF)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Major allele frequency") + 
  ylab("Variant count") -> maf_hist

n_maf <- nrow(maf)
ggplot(maf, aes(x = plink, y = MAF)) +
  geom_jitter(alpha=0.1) +
  theme_classic() + 
  ylab("Proportion of sample missing per SNP") +
  xlab(glue::glue("Variant (n = {n_maf})")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() -> maf_scatter

ggsave("maf.png", arrangeGrob(maf_hist, 
  maf_scatter, ncol=2), dpi = 300)