#!/usr/bin/env Rscript

library("tidyverse")

# Plot missingness per case/control status
# Args
# 1: plink.missing file path
# 2: missingness threshold
# 3: before or after threshold


args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> cs_missingness # plink is a dummy column

ggplot(cs_missingness, aes(x = P)) +
  geom_histogram() +
  geom_vline(xintercept = args[[2]], colour = "red")+
  theme_linedraw() +
  ylab("Variant count") +
  xlab("P-value") + 
  ggtitle("Missingness vs Case/Control status") -> cs_missingness_hist 


ggsave(paste0("cs_missingness_", args[[3]], ".png"), cs_missingness_hist, dpi = 300)