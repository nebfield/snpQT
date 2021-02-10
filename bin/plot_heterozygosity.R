#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: before .het file path
# 2: after .het file path 

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
    mutate(HET_RATE = (`N(NM)` - `O(HOM)`) / `N(NM)`) %>%
    mutate(type = "before") -> before

read_table(args[[2]]) %>%
    mutate(HET_RATE = (`N(NM)` - `O(HOM)`) / `N(NM)`) %>%
    mutate(type = "after") %>%
    bind_rows(before) %>%
    mutate(type = fct_relevel(type, "before")) -> het

het %>%
  filter(type == "before") %>%
  summarise(thresh_max = mean(HET_RATE) + (3 * sd(HET_RATE))) %>%
  pull(thresh_max) -> thresh_max

het %>%
  filter(type == "before") %>%
  summarise(thresh_min = mean(HET_RATE) - (3 * sd(HET_RATE))) %>%
  pull(thresh_min) -> thresh_min

ggplot(het, aes(x = IID, y = HET_RATE)) + 
  geom_point() +
  geom_hline(yintercept = thresh_max, colour = "red") +
  geom_hline(yintercept = thresh_min, colour = "red") +
  facet_grid(~ type) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("Sample") +
  ylab("Heterozygosity rate") +
  ggtitle("Heterozygosity rate per sample")
ggsave("heterozygosity_rate.png")

