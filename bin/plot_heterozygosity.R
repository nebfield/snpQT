#!/usr/bin/env Rscript

library("tidyverse")
library(cowplot)

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
  group_by(type) %>%
  summarise(thresh_max = mean(HET_RATE) + (3 * sd(HET_RATE))) -> thresh_max

het %>%
  group_by(type) %>%
  summarise(thresh_min = mean(HET_RATE) - (3 * sd(HET_RATE))) %>%
  left_join(thresh_max) -> thresholds

ggplot(het, aes(x = IID, y = HET_RATE)) + 
  geom_jitter(alpha=0.3) +
  geom_hline(data = thresholds, aes(yintercept = thresh_max), colour = "red") +
  geom_hline(data = thresholds, aes(yintercept = thresh_min), colour = "red") +
  facet_grid(~ type) +
  theme_cowplot() +
  panel_border() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("Samples") +
  ylab("Heterozygosity rate") +
  ggtitle("Heterozygosity rate per sample")
ggsave("heterozygosity_rate_scatter.png")

ggplot(het, aes(x = HET_RATE)) +
    geom_histogram() +
    geom_vline(data = thresholds, aes(xintercept = thresh_max), colour = "red") +
    geom_vline(data = thresholds, aes(xintercept = thresh_min), colour = "red") +
	facet_grid(~ type) +
    theme_cowplot() +
	background_grid()+
	panel_border() +
	scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    xlab("Heterozygosity rate") +
    ylab("Samples")+
    ggtitle("Heterozygosity rate per sample")
ggsave("heterozygosity_rate_hist.png")

