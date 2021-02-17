#!/usr/bin/env Rscript

library("tidyverse")
library(cowplot)

# Plot sample missingness 
# Args
# 1: before plink imiss file path
# 2: after plink imiss file path
# 3: missingness threshold

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(plink = "plink") %>% # dummy column
  mutate(type = "before") -> before 

read_table(args[[2]]) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(plink = "plink") %>%
  mutate(type = "after") %>%
  bind_rows(before) %>%
  mutate(type = fct_relevel(type, "before")) -> sample_missingness 

sample_missingness %>%
  filter(type == "before") %>%
  count() %>%
  pull(n) -> n
  
ggplot(sample_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[3]]), colour = "red")+
  theme_cowplot() +
  background_grid()+
  panel_border() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  facet_grid(~ type ) +
  xlab("Missing call rate") + 
  ylab("Sample count") +
  ggtitle("Sample missingness rate")
ggsave("sample_missingness_hist.png", height = 7, width = 10)

sample_missingness %>%
  ggplot(aes(x = plink, y = F_MISS)) +
    geom_jitter(alpha=0.3) +
    geom_hline(yintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) + 
    theme_cowplot() +
	panel_border() +
    ylab("Missing call rate") +
    xlab(glue::glue("Samples (n = {n})")) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    coord_flip() +
    ggtitle("Sample missingness rate")

ggsave("sample_missingness_scatter.png", height = 7, width = 10)
