#!/usr/bin/env Rscript

library("tidyverse")

# Plot sample missingness 
# Args
# 1: plink.lmiss file path
# 2: missingness threshold
# 3: before or after threshold

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  mutate(plink = "plink") -> variant_missingness # plink is a dummy column

n <- nrow(variant_missingness)

ggplot(variant_missingness, aes(x = F_MISS)) +
  geom_histogram() +
  geom_vline(xintercept = as.numeric(args[[2]]), colour = "red")+
  theme_linedraw() +
  ylab("Variant count") +
  xlab("Missing call rate") + 
  ggtitle("Variant missingness rate")  
ggsave(paste0("variant_missingness_hist_", args[[3]], ".png")) 

ggplot(variant_missingness, aes(x = plink, y = F_MISS)) +
  geom_jitter(alpha=0.2) +
  geom_vline(xintercept = as.numeric(args[[2]]), colour = "red")+
  theme_linedraw() + 
  ggtitle("Variant missingness rate") + 
  ylab("Missing call rate") +
  xlab(glue::glue("Variant (n = {n})")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() 

ggsave(paste0("variant_missingness_scatter_", args[[3]], ".png"))
