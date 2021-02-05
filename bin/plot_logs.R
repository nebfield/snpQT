#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: log files parsed by parse_logs.awk
# 2: output prefix name (e.g. sample_qc)

args <- commandArgs(trailingOnly = TRUE)

df <- read_delim(args[[1]], col_names = FALSE, delim = " ")
colnames(df) <-
  c("stage",
    "variants",
    "samples",
    "pheno",
    "pheno_case",
    "pheno_control",
    "pheno_miss",
    "wd")
df %>% 
  mutate(stage = word(stage, sep = "\\.")) %>%
  ggplot(., aes(x = stage, y = variants, group = 1)) + 
  geom_point() + 
  geom_line() + 
  theme_bw() +
  labs(x = "Stage", y = "Number of variants")
ggsave(paste0(args[[2]], "_variants.png"))

df %>% 
  mutate(stage = word(stage, sep = "\\.")) %>%
  ggplot(., aes(x = stage, y = samples, group = 1)) + 
  geom_point() + 
  geom_line() + 
  theme_bw() +
  labs(x = "Stage", y = "Number of samples")
ggsave(paste0(args[[2]], "_samples.png"))
