#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: log files parsed by parse_logs.awk
# 2: output prefix name (e.g. sample_qc)

args <- commandArgs(trailingOnly = TRUE)

read_delim(args[[1]], delim = " ") %>%
    mutate(variants = as.integer(variants)) %>%
    mutate(samples = as.integer(samples)) %>%
    mutate(num = as.integer(str_extract(stage, "\\d+"))) %>%
    mutate(stage = word(stage, sep = "\\.")) %>%
    mutate(stage = fct_reorder(stage, num)) -> df

df %>% 
  ggplot(., aes(x = stage, y = variants, group = 1)) + 
  geom_point() + 
  scale_y_continuous(breaks= pretty_breaks())+
  geom_line() + 
  theme_linedraw() +
  labs(x = "Stage", y = "Number of variants")
ggsave(paste0(args[[2]], "_variants.png"))

df %>% 
  ggplot(., aes(x = stage, y = samples, group = 1)) + 
  geom_point() + 
  scale_y_continuous(breaks= pretty_breaks())+
  geom_line() + 
  theme_linedraw() +
  labs(x = "Stage", y = "Number of samples")
ggsave(paste0(args[[2]], "_samples.png"))
