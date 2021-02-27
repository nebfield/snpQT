#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: log files parsed by parse_logs.awk
# 2: output prefix name (e.g. sample_qc)

args <- commandArgs(trailingOnly = TRUE)

# A function factory for getting integer y-axis values.
integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
        breaks <- floor(pretty(x, n, ...))
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}

read_delim(args[[1]], delim = " ") %>%
    mutate(variants = as.integer(variants)) %>%
    mutate(samples = as.integer(samples)) %>%
    mutate(num = as.integer(str_extract(stage, "\\d+"))) %>%
    mutate(stage = word(stage, sep = "\\.")) %>%
    mutate(stage = fct_reorder(stage, num)) -> df

df %>% 
  ggplot(., aes(x = stage, y = as.integer(variants), group = 1)) + 
  geom_point() + 
  scale_y_continuous(breaks= integer_breaks())+
  geom_line() + 
  theme_linedraw() +
  labs(x = "Stage", y = "Number of variants")
ggsave(paste0(args[[2]], "_variants.png"))

df %>% 
  ggplot(., aes(x = stage, y = as.integer(samples), group = 1)) + 
  geom_point() + 
  scale_y_continuous(breaks= integer_breaks())+
  geom_line() + 
  theme_linedraw() +
  labs(x = "Stage", y = "Number of samples")
ggsave(paste0(args[[2]], "_samples.png"))
