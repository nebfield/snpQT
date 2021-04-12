#!/usr/bin/env Rscript

library("tidyverse")

# Plot sample missingness 
# Args
# 1: plink gwas file
# 2: lambda value

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  filter(TEST == "ADD") %>%
  arrange(P) %>%
  mutate(logobspval = -(log10(P))) %>%
  mutate(logexppval = -log10((row_number() - 0.5) / n())) %>%
  mutate(obsmax = trunc(max(logobspval))+1) %>%
  mutate(expmax = trunc(max(logexppval))+1) -> gwas

lambda <- read.table(args[[2]], header = FALSE)

	

gwas %>%
  ggplot(aes(x = logexppval, y = logobspval)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("Expected -log10 p-value") +
  ylab("Observed -log10 p-value") +
  theme_linedraw() +
  # force origin to start at 0 
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  ggtitle(bquote(QQPlot~lambda==~.(last(lambda))))

ggsave("qqplot.png", dpi = 300)