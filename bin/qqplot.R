#!/usr/bin/env Rscript

library("tidyverse")

# Plot sample missingness 
# Args
# 1: plink assoc.logistic file

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  as_tibble(.) %>%
  filter(TEST == "ADD") %>%
  arrange(P) %>%
  mutate(logobspval = -(log10(P))) %>%
  mutate(logexppval = -log10((row_number() - 0.5) / n())) %>%
  mutate(obsmax = trunc(max(logobspval))+1) %>%
  mutate(expmax = trunc(max(logexppval))+1) -> logistic

logistic %>%
  filter(L95 > 0) %>%
  mutate(beta = log(OR)) %>%
  mutate(z = (beta / SE) * (beta / SE)) %>%
  # 0.456: rounded median chi-squared distribution (qchisq(0.5, 1))
  # https://www.biostars.org/p/43328/
  summarise(median(z) / 0.456) %>% 
  pull(.) -> lambda.value

logistic %>%
  ggplot(aes(x = logexppval, y = logobspval)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("Expected -log10 p-value") +
  ylab("Observed -log10 p-value") +
  theme_linedraw() +
  # force origin to start at 0 
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  ggtitle(bquote(QQPlot~lambda==~.(lambda.value)))

ggsave("qqplot.png", dpi = 300)