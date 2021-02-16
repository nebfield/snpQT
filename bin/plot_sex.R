#!/usr/bin/env Rscript

library('tidyverse')
library('gridExtra')

# Args
# 1: before plink .sexcheck file path
# 2: after plink .sexcheck file path

args <- commandArgs(trailingOnly = TRUE)

read_table(args[[1]]) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(type = "before") -> before 

read_table(args[[2]]) %>%
  mutate(IID = as.factor(IID)) %>%
  mutate(type = "after") %>%
  bind_rows(before) %>%
  mutate(type = fct_relevel(type, "before")) -> gender 

gender %>%
  ggplot(.) +
	theme_linedraw() +
	facet_grid(~ type ) +
    geom_histogram(aes(x = F)) +
    ggtitle("All") -> sexcheck
ggsave("sexcheck_all_hist.png")

gender %>%
  filter(PEDSEX == 1) %>%
  ggplot(.) +
	theme_linedraw() + 
	facet_grid(~ type ) +
    geom_histogram(aes(x = F)) +
    ggtitle("Men") -> sexcheck_men
ggsave("sexcheck_men_hist.png")

gender %>%
  filter(PEDSEX == 2) %>%
  ggplot(.) +
	theme_linedraw() + 
	facet_grid(~ type ) +
    geom_histogram(aes(x = F)) +
    ggtitle("Women") -> sexcheck_women
ggsave("sexcheck_women_hist.png")

gender %>%
  mutate(PEDSEX = as.factor(PEDSEX)) %>%
  ggplot(.) +
	theme_linedraw() + 
	facet_grid(~ type ) +
    geom_point(aes(x=PEDSEX,y=F,colour=PEDSEX)) +
    scale_color_manual(breaks = c("0", "1", "2"), 
      values=c("black", "blue","red")) +
    scale_x_discrete(breaks=c("0", "1", "2"), 
      labels=c("Missing", "Males","Females")) +
    ggtitle("All") -> sexcheck_scatterplot
ggsave("sexcheck_all_scatter.png")

