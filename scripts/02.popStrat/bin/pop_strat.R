#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: eigenvec
# 2: racefile 

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header=TRUE) %>%
  as_tibble(.) -> eigenvec

read.table(args[[2]], header=TRUE) %>%
  as_tibble(.) %>%
  left_join(eigenvec) -> dat

dat %>%
  select(IID, race, PC1, PC2) %>%
  ggplot(., aes(x = PC1, y = PC2, colour = race)) + 
  geom_point() + 
  theme_classic()
ggsave("PCA.png", dpi = 300)