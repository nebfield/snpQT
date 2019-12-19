#!/usr/bin/env Rscript

library("tidyverse")

mds <- as_tibble(read.table("MDS_merge.mds", header = TRUE)) 
racefile <- as_tibble(read.table("racefile.txt", header = TRUE))

mds %>%
  left_join(racefile) %>%
  select(IID, race, C1, C2) %>%
  ggplot(., aes(x = C1, y = C2, colour = race)) + 
  geom_point() + 
  theme_classic()

ggsave("MDS.png")
