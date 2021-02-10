#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: racefile 
# 3: before or after removal

args <- commandArgs(trailingOnly = TRUE)

eigenvec<- read_delim(args[[1]], delim = " ")
racefile<- read_delim(args[[2]], delim = " ")
# important to override colnames to make super / sub race file input consistent
colnames(racefile) <- c("IID", "IID2", "race")

eigenvec %>%
    left_join(racefile, by = "IID") %>%
    # missing race means it's user data
    replace_na(list(race = "OWN")) -> datafile

datafile %>%
  select(IID, race, PC1, PC2) %>%
  ggplot(., aes(x = PC1, y = PC2, colour = race)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC2_",args[[3]],".png"))

datafile %>%
  select(IID, race, PC1, PC3) %>%
  ggplot(., aes(x = PC1, y = PC3, colour = race)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC3_",args[[3]],".png"))

datafile %>%
  select(IID, race, PC2, PC3) %>%
  ggplot(., aes(x = PC2, y = PC3, colour = race)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC2vsPC3_",args[[3]],".png"))

fancy_plot <- plot_ly(
  datafile,
  x =  ~ PC1,
  y =  ~ PC2,
  z =  ~ PC3,
  size = 3,
  type = "scatter3d",
  mode = "markers",
  color =  ~ race
)

saveWidget(fancy_plot, paste0("C6_3d_",args[[3]],".html"), selfcontained = T)
