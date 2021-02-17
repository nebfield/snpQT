#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec

args <- commandArgs(trailingOnly = TRUE)

eigenvec<- read_delim(args[[1]], delim = " ")

eigenvec %>%
  select(IID, PC1, PC2) %>%
  ggplot(., aes(x = PC1, y = PC2)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC2_onlyUsersData.png"))

eigenvec %>%
  select(IID, PC1, PC3) %>%
  ggplot(., aes(x = PC1, y = PC3)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC3_onlyUsersData.png"))

eigenvec %>%
  select(IID, PC2, PC3) %>%
  ggplot(., aes(x = PC2, y = PC3)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC2vsPC3_onlyUsersData.png"))

fancy_plot <- plot_ly(
  eigenvec,
  x =  ~ PC1,
  y =  ~ PC2,
  z =  ~ PC3,
  size = 3,
  type = "scatter3d",
  mode = "markers"
)
saveRDS(fancy_plot, paste0("plink_3D_PCA_onlyUsersData.rds"))
