#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: popfile 
# 3: before or after removal

args <- commandArgs(trailingOnly = TRUE)

eigenvec<- read_delim(args[[1]], delim = " ")
popfile<- read_delim(args[[2]], delim = " ")
# important to override colnames to make super / sub pop file input consistent
colnames(popfile) <- c("IID", "IID2", "pop")

eigenvec %>%
    left_join(popfile, by = "IID") %>%
    # missing pop means it's user data
    replace_na(list(pop = "OWN")) -> datafile

datafile %>%
  select(IID, pop, PC1, PC2) %>%
  ggplot(., aes(x = PC1, y = PC2, colour = pop)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC2_",args[[3]],".png"))

datafile %>%
  select(IID, pop, PC1, PC3) %>%
  ggplot(., aes(x = PC1, y = PC3, colour = pop)) + 
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC3_",args[[3]],".png"))

datafile %>%
  select(IID, pop, PC2, PC3) %>%
  ggplot(., aes(x = PC2, y = PC3, colour = pop)) + 
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
  color =  ~ pop
)
saveRDS(fancy_plot, paste0("plink_3D_pca", args[[3]], ".rds"))
