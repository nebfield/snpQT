#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: racefile 
# 3: before or after removal

args <- commandArgs(trailingOnly = TRUE)

eigenvec<- read.table(args[[1]],header=TRUE)
race<- read.table(args[[2]],header=TRUE)
    
datafile<- merge(eigenvec,race,by=c("IID"))

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
