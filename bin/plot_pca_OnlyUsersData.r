#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: case/control status

args <- commandArgs(trailingOnly = TRUE)

eigenvec<- read_delim(args[[1]], delim = " ")
status<- read_delim(args[[2]], delim = " ", col_names= FALSE)
colnames(status) <- c("FID", "IID", "status")

eigenvec %>%
    left_join(status, by = "IID") -> datafile

datafile %>%
  select(IID, PC1, PC2) %>%
  ggplot(., aes(x = PC1, y = PC2, colour = status)) + 
  scale_color_manual(breaks = c("1", "2"), 
      values=c("blue","red")) +
    scale_x_discrete(breaks=c("1", "2"), 
      labels=c("Control","Case")) +
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC2_onlyUsersData.png"))

datafile %>%
  select(IID, PC1, PC3) %>%
  ggplot(., aes(x = PC1, y = PC3, colour = status)) + 
  scale_color_manual(breaks = c("1", "2"), 
      values=c("blue","red")) +
  scale_x_discrete(breaks=c("1", "2"), 
      labels=c("Control","Case")) +
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC1vsPC3_onlyUsersData.png"))

datafile %>%
  select(IID, PC2, PC3) %>%
  ggplot(., aes(x = PC2, y = PC3, colour = status)) + 
  scale_color_manual(breaks = c("1", "2"), 
      values=c("blue","red")) +
  scale_x_discrete(breaks=c("1", "2"), 
      labels=c("Control","Case")) +
  geom_point() + 
  theme_linedraw()
ggsave(paste0("PC2vsPC3_onlyUsersData.png"))

pal <- c("red", "blue")
pal <- setNames(pal, c("Case", "Control")

fancy_plot <- plot_ly(
  datafile,
  x =  ~ PC1,
  y =  ~ PC2,
  z =  ~ PC3,
  size = 3,
  type = "scatter3d",
  mode = "markers",
  color =  ~ status,
  colors = pal
)
saveRDS(fancy_plot, paste0("plink_3D_PCA_onlyUsersData.rds"))
