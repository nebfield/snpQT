#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: case/control status

args <- commandArgs(trailingOnly = TRUE)

eigenvec<- read_delim(args[[1]], delim = " ")
statusfile<- read_delim(args[[2]], delim = " ", col_names= FALSE)
colnames(statusfile) <- c("FID", "IID", "status")

statusfile %>%
  mutate(status = fct_recode(
    as.character(status),
    "control" = "1",
    "case" = "2"
  )) -> statusfile
  
eigenvec %>%
    left_join(statusfile, by = "IID") -> datafile

datafile %>%
  ggplot(., aes(x = PC1, y = PC2, colour = status)) + 
  geom_point(alpha = 0.5) + 
  theme_linedraw()+
  scale_colour_brewer(palette = "Set1", direction = -1)
ggsave("PC1vsPC2_onlyUsersData.png")

datafile %>%
  ggplot(., aes(x = PC1, y = PC3, colour = status)) + 
  geom_point(alpha = 0.5) + 
  theme_linedraw()+
  scale_colour_brewer(palette = "Set1", direction = -1)
ggsave("PC1vsPC3_onlyUsersData.png")

datafile %>%
  ggplot(., aes(x = PC2, y = PC3, colour = status)) + 
  geom_point(alpha = 0.5) + 
  theme_linedraw()+
  scale_colour_brewer(palette = "Set1", direction = -1)
ggsave("PC2vsPC3_onlyUsersData.png")

fancy_plot <- plot_ly(
  datafile,
  x =  ~ PC1,
  y =  ~ PC2,
  z =  ~ PC3,
  size = 3,
  type = "scatter3d",
  mode = "markers",
  color = ~ status,
  colors = "Set1"
)
saveRDS(fancy_plot, "plink_3D_PCA_onlyUsersData.rds")
