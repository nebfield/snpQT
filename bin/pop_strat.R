#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: racefile 

args <- commandArgs(trailingOnly = TRUE)

eigenvals <- readr::read_table(args[[1]], skip = 1, col_names = FALSE) %>%
  mutate(id = gsub(".*:","",X1)) # colons mess up matching

df <- readr::read_delim(args[[2]], delim = " ") %>%
  left_join(eigenvals, by = c("IID" = "id")) %>%
  drop_na(X1) %>%
  select(-race)

colnames(df) <-
  c(
    "FID",
    "IID",
    "id",
    "PCA1",
    "PCA2",
    "PCA3",
    "PCA4",
    "PCA5",
    "PCA6",
    "PCA7",
    "PCA8",
    "PCA9",
    "PCA10",
    "race" 
  )

plot_pca <- function(df, ax1, ax2) {
  # {{ }} tidy evaluation for column names
  ggplot(df, aes(x = {{ ax1 }} , y = {{ ax2 }}, colour = race)) +
    geom_point() +
    theme_linedraw()
}

plot_pca(df, PCA1, PCA2)
ggsave("pca1vspca2.png")
plot_pca(df, PCA1, PCA3)
ggsave("pca1vspca3.png")
plot_pca(df, PCA2, PCA3)
ggsave("pca2vspca3.png")

fancy_plot <- plot_ly(
  df,
  x =  ~ PCA1,
  y =  ~ PCA2,
  z =  ~ PCA3,
  size = 3,
  type = "scatter3d",
  mode = "markers",
  color =  ~ race
)

saveWidget(fancy_plot, "popStrat.html", selfcontained = T)
