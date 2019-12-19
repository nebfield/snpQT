#!/usr/bin/env Rscript

library("tidyverse")

as_tibble(read.table("plink.hwe", header = TRUE)) %>%
    ggplot(.,  aes(x = P)) +
    geom_histogram() +
    xlab("P") +
    ylab("Count") +
    theme_classic() + 
    ggtitle("HWE histogram")
ggsave("hwe.png")

as_tibble(read.table("plinkzoomhwe.hwe", header = TRUE)) %>%
    ggplot(.,  aes(x = P)) +
    geom_histogram() +
    xlab("P") +
    ylab("Count") +
    theme_classic() + 
    ggtitle("HWE histogram: strongly deviating SNPs only")
ggsave("hwe_zoom.png")
