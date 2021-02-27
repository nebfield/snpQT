#!/usr/bin/env Rscript

library('tidyverse')
library('cowplot')

# Args
# 1: plink.hwe file path before
# 2: plink.hwe file path after
# 3: threshold
# 4: plot title

args <- commandArgs(trailingOnly = TRUE)

# A function factory for getting integer y-axis values.
integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
        breaks <- floor(pretty(x, n, ...))
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}

read_table(args[[1]]) %>%
    mutate(type = "before") -> before

if (nrow(before) == 0) quit(save = "no", status = 0)

read_table(args[[2]]) %>%
    mutate(type = "after") %>%
    bind_rows(before) %>%
    mutate(type = fct_relevel(type, "before")) -> hwe

fn <- ifelse(args[[4]] != "", "hwe_zoom.png", "hwe_sub.png")

ggplot(hwe, aes(x = P)) +
    geom_histogram() +
    geom_vline(xintercept = as.numeric(args[[3]]), colour = "red") +
    facet_grid(~ type) +
    theme_cowplot() +
	background_grid()+
	panel_border() +
	scale_y_continuous(breaks= integer_breaks())+
    xlab("P-value") +
    ylab("Variant count")+
    ggtitle(paste("Hardy-Weinberg Equilibrium (HWE) ", args[[4]]))
ggsave(fn, , height = 7, width = 10) 
