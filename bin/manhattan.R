#!/usr/bin/env Rscript

library(qqman)

# Args
# 1: .logistic file path

args <- commandArgs(trailingOnly = TRUE)

gwas <- read.table(file = args[[1]], sep = "", header = TRUE)

pdf("manhattan.pdf")
manhattan(gwas, main = "Manhattan Plot", cex = 0.5, cex.axis = 0.8)
dev.off()
