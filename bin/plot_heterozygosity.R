#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: .het file path

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  mutate(HET_RATE = (N.NM. - O.HOM.) / N.NM.) -> het
  
ggplot(het, aes(x = IID, y = HET_RATE)) + 
  geom_point() + 
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("Sample") +
  ylab("Heterozygosity rate") +
  ggtitle("Heterozygosity rate per sample")
ggsave("heterozygosity_rate.png", device = "png")


het %>%
  filter(HET_RATE < mean(HET_RATE) - 3 * sd(HET_RATE) |
	            HET_RATE > mean(HET_RATE) + 3 * sd(HET_RATE)) -> het_fail
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)

write.table(het_fail, "het_failed_samples.txt", row.names=FALSE, quote = FALSE)
