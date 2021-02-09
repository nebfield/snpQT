#!/usr/bin/env Rscript

library("tidyverse")

# Args
# 1: .het file path
# 2: .het threshold
# 3: before or after threshold

args <- commandArgs(trailingOnly = TRUE)

read.table(args[[1]], header = TRUE) %>%
  mutate(HET_RATE = (N.NM. - O.HOM.) / N.NM.) -> het

het %>%
  summarise(threshold = mean(HET_RATE) - 3 * sd(HET_RATE)) %>%
  pull(threshold) -> threshold
  
ggplot(het, aes(x = IID, y = HET_RATE)) + 
  geom_point() + 
  geom_vline(xintercept = threshold, colour = "red")+
  geom_vline(xintercept = -threshold, colour = "red")+
  theme_linedraw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("Sample") +
  ylab("Heterozygosity rate") +
  ggtitle("Heterozygosity rate per sample")
ggsave(paste0("heterozygosity_rate_",args[[3]],".png", device = "png")


het %>%
  filter(HET_RATE < mean(HET_RATE) - 3 * sd(HET_RATE) |
	            HET_RATE > mean(HET_RATE) + 3 * sd(HET_RATE)) -> het_fail
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)

write.table(het_fail, "het_failed_samples.txt", row.names=FALSE, quote = FALSE)
