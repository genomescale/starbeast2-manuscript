library(ggplot2)
library(cowplot)
library(scales)

branch_rates = read.csv("branch_rates.csv")
branch_rates = subset(branch_rates, clock == "st-ucln")
branch_rates = subset(branch_rates, method != "beast" | sequence == "diplotype")
branch_rates = subset(branch_rates, estimated_rate >= 0.0)

branch_rates$method = factor(branch_rates$method, levels = c("starbeast2", "beast"))
levels(branch_rates$method) = c("StarBEAST2", "Concatenation")

branch_rates$n_loci = factor(branch_rates$n_loci, levels = c("1x", "2x", "5x"))
levels(branch_rates$n_loci) = c("26 loci", "52 loci", "130 loci")

branch_rates$analysis = factor(paste(branch_rates$method, branch_rates$n_loci),
                                   levels = paste(rep(levels(branch_rates$method), c(3, 3)), levels(branch_rates$n_loci)))

concat_1x_rates = subset(branch_rates, analysis == "Concatenation 26 loci")
concat_5x_rates = subset(branch_rates, analysis == "Concatenation 130 loci")
starbeast2_1x_rates = subset(branch_rates, analysis == "StarBEAST2 26 loci")
starbeast2_2x_rates = subset(branch_rates, analysis == "StarBEAST2 52 loci")

branch_rates$rate_recovered = ifelse(branch_rates$rate > branch_rates$estimated_rate_low & branch_rates$rate < branch_rates$estimated_rate_high, 1, 0)

sink("branch_rates_unphased.txt")
aggregate(rate_recovered ~ n_loci * method, data = branch_rates, mean)
aggregate(rate_recovered ~ n_loci * method, data = branch_rates, sd)
summary(lm(log(estimated_rate) ~ log(rate), data = concat_1x_rates))
summary(lm(log(estimated_rate) ~ log(rate), data = concat_5x_rates))
summary(lm(log(estimated_rate) ~ log(rate), data = starbeast2_1x_rates))
summary(lm(log(estimated_rate) ~ log(rate), data = starbeast2_2x_rates))
sink()

ggplot(branch_rates, aes(x = rate, y = estimated_rate)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0.25, 2.75)) +
  scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0.25, 2.0)) +
  labs(y = "Estimated branch rate", x = "Simulated branch rate") +
  facet_grid(analysis ~ .) +
  background_grid(major = "xy", minor = "none")
ggsave("branch_rates_unphased.pdf", units = "mm", height = 210, width = 95)
ggsave("branch_rates_unphased.png", units = "mm", height = 210, width = 95)
