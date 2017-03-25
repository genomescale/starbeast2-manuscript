library(cowplot)
library(ggplot2)
library(scales)

branch_rates = read.csv("branch_rates.csv")
branch_rates$analysis = factor(branch_rates$analysis, levels = c("u", "v"))
branch_rates$branch = factor(branch_rates$branch, levels = c("s0", "s1", "s2", "s3", "s4", "s0,s1", "(s0,s1),s2", "((s0,s1),s2),s3"))
levels(branch_rates$analysis) = c("StarBEAST2", "Concatenation")
levels(branch_rates$branch) = c("A", "B", "C", "D", "E", "AB", "ABC", "ABCD")

ggplot(branch_rates, aes(color = analysis, x = rate_mean)) +
  geom_freqpoly(aes(y = 16*(..count..)/sum(..count..)), binwidth = 0.1) +
  # geom_vline(xintercept = 0) +
  scale_y_continuous(labels = percent) +
  scale_x_continuous(breaks = seq(0.0, 2.0, 0.2), limits = c(0.5, 1.5)) +
  labs(x = "Estimated branch rate", y = "Proportion of all branches") +
  facet_wrap(~ branch, nrow = 2) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5)) +
  background_grid(major = "xy", minor = "none")
ggsave("spils.pdf", units = "mm", height = 120, width = 200)
ggsave("spils.png", units = "mm", height = 120, width = 200)

branch_rates$length_diff = branch_rates$estimated_length - branch_rates$true_length
branch_rates$rate_diff = branch_rates$rate_mean - 1.0
ggplot(branch_rates, aes(x = rate_diff, y = length_diff, colour = analysis)) +
  geom_point(size = 0.1) +
  geom_point(x = 0, y = 0, colour = "black", shape = 3, size = 2) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2), limits = c(-0.5, 0.5)) +
  labs(y = "Branch length deviation", x = "Per-branch substitution rate deviation") +
  facet_wrap(~ branch, nrow = 2) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5)) +
  background_grid(major = "xy", minor = "none")
ggsave("scatter.pdf", units = "mm", height = 120, width = 200)
ggsave("scatter.png", units = "mm", height = 120, width = 200)

concat_rates = subset(branch_rates, analysis == "Concatenation")

concat_rates$length_diff = concat_rates$estimated_length - concat_rates$true_length
concat_rates$rate_diff = concat_rates$rate_mean - 1.0
ggplot(concat_rates, aes(x = rate_diff, y = length_diff)) +
  geom_point(size = 0.1, colour = "#00BFC4") +
  geom_point(x = 0, y = 0, colour = "black", shape = 3, size = 2) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2), limits = c(-0.5, 0.5)) +
  labs(y = "Branch length deviation", x = "Per-branch substitution rate deviation") +
  facet_wrap(~ branch, nrow = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5)) +
  background_grid(major = "xy", minor = "none")
ggsave("concat_spils.pdf", units = "mm", height = 120, width = 200)
ggsave("concat_spils.png", units = "mm", height = 120, width = 200)
