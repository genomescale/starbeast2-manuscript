library(ggplot2)
library(scales)
library(cowplot)

branch_length_accuracy = read.csv("branch_length_accuracy.csv")
branch_length_accuracy = subset(branch_length_accuracy, sequence == "haplotype")

branch_length_accuracy$method = factor(branch_length_accuracy$method, levels = c("beast", "starbeast1", "starbeast2"))
levels(branch_length_accuracy$method) = c("Concatenation", "*BEAST", "StarBEAST2")

branch_length_accuracy$n_loci = factor(branch_length_accuracy$n_loci, levels = c("1x", "2x", "5x"))
levels(branch_length_accuracy$n_loci) = c("\n26 loci", "\n52 loci", "\n130 loci")

branch_length_accuracy$clock = factor(branch_length_accuracy$clock, levels = c("strict", "gt-ucln", "st-ucln"))
levels(branch_length_accuracy$clock) = c("strict clock", "GT-UCLN", "ST-UCLN")

branch_length_accuracy$analysis = factor(paste(branch_length_accuracy$method, branch_length_accuracy$n_loci),
levels = paste(rep(levels(branch_length_accuracy$method), c(3, 3, 3)), levels(branch_length_accuracy$n_loci)))

tip_error_plot = ggplot(branch_length_accuracy, aes(x = analysis, y = tip_length_error, fill = clock)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Tip length error") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

tip_bias_plot = ggplot(branch_length_accuracy, aes(x = analysis, y = tip_length_bias, fill = clock)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Total tip length bias") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0), axis.title.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

tip_hpd_plot = ggplot(branch_length_accuracy, aes(x = analysis, y = tip_truth_in_hpd, fill = clock)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  ylab("Tip lengths in HPD") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

internal_error_plot = ggplot(branch_length_accuracy, aes(x = analysis, y = internal_length_error, fill = clock)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Internal length error") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

internal_bias_plot = ggplot(branch_length_accuracy, aes(x = analysis, y = internal_length_bias, fill = clock)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Total internal length bias") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0), axis.title.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

internal_hpd_plot = ggplot(branch_length_accuracy, aes(x = analysis, y = internal_truth_in_hpd, fill = clock)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  ylab("Internal lengths in HPD") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

prow = plot_grid(tip_hpd_plot + theme(legend.position = "none"),
  internal_hpd_plot + theme(legend.position = "none"),
  tip_error_plot + theme(legend.position = "none"),
  internal_error_plot + theme(legend.position = "none"),
  tip_bias_plot + theme(legend.position = "none"),
  internal_bias_plot + theme(legend.position = "none"),
  labels = c("A", "B", "C", "D", "E", "F"), nrow = 3, rel_heights = c(0.8, 0.8, 1.1), align = "v", vjust = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
legend = get_legend(tip_error_plot)
plot_grid(legend, prow, ncol = 1, rel_heights = c(0.07, 1))

ggsave("branch_length_accuracy_phased.pdf", units = "mm", height = 250, width = 200)
ggsave("branch_length_accuracy_phased.png", units = "mm", height = 250, width = 200)
