library(ggplot2)
library(scales)
library(cowplot)

topology_accuracy = read.csv("topology_accuracy.csv")
topology_accuracy = subset(topology_accuracy, method != "beast" | sequence == "diplotype")

topology_accuracy$method = factor(topology_accuracy$method, levels = c("beast", "starbeast1", "starbeast2"))
levels(topology_accuracy$method) = c("Concatenation", "*BEAST", "StarBEAST2")

topology_accuracy$n_loci = factor(topology_accuracy$n_loci, levels = c("1x", "2x", "5x"))
levels(topology_accuracy$n_loci) = c("\n26 loci", "\n52 loci", "\n130 loci")

topology_accuracy$clock = factor(topology_accuracy$clock, levels = c("strict", "gt-ucln", "st-ucln"))
levels(topology_accuracy$clock) = c("strict clock", "GT-UCLN", "ST-UCLN")

topology_accuracy$analysis = factor(paste(topology_accuracy$method, topology_accuracy$n_loci),
levels = paste(rep(levels(topology_accuracy$method), c(3, 3, 3)), levels(topology_accuracy$n_loci)))

dodge = position_dodge(width = 0.9)

mcc_rf_plot = ggplot(topology_accuracy, aes(x = analysis, y = mcc_rf_mean, color = clock)) +
  geom_point(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = mcc_rf_low, ymax = mcc_rf_high), width = 0.5, position = dodge) +
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  ylab("Mean MCC RF error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0), axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

in_cs_plot = ggplot(topology_accuracy, aes(x = analysis, y = in_cs_prop, color = clock)) +
  geom_point(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = in_cs_low, ymax = in_cs_high), width = 0.5, position = dodge) +
  ylab("True topology in credible set") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

plot_grid(in_cs_plot, mcc_rf_plot, labels = c("A", "B"), nrow = 2, rel_heights = c(1.0, 1.1), align = "v", vjust = c(4.0, 0.5))
ggsave("topology_accuracy_unphased.pdf", units = "mm", height = 200, width = 120)
ggsave("topology_accuracy_unphased.png", units = "mm", height = 200, width = 120)
