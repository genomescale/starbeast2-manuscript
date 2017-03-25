library(scales)
library(ggplot2)
library(cowplot)

simulation_rates = read.csv("21-species/minimum-ess_rates.csv")
pseudacris_rates = read.csv("pseudacris/minimum-ess_rates.csv")
crocidura_rates = read.csv("crocidura/minimum-ess_rates.csv")

simulation_rates = subset(simulation_rates, sequence == "haploid" & (method == "starbeast1" | method == "beast" | (operators == "height_changing" & popsize_integration == "analytical")))
pseudacris_rates = subset(pseudacris_rates, sequence == "haploid" & (method == "starbeast1" | method == "beast" | (operators == "height_changing" & popsize_integration == "analytical")))
crocidura_rates = subset(crocidura_rates, sequence == "haploid" & (method == "starbeast1" | method == "beast" | (operators == "height_changing" & popsize_integration == "analytical")))

simulation_rates$clock = factor(simulation_rates$clock, levels = c("strict", "gt-ucln", "st-ucln"))
pseudacris_rates$clock = factor(pseudacris_rates$clock, levels = c("strict", "gt-ucln", "st-ucln"))
crocidura_rates$clock = factor(crocidura_rates$clock, levels = c("strict", "gt-ucln", "st-ucln"))
levels(simulation_rates$clock) = c("Strict clock", "GT-UCLN", "ST-UCLN")
levels(pseudacris_rates$clock) = c("Strict clock", "GT-UCLN", "ST-UCLN")
levels(crocidura_rates$clock) = c("Strict clock", "GT-UCLN", "ST-UCLN")

simulation_rates$method = factor(simulation_rates$method, levels = c("starbeast1", "starbeast2", "beast"))
pseudacris_rates$method = factor(pseudacris_rates$method, levels = c("starbeast1", "starbeast2", "beast"))
crocidura_rates$method = factor(crocidura_rates$method, levels = c("starbeast1", "starbeast2", "beast"))
levels(simulation_rates$method) = c("*BEAST", "StarBEAST2", "Concatenation")
levels(pseudacris_rates$method) = c("*BEAST", "StarBEAST2", "Concatenation")
levels(crocidura_rates$method) = c("*BEAST", "StarBEAST2", "Concatenation")

simulation_rates$loci = factor(simulation_rates$loci, levels = c("1x", "2x", "5x"))
pseudacris_rates$loci = factor(pseudacris_rates$loci, levels = c("1x"))
crocidura_rates$loci = factor(crocidura_rates$loci, levels = c("1x", "2x"))

levels(simulation_rates$loci) = c("26 loci", "52 loci", "130 loci")
levels(pseudacris_rates$loci) = c("26 loci")
levels(crocidura_rates$loci) = c("50 loci", "100 loci")

simulation_rates$category = paste0(simulation_rates$method, "
", simulation_rates$loci)
pseudacris_rates$category = paste0(pseudacris_rates$method, "
", pseudacris_rates$loci)
crocidura_rates$category = paste0(crocidura_rates$method, "
", crocidura_rates$loci)

simulation_rates$category = factor(simulation_rates$category, levels = paste0(rep(levels(simulation_rates$method), c(3, 3, 3)), "
", levels(simulation_rates$loci)))
pseudacris_rates$category = factor(pseudacris_rates$category, levels = paste0(levels(pseudacris_rates$method), "
", levels(pseudacris_rates$loci)))
crocidura_rates$category = factor(crocidura_rates$category, levels = paste0(rep(levels(crocidura_rates$method), c(2, 2, 2)), "
", levels(crocidura_rates$loci)))

simulation_plot = ggplot(simulation_rates, aes(y = ess_per_hour, x = category, fill = clock)) +
	geom_boxplot() +
	ylab("ESS per hour") +
	scale_y_log10(limits = c(min(simulation_rates$ess_per_hour), 1000.0), breaks = 10^seq(-10, 10, 1), labels = comma) +
	facet_grid(. ~ dataset) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		axis.title.x = element_blank(),
		plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "lines"),
		legend.position = "top",
		legend.title = element_blank()) +
	background_grid(major = "y", minor = "none")

pseudacris_plot = ggplot(pseudacris_rates, aes(y = ess_per_hour, x = category, fill = clock)) +
	geom_boxplot() +
	scale_y_log10(limits = c(min(simulation_rates$ess_per_hour), 1000.0), breaks = 10^seq(-10, 10, 1), labels = comma) +
	facet_grid(. ~ dataset) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		axis.text.y = element_blank(),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.y = element_blank(),
		plot.margin = unit(c(0.5, 0.1, 0.5, 0.1), "lines"),
		legend.position = "top",
		legend.title = element_blank()) +
	background_grid(major = "y", minor = "none")

crocidura_plot = ggplot(crocidura_rates, aes(y = ess_per_hour, x = category, fill = clock)) +
	geom_boxplot() +
	scale_y_log10(limits = c(min(simulation_rates$ess_per_hour), 1000.0), breaks = 10^seq(-10, 10, 1), labels = comma) +
	facet_grid(. ~ dataset) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		axis.text.y = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.y = element_blank(),
		plot.margin = unit(c(0.5, 0.5, 0.5, 0.1), "lines"),
		legend.position = "top",
		legend.title = element_blank()) +
	background_grid(major = "y", minor = "none")

prow = plot_grid(simulation_plot + theme(legend.position = "none"), pseudacris_plot + theme(legend.position = "none"), crocidura_plot + theme(legend.position = "none"), nrow = 1, align = "h", rel_widths = c(1.0, 0.35, 0.5))
legend = get_legend(simulation_plot)
plot_grid(legend, prow, ncol = 1, rel_heights = c(0.15, 1))

ggsave("minimum-ess_per_hour-comparison.pdf", units = "mm", height = 100, width = 200)
ggsave("minimum-ess_per_hour-comparison.png", units = "mm", height = 100, width = 200)
