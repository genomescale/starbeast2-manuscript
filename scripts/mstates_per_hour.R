library(ggplot2)
library(cowplot)
library(scales)

sim_ess_table = read.csv("21-species/minimum-ess_rates.csv")
emp_ess_table = read.csv("pseudacris/minimum-ess_rates.csv")

ess_rates = rbind(sim_ess_table, emp_ess_table)

ess_rates$mstates_rate = ess_rates$ess_per_hour / ess_rates$ess_per_mstates

ess_rates$config = as.factor(ess_rates$config)
ess_rates$popsize_integration = factor(ess_rates$popsize_integration, levels = c("mcmc", "analytical"))
ess_rates$clock = factor(ess_rates$clock, levels = c("strict", "gt-ucln", "st-ucln"))
ess_rates$operators = factor(ess_rates$operators, levels = c("neither", "topology_changing", "height_changing", "both"))

starbeast2_rates = subset(ess_rates, method == "starbeast2" & loci == "1x")

levels(starbeast2_rates$popsize_integration) = c("MCMC population size integration", "Analytical population size integration")
levels(starbeast2_rates$clock) = c("Strict clock", "Gene tree relaxed clocks", "Species tree relaxed clock")
levels(starbeast2_rates$operators) = c("Neither", "Topology changing", "Height changing", "Both")

ggplot(starbeast2_rates, aes(y = mstates_rate, x = operators, fill = popsize_integration)) +
	geom_boxplot(lwd = 0.3) +
	xlab("Coordinated MCMC operators") +
	ylab("Million states per hour") +
	scale_y_continuous(limits = c(0, 25)) +
	facet_grid(dataset ~ clock) +
	geom_vline (xintercept = seq(1.5, length(unique(starbeast2_rates$operators))-0.5, 1), lwd = 0.5, colour = "white") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		legend.position = "top",
		legend.title = element_blank(),
		panel.grid.major.x = element_blank()) +
	background_grid(major = "y", minor = "none")

ggsave("mstates_per_hour.pdf", units = "mm", height = 150, width = 200)
ggsave("mstates_per_hour.png", units = "mm", height = 150, width = 200)
