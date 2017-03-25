library(scales)
library(ggplot2)
library(cowplot)

simulation_rates = read.csv("21-species/minimum-ess_rates.csv")
pseudacris_rates = read.csv("pseudacris/minimum-ess_rates.csv")

ess_rates = rbind(simulation_rates, pseudacris_rates)

ess_rates$popsize_integration = factor(ess_rates$popsize_integration, levels = c("mcmc", "analytical"))
ess_rates$clock = factor(ess_rates$clock, levels = c("strict", "gt-ucln", "st-ucln"))
ess_rates$operators = factor(ess_rates$operators, levels = c("neither", "topology_changing", "height_changing", "both"))

ess_rates$coordinated_height_operators = ess_rates$operators == "both" | ess_rates$operators == "height_changing"
ess_rates$coordinated_topology_operators = ess_rates$operators == "both" | ess_rates$operators == "topology_changing"

starbeast2_rates = subset(ess_rates, method == "starbeast2" & loci == "1x")

sink("minimum-ess_per_hour-starbeast2.txt")
summary(lm(log(ess_per_hour) ~ coordinated_topology_operators + coordinated_height_operators + popsize_integration, data = subset(starbeast2_rates, clock == "strict")))
summary(lm(log(ess_per_hour) ~ coordinated_topology_operators + coordinated_height_operators + popsize_integration, data = subset(starbeast2_rates, clock == "gt-ucln")))
summary(lm(log(ess_per_hour) ~ coordinated_topology_operators + coordinated_height_operators + popsize_integration, data = subset(starbeast2_rates, clock == "st-ucln")))
sink()

levels(starbeast2_rates$popsize_integration) = c("MCMC population size integration", "Analytical population size integration")
levels(starbeast2_rates$clock) = c("Strict clock", "Gene tree relaxed clocks", "Species tree relaxed clock")
levels(starbeast2_rates$operators) = c("Neither", "Topology changing", "Height changing", "Both")

ggplot(starbeast2_rates, aes(y = ess_per_hour, x = operators, fill = popsize_integration)) +
	geom_boxplot() +
	xlab("Coordinated MCMC operators") +
	ylab("ESS per hour") +
	scale_y_log10(breaks = 2^seq(-20, 20, 2)) +
	facet_grid(dataset ~ clock) +
	geom_vline (xintercept = seq(1.5, length(unique(starbeast2_rates$operators))-0.5, 1), lwd = 0.5, colour = "white") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		legend.position = "top",
		legend.title = element_blank(),
		panel.grid.major.x = element_blank()) +
	background_grid(major = "y", minor = "none")

ggsave("minimum-ess_per_hour-starbeast2.pdf", units = "mm", height = 150, width = 200)
ggsave("minimum-ess_per_hour-starbeast2.png", units = "mm", height = 150, width = 200)
