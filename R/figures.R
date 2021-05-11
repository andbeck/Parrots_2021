##Figures
library("RColorBrewer")

## Immature stage 1:5, constant adult stage 10 ----
source("R/stochastic_projection.R")

stoch_plot + plot_annotation(tag_levels = 'A')
#ggsave("stoch_plot.png", width = 10, height = 6)

stoch_mean_plot + plot_annotation(tag_levels = 'A')
#ggsave("stoch_mean_plot.png", width = 10, height = 6)

stoch_quasi_plot
#ggsave("stoch_quasi_plot.png", width = 8, height = 6)

stoch_elas_plot
#ggsave("stoch_elas_plot.png", width = 8, height = 6)

stoch_elas_plots <- D1_elas_plot + D2_elas_plot + D3_elas_plot + D4_elas_plot + D5_elas_plot
stoch_elas_plots


## Immature stage 1:5, adult 10:6 ----
source("R/duration_10.R")

A10_plot + plot_annotation(tag_levels = 'A')
#ggsave("A10_change_plot.png", width = 10, height = 6)

A10_mean_plot + plot_annotation(tag_levels = 'A')
#gsave("A10_change_mean_plot.png", width = 10, height = 6)

A10_elas_plots <- D1A10_elas_plot + D2A9_elas_plot + D3A8_elas_plot + D4A7_elas_plot + D5A6_elas_plot
A10_elas_plots


## Immature stage 1:5, adult 20:16 ----
source("R/duration_20.R")

A20_plot + plot_annotation(tag_levels = 'A')
#ggsave("A20_plot.png", width = 10, height = 6)

A20_mean_plot + plot_annotation(tag_levels = 'A')
#ggsave("A20_mean_plot.png", width = 10, height = 6)

A20_elas_plots <- D1A20_elas_plot + D2A19_elas_plot + D3A18_elas_plot + D4A17_elas_plot + D5A16_elas_plot
A20_elas_plots


## Immature stage 1:5, adult 30:26 ----
source("R/duration_30.R")

A30_plot + plot_annotation(tag_levels = 'A')
#ggsave("A30_plot.png", width = 10, height = 6)

A30_mean_plot + plot_annotation(tag_levels = 'A')
#ggsave("A30_mean_plot.png", width = 10, height = 6)

A30_elas_plots <- D1A30_elas_plot + D2A29_elas_plot + D3A28_elas_plot + D4A27_elas_plot + D5A26_elas_plot
A30_elas_plots


## Immature stage 1:5, adult 20:26 with 10% reduction in fledgling survival from poaching ----
source("R/poached10_A20.R")

P10_A20_plot + plot_annotation(tag_levels = 'A')
#ggsave("P10_plot.png", width = 10, height = 6)

P10_A20_mean_plot + plot_annotation(tag_levels = 'A')
#ggsave("P10_mean_plot.png", width = 10, height = 6)

P10_A20_elas_plots <- P10D1A20_elas_plot + P10D2A19_elas_plot + P10D3A18_elas_plot + P10D4A17_elas_plot + P10D5A16_elas_plot
P10_A20_elas_plots


## Immature stage 1:5, adult 30:26 with 20% reduction in fledgling survival from poaching ----
source("R/poached20_A20.R")

P20_A20_plot + plot_annotation(tag_levels = 'A')
#ggsave("P20_plot.png", width = 10, height = 6)

P20_A20_mean_plot + plot_annotation(tag_levels = 'A')
#ggsave("P20_mean_plot.png", width = 10, height = 6)

P20_A20_elas_plots <- P20D1A20_elas_plot + P20D2A19_elas_plot + P20D3A18_elas_plot + P20D4A17_elas_plot + P20D5A16_elas_plot
P20_A20_elas_plots


## Combined Plots ----

# Stochastic lambda without poaching
A10_lambda_plot + ylim(0.9, 1.15) + A20_lambda_plot + ylim(0.9, 1.15) + A30_lambda_plot + ylim(0.9, 1.15) + plot_annotation(tag_levels = 'A')
#ggsave("A10-30_lambda_plots.png", width = 14, height = 5)

# Stochastic lambda with poaching 0%, 10%, and 20%
A20_lambda_plot + ylim(0.9, 1.15) + P10_lambda_plot + ylim(0.9, 1.15) + P20_lambda_plot + ylim(0.9, 1.15) + plot_annotation(tag_levels = 'A')
#ggsave("poached_lambda_plots.png", width = 14, height = 5)


# Quasi-extinction without poaching
A10_quasi_plot + theme(legend.position = "NONE") + A20_quasi_plot + theme(legend.position = "NONE") + A30_quasi_plot + plot_annotation(tag_levels = 'A')
#ggsave("A10-30_quasi_plots.png", width = 14, height = 5)

# Quasi-extinction with poaching 0%, 10%, and 20%
A20_quasi_plot + theme(legend.position = "NONE") + P10A20_quasi_plot + theme(legend.position = "NONE") + P20A20_quasi_plot + plot_annotation(tag_levels = 'A')
#ggsave("Poaching_quasi_plots.png", width = 14, height = 5)


# facet grid elasticity plot without poaching
facet_elas_df <- rbind.data.frame(A10_elas_df, A20_elas_df, A30_elas_df)
facet_elas_df$adult <- c(rep("10 adult years", 45), rep("20 adult years", 45), rep("30 adult years", 45))

facet_elas_df$immature <- c(rep("1 immature year", 9), rep("2 immature years", 9), rep("3 immature years", 9), rep("4 immature years", 9), rep("5 immature years", 9), rep("1 immature year", 9), rep("2 immature years", 9), rep("3 immature years", 9), rep("4 immature years", 9), rep("5 immature years", 9), rep("1 immature year", 9), rep("2 immature years", 9), rep("3 immature years", 9), rep("4 immature years", 9), rep("5 immature years", 9))

facet_elas_df <- subset(facet_elas_df, stage!= "g1") 
facet_elas_df <- subset(facet_elas_df, stage!= "m1")

facet_elas_plot <- ggplot(facet_elas_df, aes(x = stage, y= elasticity, fill = stage)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Greens", name = "Matrix \nelement") +
  facet_grid(adult ~ immature)
#ggsave("facet_elas_plot.png", width = 10, height = 6)

# facet grid elasticity plot with poaching 0%, 10%, and 20%
facet_poaching_elas_df <- rbind.data.frame(A20_elas_df, P10A20_elas_df, P20A20_elas_df)

facet_poaching_elas_df$poaching <- c(rep("0% reduction from poaching", 45), rep("10% reduction from poaching", 45), rep("20% reduction from poaching", 45))

facet_poaching_elas_df$immature <- c(rep("1 immature year", 9), rep("2 immature years", 9), rep("3 immature years", 9), rep("4 immature years", 9), rep("5 immature years", 9), rep("1 immature year", 9), rep("2 immature years", 9), rep("3 immature years", 9), rep("4 immature years", 9), rep("5 immature years", 9), rep("1 immature year", 9), rep("2 immature years", 9), rep("3 immature years", 9), rep("4 immature years", 9), rep("5 immature years", 9))

facet_poaching_elas_df <- subset(facet_poaching_elas_df, stage!= "g1") 
facet_poaching_elas_df <- subset(facet_poaching_elas_df, stage!= "m1")

facet_poaching_elas_plot <- ggplot(facet_poaching_elas_df, aes(x = stage, y= elasticity, fill = stage)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Blues", name = "Matrix \nelement") +
  facet_grid(poaching ~ immature)
facet_poaching_elas_plot
#ggsave("facet_poaching_elas_plot.png", width = 10, height = 6)


# Elasticity without poaching
#A10_elas_plot + theme(legend.position = "NONE") + A20_elas_plot + theme(legend.position = "NONE") + A30_elas_plot + plot_annotation(tag_levels = 'A')
#ggsave("A10-30_elasticity_plots.png", width = 14, height = 5)

# Elasticity with poaching 0%, 10%, and 20%
#A20_elas_plot + theme(legend.position = "NONE") + P10A20_elas_plot + theme(legend.position = "NONE") + P20A20_elas_plot + plot_annotation(tag_levels = 'A')
#ggsave("Poaching_elasticity_plots.png", width = 14, height = 5)


