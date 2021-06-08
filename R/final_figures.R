##Final figures for report
library("RColorBrewer")
display.brewer.all(colorblindFriendly = T)

source("R/duration_10.R")
source("R/duration_20.R")
source("R/duration_30.R")
source("R/poached10_A20.R")
source("R/poached20_A20.R")
source("R/poached30_A20.R")
source("R/poached40_A20.R")
source("R/poached50_A20.R")
source("R/lambda_stats.R")

### facet lambda plot with lifespans ----

 
ggplot(lambda_df, aes(x = imm_dur, y = lambda)) +
  geom_boxplot(outlier.shape = NA, aes(fill = imm_dur)) +
  facet_wrap(~lifespan) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  labs(x = "Immature stage duration (years)", y = "Stochastic lambda", col = "Lifespan") +
  theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) +
  ylim(0.65, 1.5) +
  geom_signif(comparisons = list(c("1", "2"), c("2", "3"), c("3", "4"), c("4", "5")), map_signif_level = TRUE, y_position = 1.4) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), strip.text.x = element_text(size = 14), strip.background = element_rect(fill = "grey92")) +
  scale_fill_brewer(palette = "YlGnBu")

#ggsave("colour_anova_lambda_plot.png", width = 9, height = 5)


### lambda plot with poaching ----

newX <- expand.grid(
  poaching = seq(from = 0, to = 50, length = 5), 
  imm_dur = c("1 year", "2 years", "3 years", "4 years", "5 years"))
newX

newY <- predict(mod_poach_lambda, newdata = newX, interval = "confidence")
newY

addThese <- data.frame(newX, newY) %>%
  rename(lambda = fit)

#plot
ggplot(poach_lambda_df, aes(x = poaching, y = lambda)) +
  geom_point(aes(col = imm_dur), alpha = 0.1) +
  geom_line(data = addThese, aes(col = imm_dur)) +
  geom_ribbon(data = addThese, aes(ymin = lwr, ymax = upr, fill = imm_dur), alpha = 0.3) +
  theme_classic() +
  ylim(0.75, 1.27) +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  labs(x = "Poaching pressure (% reduction in nestling survival)", y = "Stochastic lambda", col = "Immature stage duration", linetype = "Immature stage duration") +
  scale_fill_brewer(palette = "YlGnBu", name = "Immature stage duration") +
  scale_colour_brewer(palette = "YlGnBu", name = "Immature stage duration") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12))

#ggsave("anova_poaching_plot.png", width = 9, height = 6)


### facet quasi-extintion plot with lifespans ----
quasi_df <- rbind(A10_quasi_df, A20_quasi_df, A30_quasi_df) 
quasi_df$lifespan <- c(rep("10 year lifespan", 500), rep("20 year lifespan", 500), rep("30 year lifespan", 500))

quasi_df$sim <- c(rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100), rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100), rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100))

quasi_plot <- ggplot(quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  scale_colour_discrete(name = "Immature stage duration") +
  facet_wrap(~lifespan) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14), strip.text.x = element_text(size = 16), strip.background = element_rect(fill = "grey92"), panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"))

quasi_plot

#ggsave("facet_quasi_lifespan_plot.png", width = 10, height = 5.25)


### facet quasi-extintion plot with poaching ----

poach_quasi_df <- rbind(A20_quasi_df, P10A20_quasi_df, P20A20_quasi_df, P30A20_quasi_df, P40A20_quasi_df, P50A20_quasi_df) 
poach_quasi_df$poaching <- c(rep("-0% nestling survival", 500), rep("-10% nestling survival", 500), rep("-20% nestling survival", 500), rep("-30% nestling survival", 500), rep("-40% nestling survival", 500), rep("-50% nestling survival", 500))

poach_quasi_df$sim <- c(rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100),rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100), rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100), rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100), rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100), rep("1 year", 100), rep("2 years", 100), rep("3 years", 100), rep("4 years", 100), rep("5 years", 100))

poach_quasi_plot <- ggplot(poach_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  facet_wrap(~poaching) +
  theme(legend.position = "bottom", legend.text = element_text(size = 13), legend.title = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 13), strip.text.x = element_text(size = 15), strip.background = element_rect(fill = "grey92"), panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) + 
  scale_colour_discrete(name = "Immature stage duration")
poach_quasi_plot

#ggsave("facet_quasi_poaching_plot.png", width = 9, height = 8)




### facet elasticity plot with lifespans ----

new_elas_df <- rbind.data.frame(A10_elas_df, A20_elas_df, A30_elas_df)
new_elas_df$adult <- c(rep("10 year lifespan", 45), rep("20 year lifespan", 45), rep("30 year lifespan", 45))

new_elas_df$immature <- c(rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9))

new_elas_df <- subset(new_elas_df, stage!= "g1") 
new_elas_df <- subset(new_elas_df, stage!= "m1")

new_elas_df$elasticity <- new_elas_df$elasticity %>% round(3) 

new_elas_df <- new_elas_df %>% filter(elasticity!=0)

new_elas_plot <- ggplot(new_elas_df, aes(x = immature, y= elasticity, fill = stage)) + 
  geom_col() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) +
  labs(x = "Immature stage duration (years)", y = "Stochastic elasticity") +
  facet_wrap("adult") +
  scale_fill_brewer(palette = "YlOrRd", name = "Vital rate", labels = c("Immature growth", "Immature reproduction", "Adult reproduction", "Fledgling survival", "Immature survival", "Adult survival")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), strip.text.x = element_text(size = 14), strip.background = element_rect(fill = "grey92")) #+geom_text(aes(label = elasticity), size = 4, position = position_stack(vjust = 0.5)) 

new_elas_plot

#ggsave("new_facet_elas_plot.png", width = 9, height = 6.4)


### facet elasticity plot with poaching 0-50% ----

new_poaching_elas_df <- rbind.data.frame(A20_elas_df, P10A20_elas_df, P20A20_elas_df, P30A20_elas_df, P40A20_elas_df, P50A20_elas_df)

new_poaching_elas_df$poaching <- c(rep("-0% nestling survival", 45), rep("-10% nestling survival", 45), rep("-20%  nestling survival", 45), rep("-30% nestling survival", 45), rep("-40% nestling survival", 45), rep("-50% nestling survival", 45))

new_poaching_elas_df$immature <- c(rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9))

new_poaching_elas_df <- subset(new_poaching_elas_df, stage!= "g1") 
new_poaching_elas_df <- subset(new_poaching_elas_df, stage!= "m1")

new_poaching_elas_df$elasticity <- new_poaching_elas_df$elasticity %>% round(3)

new_poaching_elas_df <- new_poaching_elas_df %>% filter(elasticity!=0)

new_poaching_elas_plot <- ggplot(new_poaching_elas_df, aes(x = immature, y= elasticity, fill = stage)) + 
  geom_col() +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) +
  labs(x = "Immature stage duration (years)", y = "Stochastic elasticity") +
  facet_wrap("poaching") +
  scale_fill_brewer(palette = "YlOrRd", name = "Vital rate", labels = c("Immature growth", "Immature reproduction", "Adult reproduction", "Fledgling survival", "Immature survival", "Adult survival")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), strip.text.x = element_text(size = 14), strip.background = element_rect(fill = "grey92")) #+geom_text(aes(label = elasticity), size = 4, position = position_stack(vjust = 0.5)) 
new_poaching_elas_plot

#ggsave("new_poach_elas_plot.png", width = 9, height = 11.5)
