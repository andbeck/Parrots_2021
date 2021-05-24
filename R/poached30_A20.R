# YSA Stochastic Model of Population Growth under Poaching Pressure of 30% reduced fledgling survival
# With changes in immature stage duration of 1:5 years.
# With adult stage duration of 20:16 years.
# With imputed survival rates.
# Density independent

#### Libraries ----
library(popbio)
library(tidyverse)
library(patchwork)


#### Functions ----

## Matrix model function
source("R/make_projection_matrix.R")

## Stochastic population growth function
source("R/stochastic_proj.R")


#### YSA Data ----

## YSA breeding biology data 2006-2014 from Bonaire
source("R/YSA_life_history_data.R")

# Mean fecundity
fecundity <- c(0, 0, 1.6*total_summary$mean_hatch[1]*(total_summary$mean_nestling_surv[1]*0.7))

# Mean survival (0.73 is from Salinas-Melgoza & Renton 2007, 0.838 is survival from imputation)
survival <- c(0.73, 0.838, 0.838)

# Current population is estimated around 1000 individuals. 1:1 sex ratio means female population is 500
Nc <- 500

# Time to project to
time <- 100

#### YSA Simulated Vital Rates for LSA ----

set.seed(2021)

# Number of simulations
n_sim <- 1000

# Fledgling survival
s1 <- sapply(1:n_sim, function(x) betaval(0.73, 0.2))

# Immature survival
s2 <- sapply(1:n_sim, function(x) betaval(0.838, 0.051))

# Adult survival
s3 <- sapply(1:n_sim, function(x) betaval(0.838, 0.051))

# Fecundity
m3 <- rlnorm(n = n_sim, 
             log(1.6*total_summary$mean_hatch[1]*(total_summary$mean_nestling_surv[1]*0.7)), 
             log(1.01)) #replaced sd with small value for log

## Create lists of survival and fecundity

# Survival 
survival_df <- data.frame(s1, s2, s3)
colnames(survival_df)<- c()
survival_list <- asplit(survival_df, 1)

# Fecundity 
fecundity_df <- data.frame(0, 0, m3)
colnames(fecundity_df)<- c()
fecundity_list <- asplit(fecundity_df, 1)



#### LSA for Immature Duration of 1 and Adult Duration 0f 18 ----

## Stage duration
duration <- c(1, 1, 18)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
P30D1A18_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P30D1A18_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P30D1A18_matrices[[i]] <- mpm
}

head(P30D1A18_matrices)

## Repeat Stochastic Population Growth 

P30D1A18_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P30D1A18_matrices, n = P30D1A18_n0, time = time)
  P30D1A18_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P30D1A18_total_pop <- lapply(P30D1A18_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P30D1A18_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P30D1A18_total_pop[[i]])
  P30D1A18_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P30D1A18_plot_data <- bind_rows(P30D1A18_df_plots, .id = "id")

# Plot projection
P30D1A18_plot <- ggplot(P30D1A18_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P30D1A18_mean_plot_data <- P30D1A18_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P30D1A18_plot_pred <- P30D1A18_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P30D1A18_mean_plot <- ggplot(P30D1A18_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P30D1A18_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P30D1A18_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P30D1A18_total_pop[[i]][time]
  P30D1A18_pop_sizes[i] <- ms
}

#mean pop size
P30D1A18_pop_mean <- mean(P30D1A18_pop_sizes)

# standard deviation pop size
P30D1A18_pop_sd <- sd(P30D1A18_pop_sizes)

# standard error pop size
P30D1A18_pop_se <- sd(P30D1A18_pop_sizes)/sqrt(length(P30D1A18_pop_sizes))


#### Calculate Stochastic Growth Rate

P30D1A18_lambda_s <- stoch.growth.rate(P30D1A18_matrices, prob = NULL, maxt = time,
                                       verbose = TRUE)

#convert from log
P30D1A18_lambda_s$approx <- exp(P30D1A18_lambda_s$approx)
P30D1A18_lambda_s$sim <- exp(P30D1A18_lambda_s$sim)
P30D1A18_lambda_s$sim.CI <- exp(P30D1A18_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P30D1A18_quasi <- stoch.quasi.ext(P30D1A18_matrices, n0= P30D1A18_n0, Nx = 50, tmax = time, maxruns = 1,
                                  nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P30D1A18_quasi_df <- data.frame(P30D1A18_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P30D1A18_quasi_plot <- ggplot(P30D1A18_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P30D1A18_sens <- stoch.sens(P30D1A18_matrices, tlimit=time)

P30D1A18_elas <- P30D1A18_sens$elasticities

P30D1A18_elas_v <- c(P30D1A18_elas[1,1], P30D1A18_elas[1,2], P30D1A18_elas[1,3], P30D1A18_elas[2,1], P30D1A18_elas[2,2], P30D1A18_elas[2,3], P30D1A18_elas[3,1], P30D1A18_elas[3,2], P30D1A18_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

P30D1A18_elas_df <- data.frame(P30D1A18_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P30D1A18_elas_plot <- ggplot(P30D1A18_elas_df, aes(x = stage, y= P30D1A18_elas_v)) + 
  labs(x = "Vital rate", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 2 and Adult Duration of 17 ----

## Stage duration
duration <- c(1, 2, 17)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

## Initial population vector estimated from stable stage distribution
P30D2A17_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P30D2A17_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P30D2A17_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

P30D2A17_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P30D2A17_matrices, n = P30D2A17_n0, time = time) 
  P30D2A17_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P30D2A17_total_pop <- lapply(P30D2A17_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P30D2A17_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P30D2A17_total_pop[[i]])
  P30D2A17_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P30D2A17_plot_data <- bind_rows(P30D2A17_df_plots, .id = "id")

# Plot projection
P30D2A17_plot <- ggplot(P30D2A17_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P30D2A17_mean_plot_data <- P30D2A17_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P30D2A17_plot_pred <- P30D2A17_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P30D2A17_mean_plot <- ggplot(P30D2A17_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P30D2A17_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P30D2A17_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P30D2A17_total_pop[[i]][time]
  P30D2A17_pop_sizes[i] <- ms
}

# mean pop size
P30D2A17_pop_mean <- mean(P30D2A17_pop_sizes)

# standard deviation pop size
P30D2A17_pop_sd <- sd(P30D2A17_pop_sizes)

#standard error pop size
P30D2A17_pop_se <- sd(P30D2A17_pop_sizes)/sqrt(length(P30D2A17_pop_sizes))


#### Calculate Stochastic Growth Rate

P30D2A17_lambda_s <- stoch.growth.rate(P30D2A17_matrices, prob = NULL, maxt = time,
                                       verbose = TRUE)

#convert from log
P30D2A17_lambda_s$approx <- exp(P30D2A17_lambda_s$approx)
P30D2A17_lambda_s$sim <- exp(P30D2A17_lambda_s$sim)
P30D2A17_lambda_s$sim.CI <- exp(P30D2A17_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P30D2A17_quasi <- stoch.quasi.ext(P30D2A17_matrices, n0= P30D2A17_n0, Nx = 50, tmax = time, maxruns = 1,
                                  nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P30D2A17_quasi_df <- data.frame(P30D2A17_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P30D2A17_quasi_plot <- ggplot(P30D2A17_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P30D2A17_sens <- stoch.sens(P30D2A17_matrices, tlimit=time)

P30D2A17_elas <- P30D2A17_sens$elasticities

P30D2A17_elas_v <- c(P30D2A17_elas[1,1], P30D2A17_elas[1,2], P30D2A17_elas[1,3], P30D2A17_elas[2,1], P30D2A17_elas[2,2], P30D2A17_elas[2,3], P30D2A17_elas[3,1], P30D2A17_elas[3,2], P30D2A17_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

P30D2A17_elas_df <- data.frame(P30D2A17_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P30D2A17_elas_plot <- ggplot(P30D2A17_elas_df, aes(x = stage, y= P30D2A17_elas_v)) + 
  labs(x = "Vital rate", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 3 and Adult Duration of 16 ----

## Stage duration
duration <- c(1, 3, 16)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
P30D3A16_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P30D3A16_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P30D3A16_matrices[[i]] <- mpm
}


##Repeat Stochastic Population Growth 

P30D3A16_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P30D3A16_matrices, n = P30D3A16_n0, time = time) 
  P30D3A16_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P30D3A16_total_pop <- lapply(P30D3A16_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P30D3A16_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P30D3A16_total_pop[[i]])
  P30D3A16_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P30D3A16_plot_data <- bind_rows(P30D3A16_df_plots, .id = "id")

# Plot projection
P30D3A16_plot <- ggplot(P30D3A16_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P30D3A16_mean_plot_data <- P30D3A16_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P30D3A16_plot_pred <- P30D3A16_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P30D3A16_mean_plot <- ggplot(P30D3A16_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P30D3A16_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P30D3A16_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P30D3A16_total_pop[[i]][time]
  P30D3A16_pop_sizes[i] <- ms
}

# mean pop size
P30D3A16_pop_mean <- mean(P30D3A16_pop_sizes)

# standard deviation pop size
P30D3A16_pop_sd <- sd(P30D3A16_pop_sizes)

# standard error pop size
P30D3A16_pop_se <- sd(P30D3A16_pop_sizes)/sqrt(length(P30D3A16_pop_sizes))

#### Calculate Stochastic Growth Rate

P30D3A16_lambda_s <- stoch.growth.rate(P30D3A16_matrices, prob = NULL, maxt = time,
                                       verbose = TRUE)

#convert from log
P30D3A16_lambda_s$approx <- exp(P30D3A16_lambda_s$approx)
P30D3A16_lambda_s$sim <- exp(P30D3A16_lambda_s$sim)
P30D3A16_lambda_s$sim.CI <- exp(P30D3A16_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P30D3A16_quasi <- stoch.quasi.ext(P30D3A16_matrices, n0= P30D3A16_n0, Nx = 50, tmax = time, maxruns = 1,
                                  nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P30D3A16_quasi_df <- data.frame(P30D3A16_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P30D3A16_quasi_plot <- ggplot(P30D3A16_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P30D3A16_sens <- stoch.sens(P30D3A16_matrices, tlimit=time)

P30D3A16_elas <- P30D3A16_sens$elasticities

P30D3A16_elas_v <- c(P30D3A16_elas[1,1], P30D3A16_elas[1,2], P30D3A16_elas[1,3], P30D3A16_elas[2,1], P30D3A16_elas[2,2], P30D3A16_elas[2,3], P30D3A16_elas[3,1], P30D3A16_elas[3,2], P30D3A16_elas[3,3])

P30D3A16_elas_df <- data.frame(P30D3A16_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P30D3A16_elas_plot <- ggplot(P30D3A16_elas_df, aes(x = stage, y= P30D3A16_elas_v)) + 
  labs(x = "Vital rate", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 4 and Adult Duration of 15----

## Stage duration
duration <- c(1, 4, 15)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
P30D4A15_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P30D4A15_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P30D4A15_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

P30D4A15_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P30D4A15_matrices, n = P30D4A15_n0, time = time) 
  P30D4A15_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P30D4A15_total_pop <- lapply(P30D4A15_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P30D4A15_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P30D4A15_total_pop[[i]])
  P30D4A15_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P30D4A15_plot_data <- bind_rows(P30D4A15_df_plots, .id = "id")

# Plot projection
P30D4A15_plot <- ggplot(P30D4A15_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P30D4A15_mean_plot_data <- P30D4A15_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P30D4A15_plot_pred <- P30D4A15_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P30D4A15_mean_plot <- ggplot(P30D4A15_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P30D4A15_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P30D4A15_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P30D4A15_total_pop[[i]][time]
  P30D4A15_pop_sizes[i] <- ms
}

# mean pop size
P30D4A15_pop_mean <- mean(P30D4A15_pop_sizes)

# standard deviation pop size
P30D4A15_pop_sd <- sd(P30D4A15_pop_sizes)

# standard error pop size
P30D4A15_pop_se <- sd(P30D4A15_pop_sizes)/sqrt(length(P30D4A15_pop_sizes))

#### Calculate Stochastic Growth Rate

P30D4A15_lambda_s <- stoch.growth.rate(P30D4A15_matrices, prob = NULL, maxt = time,
                                       verbose = TRUE)

#convert from log
P30D4A15_lambda_s$approx <- exp(P30D4A15_lambda_s$approx)
P30D4A15_lambda_s$sim <- exp(P30D4A15_lambda_s$sim)
P30D4A15_lambda_s$sim.CI <- exp(P30D4A15_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P30D4A15_quasi <- stoch.quasi.ext(P30D4A15_matrices, n0= P30D4A15_n0, Nx = 50, tmax = time, maxruns = 1,
                                  nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P30D4A15_quasi_df <- data.frame(P30D4A15_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P30D4A15_quasi_plot <- ggplot(P30D4A15_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P30D4A15_sens <- stoch.sens(P30D4A15_matrices, tlimit=time)

P30D4A15_elas <- P30D4A15_sens$elasticities

P30D4A15_elas_v <- c(P30D4A15_elas[1,1], P30D4A15_elas[1,2], P30D4A15_elas[1,3], P30D4A15_elas[2,1], P30D4A15_elas[2,2], P30D4A15_elas[2,3], P30D4A15_elas[3,1], P30D4A15_elas[3,2], P30D4A15_elas[3,3])

P30D4A15_elas_df <- data.frame(P30D4A15_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P30D4A15_elas_plot <- ggplot(P30D4A15_elas_df, aes(x = stage, y= P30D4A15_elas_v)) + 
  labs(x = "Vital rate", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")



#### LSA for Immature Duration of 5 and Adult Duration of 14 ----

## Stage duration
duration <- c(1, 5, 14)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
P30D5A14_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P30D5A14_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P30D5A14_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

P30D5A14_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P30D5A14_matrices, n = P30D5A14_n0, time = time) 
  P30D5A14_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P30D5A14_total_pop <- lapply(P30D5A14_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P30D5A14_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P30D5A14_total_pop[[i]])
  P30D5A14_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P30D5A14_plot_data <- bind_rows(P30D5A14_df_plots, .id = "id")

# Plot projection
P30D5A14_plot <- ggplot(P30D5A14_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P30D5A14_mean_plot_data <- P30D5A14_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P30D5A14_plot_pred <- P30D5A14_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P30D5A14_mean_plot <- ggplot(P30D5A14_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P30D5A14_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P30D5A14_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P30D5A14_total_pop[[i]][time]
  P30D5A14_pop_sizes[i] <- ms
}

# mean pop size
P30D5A14_pop_mean <- mean(P30D5A14_pop_sizes)

# standard deviation pop size
P30D5A14_pop_sd <- sd(P30D5A14_pop_sizes)

# standard error pop size
P30D5A14_pop_se <- sd(P30D5A14_pop_sizes)/sqrt(length(P30D5A14_pop_sizes))

#### Calculate Stochastic Growth Rate

P30D5A14_lambda_s <- stoch.growth.rate(P30D5A14_matrices, prob = NULL, maxt = time,
                                       verbose = TRUE)

#convert from log
P30D5A14_lambda_s$approx <- exp(P30D5A14_lambda_s$approx)
P30D5A14_lambda_s$sim <- exp(P30D5A14_lambda_s$sim)
P30D5A14_lambda_s$sim.CI <- exp(P30D5A14_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P30D5A14_quasi <- stoch.quasi.ext(P30D5A14_matrices, n0= P30D5A14_n0, Nx = 50, tmax = time, maxruns = 1,
                                  nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P30D5A14_quasi_df <- data.frame(P30D5A14_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P30D5A14_quasi_plot <- ggplot(P30D5A14_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P30D5A14_sens <- stoch.sens(P30D5A14_matrices, tlimit=time)

P30D5A14_elas <- P30D5A14_sens$elasticities

P30D5A14_elas_v <- c(P30D5A14_elas[1,1], P30D5A14_elas[1,2], P30D5A14_elas[1,3], P30D5A14_elas[2,1], P30D5A14_elas[2,2], P30D5A14_elas[2,3], P30D5A14_elas[3,1], P30D5A14_elas[3,2], P30D5A14_elas[3,3])

P30D5A14_elas_df <- data.frame(P30D5A14_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P30D5A14_elas_plot <- ggplot(P30D5A14_elas_df, aes(x = stage, y= P30D5A14_elas_v)) + 
  labs(x = "Vital rate", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")



#### PLOTS ----

## Stochastic Population Projection Plot
P30_A20_plot <- P30D1A18_plot + P30D2A17_plot + P30D3A16_plot + P30D4A15_plot + P30D5A14_plot

## Mean and CI Stochastic Population Plot
P30_A20_mean_plot <- P30D1A18_mean_plot + P30D2A17_mean_plot + P30D3A16_mean_plot + P30D4A15_mean_plot + P30D5A14_mean_plot


## Stochastic Population Growth (Lambda s)

P30_lambda_approx <- c(P30D1A18_lambda_s$approx, P30D2A17_lambda_s$approx, P30D3A16_lambda_s$approx, P30D4A15_lambda_s$approx, P30D5A14_lambda_s$approx)

P30_lambda_sim <- c(P30D1A18_lambda_s$sim, P30D2A17_lambda_s$sim, P30D3A16_lambda_s$sim, P30D4A15_lambda_s$sim, P30D5A14_lambda_s$sim)

P30_lower_CI <- c(P30D1A18_lambda_s$sim.CI[1], P30D2A17_lambda_s$sim.CI[1], P30D3A16_lambda_s$sim.CI[1], P30D4A15_lambda_s$sim.CI[1], P30D5A14_lambda_s$sim.CI[1])

P30_upper_CI <- c(P30D1A18_lambda_s$sim.CI[2], P30D2A17_lambda_s$sim.CI[2], P30D3A16_lambda_s$sim.CI[2], P30D4A15_lambda_s$sim.CI[2], P30D5A14_lambda_s$sim.CI[2])

stage_duration <- c("1 year", "2 years", "3 years", "4 years", "5 years")

P30_lambda_df <- data.frame(stage_duration, P30_lambda_approx, P30_lambda_sim, P30_upper_CI, P30_lower_CI)

P30_lambda_plot <- ggplot(P30_lambda_df) +
  geom_point(aes(x = stage_duration, y = P30_lambda_sim), fill = "grey20", size = 2) +
  geom_errorbar(aes(x = stage_duration, ymin = P30_lower_CI, ymax = P30_upper_CI), width = 0.2) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  scale_x_discrete(labels=c("1 year" = "1", "2 years" = "2", "3 years" = "3", "4 years" = "4", "5 years" = "5")) +
  labs(x = "Immature stage duration (years)", y = "Lambda for stochastic population growth")


## Quasi-extinction Threshold Plots

P30A20_quasi_df <- rbind.data.frame(P30D1A18_quasi_df, P30D2A17_quasi_df, P30D3A16_quasi_df, P30D4A15_quasi_df, P30D5A14_quasi_df)

P30A20_quasi_plot <- ggplot(P30A20_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  scale_colour_discrete(name = "Immature stage \nduration", breaks = c("P30D1A18_quasi", "P30D2A17_quasi", "P30D3A16_quasi", "P30D4A15_quasi", "P30D5A14_quasi"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"))

#P30_A20_quasi_plots <- P30D1A18_quasi_plot + P30D2A17_quasi_plot + P30D3A16_quasi_plot + P30D4A15_quasi_plot + P30D5A14_quasi_plot


#Elasticity analysis plots

P30A20_elas_df<- rbind.data.frame(P30D1A18_elas_df, P30D2A17_elas_df, P30D3A16_elas_df, P30D4A15_elas_df, P30D5A14_elas_df)

P30A20_elas_plot <- ggplot(P30A20_elas_df, aes(x = stage, y= elasticity, fill = duration)) + 
  labs(x = "Vital rate", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(name = "Immature stage  \nduration", breaks = c("P30D1A18_elas_v", "P30D2A17_elas_v", "P30D3A16_elas_v", "P30D4A15_elas_v", "P30D5A14_elas_v"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"), values = c("grey65", "grey40", "grey35", "grey15", "grey0"))

#image2(P30D1A18_elas)
#image2(P30D2A17_elas)
#image2(P30D3A16_elas)
#image2(P30D4A15_elas)
#image2(P30D5A14_elas)

