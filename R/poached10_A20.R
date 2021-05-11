# YSA Stochastic Model of Population Growth under Poaching Pressure of 10% reduced fledgling survival
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
fecundity <- c(0, 0, 1.6*total_summary$mean_hatch[1]*total_summary$mean_nestling_surv[1])

# Mean survival (0.71 is from Salinas-Melgoza & Renton 2007, 0.838 is survival from imputation)
survival <- c(0.71*0.9, 0.838, 0.838)

# Current population is estimated around 1000 individuals. 1:1 sex ratio means female population is 500
Nc <- 500

# Time to project to
time <- 100

#### YSA Simulated Vital Rates for LSA ----

set.seed(2021)

# Number of simulations
n_sim <- 1000

# Fledgling survival
s1 <- sapply(1:n_sim, function(x) betaval(0.71*0.9, 0.2))

# Immature survival
s2 <- sapply(1:n_sim, function(x) betaval(0.838, 0.051))

# Adult survival
s3 <- sapply(1:n_sim, function(x) betaval(0.838, 0.051))

# Fecundity
m3 <- rlnorm(n = n_sim, 
             log(1.6*total_summary$mean_hatch[1]*total_summary$mean_nestling_surv[1]), 
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



#### LSA for Immature Duration of 1 and Adult Duration 0f 20----

## Stage duration
duration <- c(1, 1, 20)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
P10D1A20_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P10D1A20_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P10D1A20_matrices[[i]] <- mpm
}

head(P10D1A20_matrices)

## Repeat Stochastic Population Growth 

P10D1A20_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P10D1A20_matrices, n = P10D1A20_n0, time = time)
  P10D1A20_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P10D1A20_total_pop <- lapply(P10D1A20_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P10D1A20_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P10D1A20_total_pop[[i]])
  P10D1A20_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P10D1A20_plot_data <- bind_rows(P10D1A20_df_plots, .id = "id")

# Plot projection
P10D1A20_plot <- ggplot(P10D1A20_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P10D1A20_mean_plot_data <- P10D1A20_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P10D1A20_plot_pred <- P10D1A20_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P10D1A20_mean_plot <- ggplot(P10D1A20_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P10D1A20_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P10D1A20_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P10D1A20_total_pop[[i]][time]
  P10D1A20_pop_sizes[i] <- ms
}

#mean pop size
P10D1A20_pop_mean <- mean(P10D1A20_pop_sizes)

# standard deviation pop size
P10D1A20_pop_sd <- sd(P10D1A20_pop_sizes)

# standard error pop size
P10D1A20_pop_se <- sd(P10D1A20_pop_sizes)/sqrt(length(P10D1A20_pop_sizes))


#### Calculate Stochastic Growth Rate

P10D1A20_lambda_s <- stoch.growth.rate(P10D1A20_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
P10D1A20_lambda_s$approx <- exp(P10D1A20_lambda_s$approx)
P10D1A20_lambda_s$sim <- exp(P10D1A20_lambda_s$sim)
P10D1A20_lambda_s$sim.CI <- exp(P10D1A20_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P10D1A20_quasi <- stoch.quasi.ext(P10D1A20_matrices, n0= P10D1A20_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P10D1A20_quasi_df <- data.frame(P10D1A20_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P10D1A20_quasi_plot <- ggplot(P10D1A20_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P10D1A20_sens <- stoch.sens(P10D1A20_matrices, tlimit=time)

P10D1A20_elas <- P10D1A20_sens$elasticities

P10D1A20_elas_v <- c(P10D1A20_elas[1,1], P10D1A20_elas[1,2], P10D1A20_elas[1,3], P10D1A20_elas[2,1], P10D1A20_elas[2,2], P10D1A20_elas[2,3], P10D1A20_elas[3,1], P10D1A20_elas[3,2], P10D1A20_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

P10D1A20_elas_df <- data.frame(P10D1A20_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P10D1A20_elas_plot <- ggplot(P10D1A20_elas_df, aes(x = stage, y= P10D1A20_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 2 and Adult Duration of 19 ----

## Stage duration
duration <- c(1, 2, 19)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

## Initial population vector estimated from stable stage distribution
P10D2A19_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P10D2A19_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P10D2A19_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

P10D2A19_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P10D2A19_matrices, n = P10D2A19_n0, time = time) 
  P10D2A19_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P10D2A19_total_pop <- lapply(P10D2A19_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P10D2A19_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P10D2A19_total_pop[[i]])
  P10D2A19_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P10D2A19_plot_data <- bind_rows(P10D2A19_df_plots, .id = "id")

# Plot projection
P10D2A19_plot <- ggplot(P10D2A19_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P10D2A19_mean_plot_data <- P10D2A19_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P10D2A19_plot_pred <- P10D2A19_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P10D2A19_mean_plot <- ggplot(P10D2A19_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P10D2A19_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P10D2A19_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P10D2A19_total_pop[[i]][time]
  P10D2A19_pop_sizes[i] <- ms
}

# mean pop size
P10D2A19_pop_mean <- mean(P10D2A19_pop_sizes)

# standard deviation pop size
P10D2A19_pop_sd <- sd(P10D2A19_pop_sizes)

#standard error pop size
P10D2A19_pop_se <- sd(P10D2A19_pop_sizes)/sqrt(length(P10D2A19_pop_sizes))


#### Calculate Stochastic Growth Rate

P10D2A19_lambda_s <- stoch.growth.rate(P10D2A19_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
P10D2A19_lambda_s$approx <- exp(P10D2A19_lambda_s$approx)
P10D2A19_lambda_s$sim <- exp(P10D2A19_lambda_s$sim)
P10D2A19_lambda_s$sim.CI <- exp(P10D2A19_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P10D2A19_quasi <- stoch.quasi.ext(P10D2A19_matrices, n0= P10D2A19_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P10D2A19_quasi_df <- data.frame(P10D2A19_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P10D2A19_quasi_plot <- ggplot(P10D2A19_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P10D2A19_sens <- stoch.sens(P10D2A19_matrices, tlimit=time)

P10D2A19_elas <- P10D2A19_sens$elasticities

P10D2A19_elas_v <- c(P10D2A19_elas[1,1], P10D2A19_elas[1,2], P10D2A19_elas[1,3], P10D2A19_elas[2,1], P10D2A19_elas[2,2], P10D2A19_elas[2,3], P10D2A19_elas[3,1], P10D2A19_elas[3,2], P10D2A19_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

P10D2A19_elas_df <- data.frame(P10D2A19_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P10D2A19_elas_plot <- ggplot(P10D2A19_elas_df, aes(x = stage, y= P10D2A19_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 3 and Adult Duration of 18 ----

## Stage duration
duration <- c(1, 3, 18)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
P10D3A18_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P10D3A18_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P10D3A18_matrices[[i]] <- mpm
}


##Repeat Stochastic Population Growth 

P10D3A18_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P10D3A18_matrices, n = P10D3A18_n0, time = time) 
  P10D3A18_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P10D3A18_total_pop <- lapply(P10D3A18_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P10D3A18_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P10D3A18_total_pop[[i]])
  P10D3A18_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P10D3A18_plot_data <- bind_rows(P10D3A18_df_plots, .id = "id")

# Plot projection
P10D3A18_plot <- ggplot(P10D3A18_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P10D3A18_mean_plot_data <- P10D3A18_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P10D3A18_plot_pred <- P10D3A18_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P10D3A18_mean_plot <- ggplot(P10D3A18_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P10D3A18_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P10D3A18_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P10D3A18_total_pop[[i]][time]
  P10D3A18_pop_sizes[i] <- ms
}

# mean pop size
P10D3A18_pop_mean <- mean(P10D3A18_pop_sizes)

# standard deviation pop size
P10D3A18_pop_sd <- sd(P10D3A18_pop_sizes)

# standard error pop size
P10D3A18_pop_se <- sd(P10D3A18_pop_sizes)/sqrt(length(P10D3A18_pop_sizes))

#### Calculate Stochastic Growth Rate

P10D3A18_lambda_s <- stoch.growth.rate(P10D3A18_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
P10D3A18_lambda_s$approx <- exp(P10D3A18_lambda_s$approx)
P10D3A18_lambda_s$sim <- exp(P10D3A18_lambda_s$sim)
P10D3A18_lambda_s$sim.CI <- exp(P10D3A18_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P10D3A18_quasi <- stoch.quasi.ext(P10D3A18_matrices, n0= P10D3A18_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P10D3A18_quasi_df <- data.frame(P10D3A18_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P10D3A18_quasi_plot <- ggplot(P10D3A18_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P10D3A18_sens <- stoch.sens(P10D3A18_matrices, tlimit=time)

P10D3A18_elas <- P10D3A18_sens$elasticities

P10D3A18_elas_v <- c(P10D3A18_elas[1,1], P10D3A18_elas[1,2], P10D3A18_elas[1,3], P10D3A18_elas[2,1], P10D3A18_elas[2,2], P10D3A18_elas[2,3], P10D3A18_elas[3,1], P10D3A18_elas[3,2], P10D3A18_elas[3,3])

P10D3A18_elas_df <- data.frame(P10D3A18_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P10D3A18_elas_plot <- ggplot(P10D3A18_elas_df, aes(x = stage, y= P10D3A18_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 4 and Adult Duration of 17 ----

## Stage duration
duration <- c(1, 4, 17)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
P10D4A17_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P10D4A17_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P10D4A17_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

P10D4A17_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P10D4A17_matrices, n = P10D4A17_n0, time = time) 
  P10D4A17_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P10D4A17_total_pop <- lapply(P10D4A17_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P10D4A17_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P10D4A17_total_pop[[i]])
  P10D4A17_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P10D4A17_plot_data <- bind_rows(P10D4A17_df_plots, .id = "id")

# Plot projection
P10D4A17_plot <- ggplot(P10D4A17_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P10D4A17_mean_plot_data <- P10D4A17_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P10D4A17_plot_pred <- P10D4A17_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P10D4A17_mean_plot <- ggplot(P10D4A17_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P10D4A17_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P10D4A17_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P10D4A17_total_pop[[i]][time]
  P10D4A17_pop_sizes[i] <- ms
}

# mean pop size
P10D4A17_pop_mean <- mean(P10D4A17_pop_sizes)

# standard deviation pop size
P10D4A17_pop_sd <- sd(P10D4A17_pop_sizes)

# standard error pop size
P10D4A17_pop_se <- sd(P10D4A17_pop_sizes)/sqrt(length(P10D4A17_pop_sizes))

#### Calculate Stochastic Growth Rate

P10D4A17_lambda_s <- stoch.growth.rate(P10D4A17_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
P10D4A17_lambda_s$approx <- exp(P10D4A17_lambda_s$approx)
P10D4A17_lambda_s$sim <- exp(P10D4A17_lambda_s$sim)
P10D4A17_lambda_s$sim.CI <- exp(P10D4A17_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P10D4A17_quasi <- stoch.quasi.ext(P10D4A17_matrices, n0= P10D4A17_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P10D4A17_quasi_df <- data.frame(P10D4A17_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P10D4A17_quasi_plot <- ggplot(P10D4A17_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P10D4A17_sens <- stoch.sens(P10D4A17_matrices, tlimit=time)

P10D4A17_elas <- P10D4A17_sens$elasticities

P10D4A17_elas_v <- c(P10D4A17_elas[1,1], P10D4A17_elas[1,2], P10D4A17_elas[1,3], P10D4A17_elas[2,1], P10D4A17_elas[2,2], P10D4A17_elas[2,3], P10D4A17_elas[3,1], P10D4A17_elas[3,2], P10D4A17_elas[3,3])

P10D4A17_elas_df <- data.frame(P10D4A17_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P10D4A17_elas_plot <- ggplot(P10D4A17_elas_df, aes(x = stage, y= P10D4A17_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")



#### LSA for Immature Duration of 5 and Adult Duration of 16 ----

## Stage duration
duration <- c(1, 5, 16)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
P10D5A16_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

P10D5A16_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  P10D5A16_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

P10D5A16_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(P10D5A16_matrices, n = P10D5A16_n0, time = time) 
  P10D5A16_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
P10D5A16_total_pop <- lapply(P10D5A16_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

P10D5A16_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = P10D5A16_total_pop[[i]])
  P10D5A16_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
P10D5A16_plot_data <- bind_rows(P10D5A16_df_plots, .id = "id")

# Plot projection
P10D5A16_plot <- ggplot(P10D5A16_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

P10D5A16_mean_plot_data <- P10D5A16_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
P10D5A16_plot_pred <- P10D5A16_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
P10D5A16_mean_plot <- ggplot(P10D5A16_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = P10D5A16_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

P10D5A16_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- P10D5A16_total_pop[[i]][time]
  P10D5A16_pop_sizes[i] <- ms
}

# mean pop size
P10D5A16_pop_mean <- mean(P10D5A16_pop_sizes)

# standard deviation pop size
P10D5A16_pop_sd <- sd(P10D5A16_pop_sizes)

# standard error pop size
P10D5A16_pop_se <- sd(P10D5A16_pop_sizes)/sqrt(length(P10D5A16_pop_sizes))

#### Calculate Stochastic Growth Rate

P10D5A16_lambda_s <- stoch.growth.rate(P10D5A16_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
P10D5A16_lambda_s$approx <- exp(P10D5A16_lambda_s$approx)
P10D5A16_lambda_s$sim <- exp(P10D5A16_lambda_s$sim)
P10D5A16_lambda_s$sim.CI <- exp(P10D5A16_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

P10D5A16_quasi <- stoch.quasi.ext(P10D5A16_matrices, n0= P10D5A16_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
P10D5A16_quasi_df <- data.frame(P10D5A16_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

P10D5A16_quasi_plot <- ggplot(P10D5A16_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

P10D5A16_sens <- stoch.sens(P10D5A16_matrices, tlimit=time)

P10D5A16_elas <- P10D5A16_sens$elasticities

P10D5A16_elas_v <- c(P10D5A16_elas[1,1], P10D5A16_elas[1,2], P10D5A16_elas[1,3], P10D5A16_elas[2,1], P10D5A16_elas[2,2], P10D5A16_elas[2,3], P10D5A16_elas[3,1], P10D5A16_elas[3,2], P10D5A16_elas[3,3])

P10D5A16_elas_df <- data.frame(P10D5A16_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

P10D5A16_elas_plot <- ggplot(P10D5A16_elas_df, aes(x = stage, y= P10D5A16_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")



#### PLOTS ----

## Stochastic Population Projection Plot
P10_A20_plot <- P10D1A20_plot + P10D2A19_plot + P10D3A18_plot + P10D4A17_plot + P10D5A16_plot

## Mean and CI Stochastic Population Plot
P10_A20_mean_plot <- P10D1A20_mean_plot + P10D2A19_mean_plot + P10D3A18_mean_plot + P10D4A17_mean_plot + P10D5A16_mean_plot


## Stochastic Population Growth (Lambda s)

P10_lambda_approx <- c(P10D1A20_lambda_s$approx, P10D2A19_lambda_s$approx, P10D3A18_lambda_s$approx, P10D4A17_lambda_s$approx, P10D5A16_lambda_s$approx)

P10_lambda_sim <- c(P10D1A20_lambda_s$sim, P10D2A19_lambda_s$sim, P10D3A18_lambda_s$sim, P10D4A17_lambda_s$sim, P10D5A16_lambda_s$sim)

P10_lower_CI <- c(P10D1A20_lambda_s$sim.CI[1], P10D2A19_lambda_s$sim.CI[1], P10D3A18_lambda_s$sim.CI[1], P10D4A17_lambda_s$sim.CI[1], P10D5A16_lambda_s$sim.CI[1])

P10_upper_CI <- c(P10D1A20_lambda_s$sim.CI[2], P10D2A19_lambda_s$sim.CI[2], P10D3A18_lambda_s$sim.CI[2], P10D4A17_lambda_s$sim.CI[2], P10D5A16_lambda_s$sim.CI[2])

stage_duration <- c("1 year", "2 years", "3 years", "4 years", "5 years")

P10_lambda_df <- data.frame(stage_duration, P10_lambda_approx, P10_lambda_sim, P10_upper_CI, P10_lower_CI)

P10_lambda_plot <- ggplot(P10_lambda_df) +
  geom_point(aes(x = stage_duration, y = P10_lambda_sim), fill = "grey20", size = 2) +
  geom_errorbar(aes(x = stage_duration, ymin = P10_lower_CI, ymax = P10_upper_CI), width = 0.2) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  scale_x_discrete(labels=c("1 year" = "1", "2 years" = "2", "3 years" = "3", "4 years" = "4", "5 years" = "5")) +
  labs(x = "Immature stage duration (years)", y = "Lambda for stochastic population growth")


## Quasi-extinction Threshold Plots

P10A20_quasi_df <- rbind.data.frame(P10D1A20_quasi_df, P10D2A19_quasi_df, P10D3A18_quasi_df, P10D4A17_quasi_df, P10D5A16_quasi_df)

P10A20_quasi_plot <- ggplot(P10A20_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  scale_colour_discrete(name = "Immature stage \nduration", breaks = c("P10D1A20_quasi", "P10D2A19_quasi", "P10D3A18_quasi", "P10D4A17_quasi", "P10D5A16_quasi"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"))

#P10_A20_quasi_plots <- P10D1A20_quasi_plot + P10D2A19_quasi_plot + P10D3A18_quasi_plot + P10D4A17_quasi_plot + P10D5A16_quasi_plot


#Elasticity analysis plots

P10A20_elas_df<- rbind.data.frame(P10D1A20_elas_df, P10D2A19_elas_df, P10D3A18_elas_df, P10D4A17_elas_df, P10D5A16_elas_df)

P10A20_elas_plot <- ggplot(P10A20_elas_df, aes(x = stage, y= elasticity, fill = duration)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(name = "Immature stage  \nduration", breaks = c("P10D1A20_elas_v", "P10D2A19_elas_v", "P10D3A18_elas_v", "P10D4A17_elas_v", "P10D5A16_elas_v"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"), values = c("grey65", "grey40", "grey35", "grey15", "grey0"))

#image2(P10D1A20_elas)
#image2(P10D2A19_elas)
#image2(P10D3A18_elas)
#image2(P10D4A17_elas)
#image2(P10D5A16_elas)

