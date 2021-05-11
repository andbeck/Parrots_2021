## YSA Stochastic Model of Population Growth
# With changes in immature stage duration of 1:5 years.
# With adult stage duration of 10:6 years.
# With imputed survival rates.
# Density independent.

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
survival <- c(0.71, 0.838, 0.838)

# Current population is estimated around 1000 individuals. 1:1 sex ratio means female population is 500
Nc <- 500

# Time to project to
time <- 100

#### YSA Simulated Vital Rates for LSA ----

set.seed(2021)

# Number of simulations
n_sim <- 1000

# Fledgling survival
s1 <- sapply(1:n_sim, function(x) betaval(0.71, 0.2))

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


#### LSA for Immature Duration of 1 and Adult Duration 0f 10----

## Stage duration
duration <- c(1, 1, 10)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
D1A10_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D1A10_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D1A10_matrices[[i]] <- mpm
}

head(D1A10_matrices)

## Repeat Stochastic Population Growth 

D1A10_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D1A10_matrices, n = D1A10_n0, time = time)
  D1A10_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D1A10_total_pop <- lapply(D1A10_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D1A10_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D1A10_total_pop[[i]])
  D1A10_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D1A10_plot_data <- bind_rows(D1A10_df_plots, .id = "id")

# Plot projection
D1A10_plot <- ggplot(D1A10_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D1A10_mean_plot_data <- D1A10_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D1A10_plot_pred <- D1A10_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D1A10_mean_plot <- ggplot(D1A10_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D1A10_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D1A10_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D1A10_total_pop[[i]][time]
  D1A10_pop_sizes[i] <- ms
}

#mean pop size
D1A10_pop_mean <- mean(D1A10_pop_sizes)

# standard deviation pop size
D1A10_pop_sd <- sd(D1A10_pop_sizes)

# standard error pop size
D1A10_pop_se <- sd(D1A10_pop_sizes)/sqrt(length(D1A10_pop_sizes))


#### Calculate Stochastic Growth Rate

D1A10_lambda_s <- stoch.growth.rate(D1A10_matrices, prob = NULL, maxt = time,
                                   verbose = TRUE)

#convert from log
D1A10_lambda_s$approx <- exp(D1A10_lambda_s$approx)
D1A10_lambda_s$sim <- exp(D1A10_lambda_s$sim)
D1A10_lambda_s$sim.CI <- exp(D1A10_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D1A10_quasi <- stoch.quasi.ext(D1A10_matrices, n0= D1A10_n0, Nx = 50, tmax = time, maxruns = 1,
                              nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D1A10_quasi_df <- data.frame(D1A10_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D1A10_quasi_plot <- ggplot(D1A10_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D1A10_sens <- stoch.sens(D1A10_matrices, tlimit=time)

D1A10_elas <- D1A10_sens$elasticities

D1A10_elas_v <- c(D1A10_elas[1,1], D1A10_elas[1,2], D1A10_elas[1,3], D1A10_elas[2,1], D1A10_elas[2,2], D1A10_elas[2,3], D1A10_elas[3,1], D1A10_elas[3,2], D1A10_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D1A10_elas_df <- data.frame(D1A10_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D1A10_elas_plot <- ggplot(D1A10_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 2 and Adult Duration 0f 9----

## Stage duration
duration <- c(1, 2, 9)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
D2A9_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D2A9_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D2A9_matrices[[i]] <- mpm
}

head(D2A9_matrices)

## Repeat Stochastic Population Growth 

D2A9_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D2A9_matrices, n = D2A9_n0, time = time)
  D2A9_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D2A9_total_pop <- lapply(D2A9_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D2A9_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D2A9_total_pop[[i]])
  D2A9_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D2A9_plot_data <- bind_rows(D2A9_df_plots, .id = "id")

# Plot projection
D2A9_plot <- ggplot(D2A9_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D2A9_mean_plot_data <- D2A9_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D2A9_plot_pred <- D2A9_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D2A9_mean_plot <- ggplot(D2A9_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D2A9_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D2A9_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D2A9_total_pop[[i]][time]
  D2A9_pop_sizes[i] <- ms
}

#mean pop size
D2A9_pop_mean <- mean(D2A9_pop_sizes)

# standard deviation pop size
D2A9_pop_sd <- sd(D2A9_pop_sizes)

# standard error pop size
D2A9_pop_se <- sd(D2A9_pop_sizes)/sqrt(length(D2A9_pop_sizes))


#### Calculate Stochastic Growth Rate

D2A9_lambda_s <- stoch.growth.rate(D2A9_matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)

#convert from log
D2A9_lambda_s$approx <- exp(D2A9_lambda_s$approx)
D2A9_lambda_s$sim <- exp(D2A9_lambda_s$sim)
D2A9_lambda_s$sim.CI <- exp(D2A9_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D2A9_quasi <- stoch.quasi.ext(D2A9_matrices, n0= D2A9_n0, Nx = 50, tmax = time, maxruns = 1,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D2A9_quasi_df <- data.frame(D2A9_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D2A9_quasi_plot <- ggplot(D2A9_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D2A9_sens <- stoch.sens(D2A9_matrices, tlimit=time)

D2A9_elas <- D2A9_sens$elasticities

D2A9_elas_v <- c(D2A9_elas[1,1], D2A9_elas[1,2], D2A9_elas[1,3], D2A9_elas[2,1], D2A9_elas[2,2], D2A9_elas[2,3], D2A9_elas[3,1], D2A9_elas[3,2], D2A9_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D2A9_elas_df <- data.frame(D2A9_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D2A9_elas_plot <- ggplot(D2A9_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 3 and Adult Duration of 8 ----

## Stage duration
duration <- c(1, 3, 8)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

## Initial population vector estimated from stable stage distribution
D3A8_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D3A8_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D3A8_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D3A8_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D3A8_matrices, n = D3A8_n0, time = time) 
  D3A8_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D3A8_total_pop <- lapply(D3A8_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D3A8_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D3A8_total_pop[[i]])
  D3A8_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D3A8_plot_data <- bind_rows(D3A8_df_plots, .id = "id")

# Plot projection
D3A8_plot <- ggplot(D3A8_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D3A8_mean_plot_data <- D3A8_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D3A8_plot_pred <- D3A8_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D3A8_mean_plot <- ggplot(D3A8_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D3A8_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D3A8_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D3A8_total_pop[[i]][time]
  D3A8_pop_sizes[i] <- ms
}

# mean pop size
D3A8_pop_mean <- mean(D3A8_pop_sizes)

# standard deviation pop size
D3A8_pop_sd <- sd(D3A8_pop_sizes)

#standard error pop size
D3A8_pop_se <- sd(D3A8_pop_sizes)/sqrt(length(D3A8_pop_sizes))


#### Calculate Stochastic Growth Rate

D3A8_lambda_s <- stoch.growth.rate(D3A8_matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)

#convert from log
D3A8_lambda_s$approx <- exp(D3A8_lambda_s$approx)
D3A8_lambda_s$sim <- exp(D3A8_lambda_s$sim)
D3A8_lambda_s$sim.CI <- exp(D3A8_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D3A8_quasi <- stoch.quasi.ext(D3A8_matrices, n0= D3A8_n0, Nx = 50, tmax = time, maxruns = 1,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D3A8_quasi_df <- data.frame(D3A8_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D3A8_quasi_plot <- ggplot(D3A8_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D3A8_sens <- stoch.sens(D3A8_matrices, tlimit=time)

D3A8_elas <- D3A8_sens$elasticities

D3A8_elas_v <- c(D3A8_elas[1,1], D3A8_elas[1,2], D3A8_elas[1,3], D3A8_elas[2,1], D3A8_elas[2,2], D3A8_elas[2,3], D3A8_elas[3,1], D3A8_elas[3,2], D3A8_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D3A8_elas_df <- data.frame(D3A8_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D3A8_elas_plot <- ggplot(D3A8_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 4 and Adult Duration of 7 ----

## Stage duration
duration <- c(1, 4, 7)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D4A7_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D4A7_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D4A7_matrices[[i]] <- mpm
}


##Repeat Stochastic Population Growth 

D4A7_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D4A7_matrices, n = D4A7_n0, time = time) 
  D4A7_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D4A7_total_pop <- lapply(D4A7_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D4A7_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D4A7_total_pop[[i]])
  D4A7_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D4A7_plot_data <- bind_rows(D4A7_df_plots, .id = "id")

# Plot projection
D4A7_plot <- ggplot(D4A7_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D4A7_mean_plot_data <- D4A7_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D4A7_plot_pred <- D4A7_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D4A7_mean_plot <- ggplot(D4A7_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D4A7_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D4A7_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D4A7_total_pop[[i]][time]
  D4A7_pop_sizes[i] <- ms
}

# mean pop size
D4A7_pop_mean <- mean(D4A7_pop_sizes)

# standard deviation pop size
D4A7_pop_sd <- sd(D4A7_pop_sizes)

# standard error pop size
D4A7_pop_se <- sd(D4A7_pop_sizes)/sqrt(length(D4A7_pop_sizes))

#### Calculate Stochastic Growth Rate

D4A7_lambda_s <- stoch.growth.rate(D4A7_matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)

#convert from log
D4A7_lambda_s$approx <- exp(D4A7_lambda_s$approx)
D4A7_lambda_s$sim <- exp(D4A7_lambda_s$sim)
D4A7_lambda_s$sim.CI <- exp(D4A7_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D4A7_quasi <- stoch.quasi.ext(D4A7_matrices, n0= D4A7_n0, Nx = 50, tmax = time, maxruns = 1,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D4A7_quasi_df <- data.frame(D4A7_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D4A7_quasi_plot <- ggplot(D4A7_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D4A7_sens <- stoch.sens(D4A7_matrices, tlimit=time)

D4A7_elas <- D4A7_sens$elasticities

D4A7_elas_v <- c(D4A7_elas[1,1], D4A7_elas[1,2], D4A7_elas[1,3], D4A7_elas[2,1], D4A7_elas[2,2], D4A7_elas[2,3], D4A7_elas[3,1], D4A7_elas[3,2], D4A7_elas[3,3])

D4A7_elas_df <- data.frame(D4A7_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D4A7_elas_plot <- ggplot(D4A7_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 5 and Adult Duration of 6 ----

## Stage duration
duration <- c(1, 5, 6)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D5A6_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D5A6_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D5A6_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D5A6_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D5A6_matrices, n = D5A6_n0, time = time) 
  D5A6_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D5A6_total_pop <- lapply(D5A6_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D5A6_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D5A6_total_pop[[i]])
  D5A6_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D5A6_plot_data <- bind_rows(D5A6_df_plots, .id = "id")

# Plot projection
D5A6_plot <- ggplot(D5A6_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D5A6_mean_plot_data <- D5A6_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D5A6_plot_pred <- D5A6_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D5A6_mean_plot <- ggplot(D5A6_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D5A6_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D5A6_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D5A6_total_pop[[i]][time]
  D5A6_pop_sizes[i] <- ms
}

# mean pop size
D5A6_pop_mean <- mean(D5A6_pop_sizes)

# standard deviation pop size
D5A6_pop_sd <- sd(D5A6_pop_sizes)

# standard error pop size
D5A6_pop_se <- sd(D5A6_pop_sizes)/sqrt(length(D5A6_pop_sizes))

#### Calculate Stochastic Growth Rate

D5A6_lambda_s <- stoch.growth.rate(D5A6_matrices, prob = NULL, maxt = time,
                                   verbose = TRUE)

#convert from log
D5A6_lambda_s$approx <- exp(D5A6_lambda_s$approx)
D5A6_lambda_s$sim <- exp(D5A6_lambda_s$sim)
D5A6_lambda_s$sim.CI <- exp(D5A6_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D5A6_quasi <- stoch.quasi.ext(D5A6_matrices, n0= D5A6_n0, Nx = 50, tmax = time, maxruns = 1,
                              nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D5A6_quasi_df <- data.frame(D5A6_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D5A6_quasi_plot <- ggplot(D5A6_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D5A6_sens <- stoch.sens(D5A6_matrices, tlimit=time)

D5A6_elas <- D5A6_sens$elasticities

D5A6_elas_v <- c(D5A6_elas[1,1], D5A6_elas[1,2], D5A6_elas[1,3], D5A6_elas[2,1], D5A6_elas[2,2], D5A6_elas[2,3], D5A6_elas[3,1], D5A6_elas[3,2], D5A6_elas[3,3])

D5A6_elas_df <- data.frame(D5A6_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D5A6_elas_plot <- ggplot(D5A6_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### PLOTS ----

## Stochastic Population Projection Plot
A10_plot <- D1A10_plot + D2A9_plot + D3A8_plot + D4A7_plot + D5A6_plot

#A10_plot_same_lims <- D1A10_plot + ylim(0, 1.0e+07) + D2A9_plot + ylim(0, 1.0e+07) + D3A8_plot + ylim(0, 1.0e+07) + D4A7_plot + ylim(0, 1.0e+07) + D5A6_plot + ylim(0, 1.0e+07)

## Mean and CI Stochastic Population Plot
A10_mean_plot <- D1A10_mean_plot + D2A9_mean_plot + D3A8_mean_plot + D4A7_mean_plot + D5A6_mean_plot


# Stochastic Population Growth (Lambda s)

A10_lambda_approx <- c(D1A10_lambda_s$approx, D2A9_lambda_s$approx, D3A8_lambda_s$approx, D4A7_lambda_s$approx, D5A6_lambda_s$approx)

A10_lambda_sim <- c(D1A10_lambda_s$sim, D2A9_lambda_s$sim, D3A8_lambda_s$sim, D4A7_lambda_s$sim, D5A6_lambda_s$sim)

A10_lower_CI <- c(D1A10_lambda_s$sim.CI[1], D2A9_lambda_s$sim.CI[1], D3A8_lambda_s$sim.CI[1], D4A7_lambda_s$sim.CI[1], D5A6_lambda_s$sim.CI[1])

A10_upper_CI <- c(D1A10_lambda_s$sim.CI[2], D2A9_lambda_s$sim.CI[2], D3A8_lambda_s$sim.CI[2], D4A7_lambda_s$sim.CI[2], D5A6_lambda_s$sim.CI[2])

stage_duration <- c("1 year", "2 years", "3 years", "4 years", "5 years")

A10_lambda_df <- data.frame(stage_duration, A10_lambda_approx, A10_lambda_sim, A10_upper_CI, A10_lower_CI)

A10_lambda_plot <- ggplot(A10_lambda_df) +
  geom_point(aes(x = stage_duration, y = A10_lambda_sim), fill = "grey20", size = 2) +
  geom_errorbar(aes(x = stage_duration, ymin = A10_lower_CI, ymax = A10_upper_CI), width = 0.2) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  scale_x_discrete(labels=c("1 year" = "1", "2 years" = "2", "3 years" = "3", "4 years" = "4", "5 years" = "5")) +
  labs(x = "Immature stage duration (years)", y = "Lambda for stochastic population growth")


## Quasi-extinction Threshold Plots

A10_quasi_df<- rbind.data.frame(D1A10_quasi_df, D2A9_quasi_df, D3A8_quasi_df, D4A7_quasi_df, D5A6_quasi_df)

A10_quasi_plot <- ggplot(A10_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  scale_colour_discrete(name = "Immature stage \nduration", breaks = c("D1A10_quasi", "D2A9_quasi", "D3A8_quasi", "D4A7_quasi", "D5A6_quasi"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"))

#A10_quasi_plots <- D1A10_quasi_plot + D2A9_quasi_plot + D3A8_quasi_plot + D4A7_quasi_plot + D5A6_quasi_plot


#Elasticity analysis plots

A10_elas_df<- rbind.data.frame(D1A10_elas_df, D2A9_elas_df, D3A8_elas_df, D4A7_elas_df, D5A6_elas_df)

A10_elas_plot <- ggplot(A10_elas_df, aes(x = stage, y= elasticity, fill = duration)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(name = "Immature stage  \nduration", breaks = c("D1A10_elas_v", "D2A9_elas_v", "D3A8_elas_v", "D4A7_elas_v", "D5A6_elas_v"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"), values = c("grey65", "grey40", "grey35", "grey15", "grey0"))


#image2(D1A10_elas)
#image2(D2A9_elas)
#image2(D3A8_elas)
#image2(D4A7_elas)
#image2(D5A6_elas)

