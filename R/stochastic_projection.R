## YSA Stochastic Model of Population Growth
# With changes in immature stage duration of 1:5 years.
# Constant adult stage duration of 10 years.
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

# Mean survival (0.71 is from Salinas-Melgoza & Renton 2007, 0.838 is adult survival from imputation)
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


#### LSA for Immature Duration of 1 ----

## Stage duration
duration <- c(1, 1, 10)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
D1_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D1_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D1_matrices[[i]] <- mpm
}

head(D1_matrices)

## Repeat Stochastic Population Growth 

D1_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D1_matrices, n = D1_n0, time = time)
  D1_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D1_total_pop <- lapply(D1_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D1_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D1_total_pop[[i]])
  D1_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D1_plot_data <- bind_rows(D1_df_plots, .id = "id")

# Plot projection
D1_plot <- ggplot(D1_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D1_mean_plot_data <- D1_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D1_plot_pred <- D1_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D1_mean_plot <- ggplot(D1_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D1_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D1_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D1_total_pop[[i]][time]
  D1_pop_sizes[i] <- ms
}

#mean pop size
D1_pop_mean <- mean(D1_pop_sizes)

# standard deviation pop size
D1_pop_sd <- sd(D1_pop_sizes)

# standard error pop size
D1_pop_se <- sd(D1_pop_sizes)/sqrt(length(D1_pop_sizes))


#### Calculate Stochastic Growth Rate

D1_lambda_s <- stoch.growth.rate(D1_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
D1_lambda_s$approx <- exp(D1_lambda_s$approx)
D1_lambda_s$sim <- exp(D1_lambda_s$sim)
D1_lambda_s$sim.CI <- exp(D1_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D1_quasi <- stoch.quasi.ext(D1_matrices, n0= D1_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D1_quasi_df <- data.frame(D1_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D1_quasi_plot <- ggplot(D1_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D1_sens <- stoch.sens(D1_matrices, tlimit=time)

D1_elas <- D1_sens$elasticities

D1_elas_v <- c(D1_elas[1,1], D1_elas[1,2], D1_elas[1,3], D1_elas[2,1], D1_elas[2,2], D1_elas[2,3], D1_elas[3,1], D1_elas[3,2], D1_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D1_elas_df <- data.frame(D1_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D1_elas_plot <- ggplot(D1_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 2 ----

## Stage duration
duration <- c(1, 2, 10)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
D2_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D2_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D2_matrices[[i]] <- mpm
}

head(D2_matrices)

## Repeat Stochastic Population Growth 

D2_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D2_matrices, n = D2_n0, time = time)
  D2_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D2_total_pop <- lapply(D2_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D2_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D2_total_pop[[i]])
  D2_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D2_plot_data <- bind_rows(D2_df_plots, .id = "id")

# Plot projection
D2_plot <- ggplot(D2_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D2_mean_plot_data <- D2_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D2_plot_pred <- D2_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D2_mean_plot <- ggplot(D2_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D2_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D2_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D2_total_pop[[i]][time]
  D2_pop_sizes[i] <- ms
}

#mean pop size
D2_pop_mean <- mean(D2_pop_sizes)

# standard deviation pop size
D2_pop_sd <- sd(D2_pop_sizes)

# standard error pop size
D2_pop_se <- sd(D2_pop_sizes)/sqrt(length(D2_pop_sizes))


#### Calculate Stochastic Growth Rate

D2_lambda_s <- stoch.growth.rate(D2_matrices, prob = NULL, maxt = time,
                                   verbose = TRUE)

#convert from log
D2_lambda_s$approx <- exp(D2_lambda_s$approx)
D2_lambda_s$sim <- exp(D2_lambda_s$sim)
D2_lambda_s$sim.CI <- exp(D2_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D2_quasi <- stoch.quasi.ext(D2_matrices, n0= D2_n0, Nx = 50, tmax = time, maxruns = 1,
                              nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D2_quasi_df <- data.frame(D2_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D2_quasi_plot <- ggplot(D2_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D2_sens <- stoch.sens(D2_matrices, tlimit=time)

D2_elas <- D2_sens$elasticities

D2_elas_v <- c(D2_elas[1,1], D2_elas[1,2], D2_elas[1,3], D2_elas[2,1], D2_elas[2,2], D2_elas[2,3], D2_elas[3,1], D2_elas[3,2], D2_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D2_elas_df <- data.frame(D2_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D2_elas_plot <- ggplot(D2_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 3 ----

## Stage duration
duration <- c(1, 3, 10)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

## Initial population vector estimated from stable stage distribution
D3_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D3_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D3_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D3_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D3_matrices, n = D3_n0, time = time) 
  D3_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D3_total_pop <- lapply(D3_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D3_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D3_total_pop[[i]])
  D3_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D3_plot_data <- bind_rows(D3_df_plots, .id = "id")

# Plot projection
D3_plot <- ggplot(D3_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D3_mean_plot_data <- D3_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D3_plot_pred <- D3_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D3_mean_plot <- ggplot(D3_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D3_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D3_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D3_total_pop[[i]][time]
  D3_pop_sizes[i] <- ms
}

# mean pop size
D3_pop_mean <- mean(D3_pop_sizes)

# standard deviation pop size
D3_pop_sd <- sd(D3_pop_sizes)

#standard error pop size
D3_pop_se <- sd(D3_pop_sizes)/sqrt(length(D3_pop_sizes))


#### Calculate Stochastic Growth Rate

D3_lambda_s <- stoch.growth.rate(D3_matrices, prob = NULL, maxt = time,
                                   verbose = TRUE)

#convert from log
D3_lambda_s$approx <- exp(D3_lambda_s$approx)
D3_lambda_s$sim <- exp(D3_lambda_s$sim)
D3_lambda_s$sim.CI <- exp(D3_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D3_quasi <- stoch.quasi.ext(D3_matrices, n0= D3_n0, Nx = 50, tmax = time, maxruns = 1,
                              nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D3_quasi_df <- data.frame(D3_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D3_quasi_plot <- ggplot(D3_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D3_sens <- stoch.sens(D3_matrices, tlimit=time)

D3_elas <- D3_sens$elasticities

D3_elas_v <- c(D3_elas[1,1], D3_elas[1,2], D3_elas[1,3], D3_elas[2,1], D3_elas[2,2], D3_elas[2,3], D3_elas[3,1], D3_elas[3,2], D3_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D3_elas_df <- data.frame(D3_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D3_elas_plot <- ggplot(D3_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 4 ----

## Stage duration
duration <- c(1, 4, 10)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D4_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D4_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D4_matrices[[i]] <- mpm
}


##Repeat Stochastic Population Growth 

D4_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D4_matrices, n = D4_n0, time = time) 
  D4_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D4_total_pop <- lapply(D4_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D4_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D4_total_pop[[i]])
  D4_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D4_plot_data <- bind_rows(D4_df_plots, .id = "id")

# Plot projection
D4_plot <- ggplot(D4_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D4_mean_plot_data <- D4_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D4_plot_pred <- D4_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D4_mean_plot <- ggplot(D4_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D4_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D4_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D4_total_pop[[i]][time]
  D4_pop_sizes[i] <- ms
}

# mean pop size
D4_pop_mean <- mean(D4_pop_sizes)

# standard deviation pop size
D4_pop_sd <- sd(D4_pop_sizes)

# standard error pop size
D4_pop_se <- sd(D4_pop_sizes)/sqrt(length(D4_pop_sizes))

#### Calculate Stochastic Growth Rate

D4_lambda_s <- stoch.growth.rate(D4_matrices, prob = NULL, maxt = time,
                                   verbose = TRUE)

#convert from log
D4_lambda_s$approx <- exp(D4_lambda_s$approx)
D4_lambda_s$sim <- exp(D4_lambda_s$sim)
D4_lambda_s$sim.CI <- exp(D4_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D4_quasi <- stoch.quasi.ext(D4_matrices, n0= D4_n0, Nx = 50, tmax = time, maxruns = 1,
                              nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D4_quasi_df <- data.frame(D4_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D4_quasi_plot <- ggplot(D4_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D4_sens <- stoch.sens(D4_matrices, tlimit=time)

D4_elas <- D4_sens$elasticities

D4_elas_v <- c(D4_elas[1,1], D4_elas[1,2], D4_elas[1,3], D4_elas[2,1], D4_elas[2,2], D4_elas[2,3], D4_elas[3,1], D4_elas[3,2], D4_elas[3,3])

D4_elas_df <- data.frame(D4_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D4_elas_plot <- ggplot(D4_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 5 ----

## Stage duration
duration <- c(1, 5, 10)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D5_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D5_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D5_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D5_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D5_matrices, n = D5_n0, time = time) 
  D5_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D5_total_pop <- lapply(D5_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D5_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D5_total_pop[[i]])
  D5_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D5_plot_data <- bind_rows(D5_df_plots, .id = "id")

# Plot projection
D5_plot <- ggplot(D5_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D5_mean_plot_data <- D5_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D5_plot_pred <- D5_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D5_mean_plot <- ggplot(D5_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D5_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D5_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D5_total_pop[[i]][time]
  D5_pop_sizes[i] <- ms
}

# mean pop size
D5_pop_mean <- mean(D5_pop_sizes)

# standard deviation pop size
D5_pop_sd <- sd(D5_pop_sizes)

# standard error pop size
D5_pop_se <- sd(D5_pop_sizes)/sqrt(length(D5_pop_sizes))

#### Calculate Stochastic Growth Rate

D5_lambda_s <- stoch.growth.rate(D5_matrices, prob = NULL, maxt = time,
                                   verbose = TRUE)

#convert from log
D5_lambda_s$approx <- exp(D5_lambda_s$approx)
D5_lambda_s$sim <- exp(D5_lambda_s$sim)
D5_lambda_s$sim.CI <- exp(D5_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D5_quasi <- stoch.quasi.ext(D5_matrices, n0= D5_n0, Nx = 50, tmax = time, maxruns = 1,
                              nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D5_quasi_df <- data.frame(D5_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D5_quasi_plot <- ggplot(D5_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D5_sens <- stoch.sens(D5_matrices, tlimit=time)

D5_elas <- D5_sens$elasticities

D5_elas_v <- c(D5_elas[1,1], D5_elas[1,2], D5_elas[1,3], D5_elas[2,1], D5_elas[2,2], D5_elas[2,3], D5_elas[3,1], D5_elas[3,2], D5_elas[3,3])

D5_elas_df <- data.frame(D5_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D5_elas_plot <- ggplot(D5_elas_df, aes(x = stage, y= elasticity)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### PLOTS ----

## Stochastic Population Projection Plot
stoch_plot <- D1_plot + D2_plot + D3_plot + D4_plot + D5_plot

#stoch_plot_same_lims <- D1_plot + ylim(0, 1.0e+07) + D2_plot + ylim(0, 1.0e+07) + D3_plot + ylim(0, 1.0e+07) + D4_plot + ylim(0, 1.0e+07) + D5_plot + ylim(0, 1.0e+07)

## Mean and CI Stochastic Population Plot
stoch_mean_plot <- D1_mean_plot + D2_mean_plot + D3_mean_plot + D4_mean_plot + D5_mean_plot


# Stochastic Population Growth (Lambda s)

stoch_lambda_approx <- c(D1_lambda_s$approx, D2_lambda_s$approx, D3_lambda_s$approx, D4_lambda_s$approx, D5_lambda_s$approx)

stoch_lambda_sim <- c(D1_lambda_s$sim, D2_lambda_s$sim, D3_lambda_s$sim, D4_lambda_s$sim, D5_lambda_s$sim)

stoch_lower_CI <- c(D1_lambda_s$sim.CI[1], D2_lambda_s$sim.CI[1], D3_lambda_s$sim.CI[1], D4_lambda_s$sim.CI[1], D5_lambda_s$sim.CI[1])

stoch_upper_CI <- c(D1_lambda_s$sim.CI[2], D2_lambda_s$sim.CI[2], D3_lambda_s$sim.CI[2], D4_lambda_s$sim.CI[2], D5_lambda_s$sim.CI[2])

stage_duration <- c("1 year", "2 years", "3 years", "4 years", "5 years")

stoch_lambda_df <- data.frame(stage_duration, stoch_lambda_approx, stoch_lambda_sim, stoch_upper_CI, stoch_lower_CI)

stoch_lambda_plot <- ggplot(stoch_lambda_df) +
  geom_point(aes(x = stage_duration, y = stoch_lambda_sim), fill = "grey20", size = 2) +
  geom_errorbar(aes(x = stage_duration, ymin = stoch_lower_CI, ymax = stoch_upper_CI), width = 0.2) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  scale_x_discrete(labels=c("1 year" = "1", "2 years" = "2", "3 years" = "3", "4 years" = "4", "5 years" = "5")) +
  labs(x = "Immature stage duration (years)", y = "Lambda for stochastic population growth")


## Quasi-extinction Threshold Plots

stoch_quasi_df<- rbind.data.frame(D1_quasi_df, D2_quasi_df, D3_quasi_df, D4_quasi_df, D5_quasi_df)

stoch_quasi_plot <- ggplot(stoch_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  scale_colour_discrete(name = "Immature stage \nduration", breaks = c("D1_quasi", "D2_quasi", "D3_quasi", "D4_quasi", "D5_quasi"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"))

#stoch_quasi_plots <- D1_quasi_plot + D2_quasi_plot + D3_quasi_plot + D4_quasi_plot + D5_quasi_plot


#Elasticity analysis plots

stoch_elas_df<- rbind.data.frame(D1_elas_df, D2_elas_df, D3_elas_df, D4_elas_df, D5_elas_df)

stoch_elas_plot <- ggplot(stoch_elas_df, aes(x = stage, y= elasticity, fill = duration)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(name = "Immature stage  \nduration", breaks = c("D1_elas_v", "D2_elas_v", "D3_elas_v", "D4_elas_v", "D5_elas_v"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"), values = c("grey65", "grey40", "grey35", "grey15", "grey0"))


#image2(D1_elas)
#image2(D2_elas)
#image2(D3_elas)
#image2(D4_elas)
#image2(D5_elas)
