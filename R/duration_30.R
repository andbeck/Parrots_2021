## YSA Stochastic Model of Population Growth
# With changes in immature stage duration of 1:5 years.
# With adult stage duration of 30:26 years.
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



#### LSA for Immature Duration of 1 and Adult Duration 0f 30----

## Stage duration
duration <- c(1, 1, 30)


## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution 
D1A30_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D1A30_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D1A30_matrices[[i]] <- mpm
}

head(D1A30_matrices)

## Repeat Stochastic Population Growth 

D1A30_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D1A30_matrices, n = D1A30_n0, time = time)
  D1A30_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D1A30_total_pop <- lapply(D1A30_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D1A30_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D1A30_total_pop[[i]])
  D1A30_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D1A30_plot_data <- bind_rows(D1A30_df_plots, .id = "id")

# Plot projection
D1A30_plot <- ggplot(D1A30_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D1A30_mean_plot_data <- D1A30_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D1A30_plot_pred <- D1A30_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D1A30_mean_plot <- ggplot(D1A30_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D1A30_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D1A30_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D1A30_total_pop[[i]][time]
  D1A30_pop_sizes[i] <- ms
}

#mean pop size
D1A30_pop_mean <- mean(D1A30_pop_sizes)

# standard deviation pop size
D1A30_pop_sd <- sd(D1A30_pop_sizes)

# standard error pop size
D1A30_pop_se <- sd(D1A30_pop_sizes)/sqrt(length(D1A30_pop_sizes))


#### Calculate Stochastic Growth Rate

D1A30_lambda_s <- stoch.growth.rate(D1A30_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
D1A30_lambda_s$approx <- exp(D1A30_lambda_s$approx)
D1A30_lambda_s$sim <- exp(D1A30_lambda_s$sim)
D1A30_lambda_s$sim.CI <- exp(D1A30_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D1A30_quasi <- stoch.quasi.ext(D1A30_matrices, n0= D1A30_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D1A30_quasi_df <- data.frame(D1A30_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D1A30_quasi_plot <- ggplot(D1A30_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D1A30_sens <- stoch.sens(D1A30_matrices, tlimit=time)

D1A30_elas <- D1A30_sens$elasticities

D1A30_elas_v <- c(D1A30_elas[1,1], D1A30_elas[1,2], D1A30_elas[1,3], D1A30_elas[2,1], D1A30_elas[2,2], D1A30_elas[2,3], D1A30_elas[3,1], D1A30_elas[3,2], D1A30_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D1A30_elas_df <- data.frame(D1A30_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D1A30_elas_plot <- ggplot(D1A30_elas_df, aes(x = stage, y= D1A30_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 2 and Adult Duration of 29 ----

## Stage duration
duration <- c(1, 2, 29)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

## Initial population vector estimated from stable stage distribution
D2A29_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D2A29_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D2A29_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D2A29_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D2A29_matrices, n = D2A29_n0, time = time) 
  D2A29_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D2A29_total_pop <- lapply(D2A29_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D2A29_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D2A29_total_pop[[i]])
  D2A29_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D2A29_plot_data <- bind_rows(D2A29_df_plots, .id = "id")

# Plot projection
D2A29_plot <- ggplot(D2A29_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D2A29_mean_plot_data <- D2A29_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D2A29_plot_pred <- D2A29_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D2A29_mean_plot <- ggplot(D2A29_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D2A29_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D2A29_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D2A29_total_pop[[i]][time]
  D2A29_pop_sizes[i] <- ms
}

# mean pop size
D2A29_pop_mean <- mean(D2A29_pop_sizes)

# standard deviation pop size
D2A29_pop_sd <- sd(D2A29_pop_sizes)

#standard error pop size
D2A29_pop_se <- sd(D2A29_pop_sizes)/sqrt(length(D2A29_pop_sizes))


#### Calculate Stochastic Growth Rate

D2A29_lambda_s <- stoch.growth.rate(D2A29_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
D2A29_lambda_s$approx <- exp(D2A29_lambda_s$approx)
D2A29_lambda_s$sim <- exp(D2A29_lambda_s$sim)
D2A29_lambda_s$sim.CI <- exp(D2A29_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D2A29_quasi <- stoch.quasi.ext(D2A29_matrices, n0= D2A29_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D2A29_quasi_df <- data.frame(D2A29_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D2A29_quasi_plot <- ggplot(D2A29_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D2A29_sens <- stoch.sens(D2A29_matrices, tlimit=time)

D2A29_elas <- D2A29_sens$elasticities

D2A29_elas_v <- c(D2A29_elas[1,1], D2A29_elas[1,2], D2A29_elas[1,3], D2A29_elas[2,1], D2A29_elas[2,2], D2A29_elas[2,3], D2A29_elas[3,1], D2A29_elas[3,2], D2A29_elas[3,3])

stage<-c("m1", "m2", "m3", "s1", "s2", "s3", "g1", "g2", "s3")

D2A29_elas_df <- data.frame(D2A29_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D2A29_elas_plot <- ggplot(D2A29_elas_df, aes(x = stage, y= D2A29_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 3 and Adult Duration of 28 ----

## Stage duration
duration <- c(1, 3, 28)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D3A28_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D3A28_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D3A28_matrices[[i]] <- mpm
}


##Repeat Stochastic Population Growth 

D3A28_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D3A28_matrices, n = D3A28_n0, time = time) 
  D3A28_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D3A28_total_pop <- lapply(D3A28_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D3A28_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D3A28_total_pop[[i]])
  D3A28_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D3A28_plot_data <- bind_rows(D3A28_df_plots, .id = "id")

# Plot projection
D3A28_plot <- ggplot(D3A28_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D3A28_mean_plot_data <- D3A28_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D3A28_plot_pred <- D3A28_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D3A28_mean_plot <- ggplot(D3A28_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D3A28_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D3A28_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D3A28_total_pop[[i]][time]
  D3A28_pop_sizes[i] <- ms
}

# mean pop size
D3A28_pop_mean <- mean(D3A28_pop_sizes)

# standard deviation pop size
D3A28_pop_sd <- sd(D3A28_pop_sizes)

# standard error pop size
D3A28_pop_se <- sd(D3A28_pop_sizes)/sqrt(length(D3A28_pop_sizes))

#### Calculate Stochastic Growth Rate

D3A28_lambda_s <- stoch.growth.rate(D3A28_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
D3A28_lambda_s$approx <- exp(D3A28_lambda_s$approx)
D3A28_lambda_s$sim <- exp(D3A28_lambda_s$sim)
D3A28_lambda_s$sim.CI <- exp(D3A28_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D3A28_quasi <- stoch.quasi.ext(D3A28_matrices, n0= D3A28_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D3A28_quasi_df <- data.frame(D3A28_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D3A28_quasi_plot <- ggplot(D3A28_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D3A28_sens <- stoch.sens(D3A28_matrices, tlimit=time)

D3A28_elas <- D3A28_sens$elasticities

D3A28_elas_v <- c(D3A28_elas[1,1], D3A28_elas[1,2], D3A28_elas[1,3], D3A28_elas[2,1], D3A28_elas[2,2], D3A28_elas[2,3], D3A28_elas[3,1], D3A28_elas[3,2], D3A28_elas[3,3])

D3A28_elas_df <- data.frame(D3A28_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D3A28_elas_plot <- ggplot(D3A28_elas_df, aes(x = stage, y= D3A28_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")


#### LSA for Immature Duration of 4 and Adult Duration of 27 ----

## Stage duration
duration <- c(1, 4, 27)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D4A27_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D4A27_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D4A27_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D4A27_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D4A27_matrices, n = D4A27_n0, time = time) 
  D4A27_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D4A27_total_pop <- lapply(D4A27_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D4A27_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D4A27_total_pop[[i]])
  D4A27_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D4A27_plot_data <- bind_rows(D4A27_df_plots, .id = "id")

# Plot projection
D4A27_plot <- ggplot(D4A27_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D4A27_mean_plot_data <- D4A27_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D4A27_plot_pred <- D4A27_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D4A27_mean_plot <- ggplot(D4A27_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D4A27_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D4A27_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D4A27_total_pop[[i]][time]
  D4A27_pop_sizes[i] <- ms
}

# mean pop size
D4A27_pop_mean <- mean(D4A27_pop_sizes)

# standard deviation pop size
D4A27_pop_sd <- sd(D4A27_pop_sizes)

# standard error pop size
D4A27_pop_se <- sd(D4A27_pop_sizes)/sqrt(length(D4A27_pop_sizes))

#### Calculate Stochastic Growth Rate

D4A27_lambda_s <- stoch.growth.rate(D4A27_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
D4A27_lambda_s$approx <- exp(D4A27_lambda_s$approx)
D4A27_lambda_s$sim <- exp(D4A27_lambda_s$sim)
D4A27_lambda_s$sim.CI <- exp(D4A27_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D4A27_quasi <- stoch.quasi.ext(D4A27_matrices, n0= D4A27_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D4A27_quasi_df <- data.frame(D4A27_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D4A27_quasi_plot <- ggplot(D4A27_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D4A27_sens <- stoch.sens(D4A27_matrices, tlimit=time)

D4A27_elas <- D4A27_sens$elasticities

D4A27_elas_v <- c(D4A27_elas[1,1], D4A27_elas[1,2], D4A27_elas[1,3], D4A27_elas[2,1], D4A27_elas[2,2], D4A27_elas[2,3], D4A27_elas[3,1], D4A27_elas[3,2], D4A27_elas[3,3])

D4A27_elas_df <- data.frame(D4A27_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D4A27_elas_plot <- ggplot(D4A27_elas_df, aes(x = stage, y= D4A27_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")



#### LSA for Immature Duration of 5 and Adult Duration of 26 ----

## Stage duration
duration <- c(1, 5, 26)

## Initial Population Vector

# Stable stage distribution of mean matrix
stable_stage <- make_projection_matrix(survival, fecundity, duration) %>%
  stable.stage() %>% as.list()

# Initial population vector estimated from stable stage distribution
D5A26_n0 <- c(stable_stage[[1]]*Nc, stable_stage[[2]]*Nc, stable_stage[[3]]*Nc)

### Life-stage Simulation Analysis for Population in Stochastic Environment

## Stage duration list - repeat so that length is the same as survival and fecundity
duration_list <- rep(list(duration), times = n_sim)


## Simulate list of matrices using the vital rates and make_projection_matrix function

D5A26_matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  D5A26_matrices[[i]] <- mpm
}


## Repeat Stochastic Population Growth 

D5A26_stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(D5A26_matrices, n = D5A26_n0, time = time) 
  D5A26_stochastic_pop[i] <- mp
}

# Multiply female population sizes by 2 to get total population size
D5A26_total_pop <- lapply(D5A26_stochastic_pop, "*", 2)


# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

D5A26_df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = D5A26_total_pop[[i]])
  D5A26_df_plots[[i]] <- mpl
}

# Add identifier for each simulation
D5A26_plot_data <- bind_rows(D5A26_df_plots, .id = "id")

# Plot projection
D5A26_plot <- ggplot(D5A26_plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Total population size")


# Mean population size time series with 95% confidence intervals from LSA

D5A26_mean_plot_data <- D5A26_plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
D5A26_plot_pred <- D5A26_mean_plot_data %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D5A26_mean_plot <- ggplot(D5A26_plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = D5A26_plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_classic() +
  labs(x = "Time (years)", y = "Mean total population size")


#### Calculate final mean population size and standard deviation from LSA

D5A26_pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- D5A26_total_pop[[i]][time]
  D5A26_pop_sizes[i] <- ms
}

# mean pop size
D5A26_pop_mean <- mean(D5A26_pop_sizes)

# standard deviation pop size
D5A26_pop_sd <- sd(D5A26_pop_sizes)

# standard error pop size
D5A26_pop_se <- sd(D5A26_pop_sizes)/sqrt(length(D5A26_pop_sizes))

#### Calculate Stochastic Growth Rate

D5A26_lambda_s <- stoch.growth.rate(D5A26_matrices, prob = NULL, maxt = time,
                                    verbose = TRUE)

#convert from log
D5A26_lambda_s$approx <- exp(D5A26_lambda_s$approx)
D5A26_lambda_s$sim <- exp(D5A26_lambda_s$sim)
D5A26_lambda_s$sim.CI <- exp(D5A26_lambda_s$sim.CI)

#### Calculate Quasi-extinction Probability

D5A26_quasi <- stoch.quasi.ext(D5A26_matrices, n0= D5A26_n0, Nx = 50, tmax = time, maxruns = 1,
                               nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D5A26_quasi_df <- data.frame(D5A26_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D5A26_quasi_plot <- ggplot(D5A26_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  labs(y = "Cumulative probability of quasi-extinction")


#### Calculate Stochastic Elasticities

D5A26_sens <- stoch.sens(D5A26_matrices, tlimit=time)

D5A26_elas <- D5A26_sens$elasticities

D5A26_elas_v <- c(D5A26_elas[1,1], D5A26_elas[1,2], D5A26_elas[1,3], D5A26_elas[2,1], D5A26_elas[2,2], D5A26_elas[2,3], D5A26_elas[3,1], D5A26_elas[3,2], D5A26_elas[3,3])

D5A26_elas_df <- data.frame(D5A26_elas_v) %>% gather("duration", "elasticity") %>% data.frame(stage)

D5A26_elas_plot <- ggplot(D5A26_elas_df, aes(x = stage, y= D5A26_elas_v)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(fill = "grey20")



#### PLOTS ----

## Stochastic Population Projection Plot
A30_plot <- D1A30_plot + D2A29_plot + D3A28_plot + D4A27_plot + D5A26_plot

## Mean and CI Stochastic Population Plot
A30_mean_plot <- D1A30_mean_plot + D2A29_mean_plot + D3A28_mean_plot + D4A27_mean_plot + D5A26_mean_plot


# Stochastic Population Growth (Lambda s)

A30_lambda_approx <- c(D1A30_lambda_s$approx, D2A29_lambda_s$approx, D3A28_lambda_s$approx, D4A27_lambda_s$approx, D5A26_lambda_s$approx)

A30_lambda_sim <- c(D1A30_lambda_s$sim, D2A29_lambda_s$sim, D3A28_lambda_s$sim, D4A27_lambda_s$sim, D5A26_lambda_s$sim)

A30_lower_CI <- c(D1A30_lambda_s$sim.CI[1], D2A29_lambda_s$sim.CI[1], D3A28_lambda_s$sim.CI[1], D4A27_lambda_s$sim.CI[1], D5A26_lambda_s$sim.CI[1])

A30_upper_CI <- c(D1A30_lambda_s$sim.CI[2], D2A29_lambda_s$sim.CI[2], D3A28_lambda_s$sim.CI[2], D4A27_lambda_s$sim.CI[2], D5A26_lambda_s$sim.CI[2])

stage_duration <- c("1 year", "2 years", "3 years", "4 years", "5 years")

A30_lambda_df <- data.frame(stage_duration, A30_lambda_approx, A30_lambda_sim, A30_upper_CI, A30_lower_CI)

A30_lambda_plot <- ggplot(A30_lambda_df) +
  geom_point(aes(x = stage_duration, y = A30_lambda_sim), fill = "grey20", size = 2) +
  geom_errorbar(aes(x = stage_duration, ymin = A30_lower_CI, ymax = A30_upper_CI), width = 0.2) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", colour = "red") +
  scale_x_discrete(labels=c("1 year" = "1", "2 years" = "2", "3 years" = "3", "4 years" = "4", "5 years" = "5")) +
  labs(x = "Immature stage duration (years)", y = "Lambda for stochastic population growth")


## Quasi-extinction Threshold Plots

A30_quasi_df<- rbind.data.frame(D1A30_quasi_df, D2A29_quasi_df, D3A28_quasi_df, D4A27_quasi_df, D5A26_quasi_df)

A30_quasi_plot <- ggplot(A30_quasi_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  labs(y = "Cumulative probability of quasi-extinction") +
  scale_colour_discrete(name = "Immature stage \nduration", breaks = c("D1A30_quasi", "D2A29_quasi", "D3A28_quasi", "D4A27_quasi", "D5A26_quasi"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"))

#A30_quasi_plots <- D1A30_quasi_plot + D2A29_quasi_plot + D3A28_quasi_plot + D4A27_quasi_plot + D5A26_quasi_plot


#Elasticity analysis plots

A30_elas_df <- rbind.data.frame(D1A30_elas_df, D2A29_elas_df, D3A28_elas_df, D4A27_elas_df, D5A26_elas_df)

A30_elas_plot <- ggplot(A30_elas_df, aes(x = stage, y= elasticity, fill = duration)) + 
  labs(x = "Matrix element", y = "Stochastic elasticity") +
  theme_bw() +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(name = "Immature stage \nduration", breaks = c("D1A30_elas_v", "D2A29_elas_v", "D3A28_elas_v", "D4A27_elas_v", "D5A26_elas_v"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"), values = c("grey65", "grey40", "grey35", "grey15", "grey0"))

#image2(D1A30_elas)
#image2(D2A29_elas)
#image2(D3A28_elas)
#image2(D4A27_elas)
#image2(D5A26_elas)

