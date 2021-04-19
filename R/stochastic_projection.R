## YSA Stochastic Model of Population Growth
#With imputed survival rates.
#Density independent

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

# Current population is estimated around 1000 individuals
Nc <- 1000

# Time to project to
time <- 100

#### YSA Simulated Vital Rates for LSA ----

set.seed(2021)

# Number of simulations
n_sim <- 1000

# Fledgling survival
s1 <- sapply(1:n_sim, function(x) betaval(0.71, 0.2))

# Juvenile survival
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



#### LSA for Juvenile Duration of 1 ----

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
duration_list <- rep(list(c(1, 1, 10)), times = n_sim)

## Simulate list of matrices using the vital rates and make_projection_matrix function

matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  matrices[[i]] <- mpm
}

head(matrices)

## Repeat Stochastic Population Growth 

stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(matrices, n = D1_n0, time = time) 
  stochastic_pop[[i]] <- mp
}

head(stochastic_pop)

# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = stochastic_pop[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

# Add identifier for each simulation
plot_data <- bind_rows(df_plots, .id = "id")

# Plot projection
D1_plot <- ggplot(plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Population size") +
  ggtitle("Immature stage: 1 year")


# Mean population size time series with 95% confidence intervals from LSA

plot <- plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
plot_pred <- plot %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D1_mean_plot <- ggplot(plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_bw() +
  labs(x = "Time (years)", y = "Mean population size") +
  ggtitle("Immature stage: 1 year")


#### Calculate final mean population size and standard deviation from LSA

pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- stochastic_pop[[i]]$pop.sizes[time]
  pop_sizes[i] <- ms
}

# mean pop size
D1_pop_mean <- mean(pop_sizes)

# standard deviation pop size
D1_pop_sd <- sd(pop_sizes)

# standard error pop size
D1_pop_se <- sd(pop_sizes)/sqrt(length(pop_sizes))


#### Calculate Stochastic Growth Rate

D1_lambda_s <- stoch.growth.rate(matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)


#### Calculate Quasi-extinction Probability

D1_quasi <- stoch.quasi.ext(matrices, n0= D1_n0, Nx = 50, tmax = time, maxruns = 10,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D1_df <- data.frame(D1_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D1_quasi_plot <- ggplot(D1_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Quasi-extinction probability") +
  ggtitle("Immature stage: 1 year")


#### Calculate Stochastic Sensitivities

#gives the sensitivties and elasticities for matrix elements, not vital rates? 
D1_elasticity <- stoch.sens(matrices, tlimit = time)


#### LSA for Juvenile Duration of 2 ----

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

matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  matrices[[i]] <- mpm
}

head(matrices)

## Repeat Stochastic Population Growth 

stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(matrices, n = D2_n0, time = time) 
  stochastic_pop[[i]] <- mp
}

head(stochastic_pop)

# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = stochastic_pop[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

# Add identifier for each simulation
plot_data <- bind_rows(df_plots, .id = "id")

# Plot projection
D2_plot <- ggplot(plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Population size") +
  ggtitle("Immature stage: 2 years")


# Mean population size time series with 95% confidence intervals from LSA

plot <- plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
plot_pred <- plot %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D2_mean_plot <- ggplot(plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_bw() +
  labs(x = "Time (years)", y = "Mean Population Size") +
  ggtitle("Immature stage: 2 years")


#### Calculate final mean population size and standard deviation from LSA

pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- stochastic_pop[[i]]$pop.sizes[time]
  pop_sizes[i] <- ms
}

#mean pop size
D2_pop_mean <- mean(pop_sizes)

# standard deviation pop size
D2_pop_sd <- sd(pop_sizes)

# standard error pop size
D2_pop_se <- sd(pop_sizes)/sqrt(length(pop_sizes))


#### Calculate Stochastic Growth Rate

D2_lambda_s <- stoch.growth.rate(matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)


#### Calculate Quasi-extinction Probability

D2_quasi <- stoch.quasi.ext(matrices, n0= D2_n0, Nx = 50, tmax = time, maxruns = 10,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D2_df <- data.frame(D2_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D2_quasi_plot <- ggplot(D2_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Quasi-extinction probability") +
  ggtitle("Immature stage: 2 years")


#### Calculate Stochastic Sensitivities

D2_elasticity <- stoch.sens(matrices, tlimit=time)



#### LSA for Juvenile Duration of 3 ----

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

matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  matrices[[i]] <- mpm
}

head(matrices)

## Repeat Stochastic Population Growth 

stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(matrices, n = D3_n0, time = time) 
  stochastic_pop[[i]] <- mp
}

head(stochastic_pop)

# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = stochastic_pop[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

# Add identifier for each simulation
plot_data <- bind_rows(df_plots, .id = "id")

# Plot projection
D3_plot <- ggplot(plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Population size") +
  ggtitle("Immature stage: 3 years")


# Mean population size time series with 95% confidence intervals from LSA

plot <- plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
plot_pred <- plot %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D3_mean_plot <- ggplot(plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_bw() +
  labs(x = "Time (years)", y = "Mean population size") + 
  ggtitle("Immature stage: 3 years")


#### Calculate final mean population size and standard deviation from LSA

pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- stochastic_pop[[i]]$pop.sizes[time]
  pop_sizes[i] <- ms
}

# mean pop size
D3_pop_mean <- mean(pop_sizes)

# standard deviation pop size
D3_pop_sd <- sd(pop_sizes)

#standard error pop size
D3_pop_se <- sd(pop_sizes)/sqrt(length(pop_sizes))


#### Calculate Stochastic Growth Rate

D3_lambda_s <- stoch.growth.rate(matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)


#### Calculate Quasi-extinction Probability

D3_quasi <- stoch.quasi.ext(matrices, n0= D3_n0, Nx = 50, tmax = time, maxruns = 10,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D3_df <- data.frame(D3_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D3_quasi_plot <- ggplot(D3_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Quasi-extinction probability") +
  ggtitle("Immature stage: 3 years")


#### Calculate Stochastic Sensitivities

D3_elasticity <- stoch.sens(matrices, tlimit = time)



#### LSA for Juvenile Duration of 4 ----

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

matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  matrices[[i]] <- mpm
}

head(matrices)

##Repeat Stochastic Population Growth 

stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(matrices, n = D4_n0, time = time) 
  stochastic_pop[[i]] <- mp
}

head(stochastic_pop)

# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = stochastic_pop[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

# Add identifier for each simulation
plot_data <- bind_rows(df_plots, .id = "id")

# Plot projection
D4_plot <- ggplot(plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Population size") +
  ggtitle("Immature stage: 4 years")


# Mean population size time series with 95% confidence intervals from LSA

plot <- plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
plot_pred <- plot %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D4_mean_plot <- ggplot(plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_bw() +
  labs(x = "Time (years)", y = "Mean population size") +
  ggtitle("Immature stage: 4 years")


#### Calculate final mean population size and standard deviation from LSA

pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- stochastic_pop[[i]]$pop.sizes[time]
  pop_sizes[i] <- ms
}

# mean pop size
D4_pop_mean <- mean(pop_sizes)

# standard deviation pop size
D4_pop_sd <- sd(pop_sizes)

# standard error pop size
D4_pop_se <- sd(pop_sizes)/sqrt(length(pop_sizes))

#### Calculate Stochastic Growth Rate

D4_lambda_s <- stoch.growth.rate(matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)


#### Calculate Quasi-extinction Probability

D4_quasi <- stoch.quasi.ext(matrices, n0= D4_n0, Nx = 50, tmax = time, maxruns = 10,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D4_df <- data.frame(D4_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D4_quasi_plot <- ggplot(D4_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Quasi-extinction probability") +
  ggtitle("Immature stage: 4 years")


#### Calculate Stochastic Sensitivities

D4_elasticity <- stoch.sens(matrices, time)



#### LSA for Juvenile Duration of 5 ----

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

matrices <- list()

for(i in 1:n_sim){
  mpm <- make_projection_matrix(survival_list[[i]], 
                                fecundity_list[[i]], 
                                duration_list[[i]])
  matrices[[i]] <- mpm
}

head(matrices)

## Repeat Stochastic Population Growth 

stochastic_pop <- list()

for(i in 1:n_sim){
  mp <- stochastic_proj(matrices, n = D5_n0, time = time) 
  stochastic_pop[[i]] <- mp
}

head(stochastic_pop)

# Create for loop for pop sizes in each projection as a data frame to plot with ggplot

df_plots <- list()

for(i in 1:n_sim){
  mpl <- data.frame(time = 1:time, pop_sizes = stochastic_pop[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

# Add identifier for each simulation
plot_data <- bind_rows(df_plots, .id = "id")

# Plot projection
D5_plot <- ggplot(plot_data, aes(time, pop_sizes, fill=id)) +
  geom_line() +
  theme_classic() +
  labs(x = "Time (years)", y = "Population size") +
  ggtitle("Immature stage: 5 years")


# Mean population size time series with 95% confidence intervals from LSA

plot <- plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# Get predictions and 95% CI
plot_pred <- plot %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

# Plot mean population projection with CIs
D5_mean_plot <- ggplot(plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Time (years)", y = "Mean population size") +
  ggtitle("Immature stage: 5 years")


#### Calculate final mean population size and standard deviation from LSA

pop_sizes <- numeric()
for (i in 1:n_sim) {
  ms <- stochastic_pop[[i]]$pop.sizes[time]
  pop_sizes[i] <- ms
}

# mean pop size
D5_pop_mean <- mean(pop_sizes)

# standard deviation pop size
D5_pop_sd <- sd(pop_sizes)

# standard error pop size
D5_pop_se <- sd(pop_sizes)/sqrt(length(pop_sizes))

#### Calculate Stochastic Growth Rate

D5_lambda_s <- stoch.growth.rate(matrices, prob = NULL, maxt = time,
                                 verbose = TRUE)


#### Calculate Quasi-extinction Probability

D5_quasi <- stoch.quasi.ext(matrices, n0= D5_n0, Nx = 50, tmax = time, maxruns = 10,
                            nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)

# Plot quasi-extinction probabilities
D5_df <- data.frame(D5_quasi, "Year" = 1:time) %>% gather("sim", "quasi", -"Year")

D5_quasi_plot <- ggplot(D5_df, aes(x = Year, y = quasi, colour = sim)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Quasi-extinction probability") +
  ggtitle("Immature stage: 5 years")


#### Calculate Stochastic Sensitivities

D5_elasticity <- stoch.sens(matrices, time)


#### OUTPUTS ----


## Stochastic Population Projection Plot
D1_plot + D2_plot + D3_plot + D4_plot + D5_plot


## Mean and CI Stochastic Population Plot
D1_mean_plot + D2_mean_plot + D3_mean_plot + D4_mean_plot + D5_mean_plot


# Quasi Extinction Threshold Plot
D1_quasi_plot + D2_quasi_plot + D3_quasi_plot + D4_quasi_plot + D5_quasi_plot


# Stochastic Population Growth (Lambda s)

lambda_approx <- c(D1_lambda_s$approx, D2_lambda_s$approx, D3_lambda_s$approx, D4_lambda_s$approx, D5_lambda_s$approx)

lambda_sim <- c(D1_lambda_s$sim, D2_lambda_s$sim, D3_lambda_s$sim, D4_lambda_s$sim, D5_lambda_s$sim)

lower_CI <- c(D1_lambda_s$sim.CI[1], D2_lambda_s$sim.CI[1], D3_lambda_s$sim.CI[1], D4_lambda_s$sim.CI[1], D5_lambda_s$sim.CI[1])

upper_CI <- c(D1_lambda_s$sim.CI[2], D2_lambda_s$sim.CI[2], D3_lambda_s$sim.CI[2], D4_lambda_s$sim.CI[2], D5_lambda_s$sim.CI[2])

stage_duration <- c("1 year", "2 years", "3 years", "4 years", "5 years")

lambda_df <- data.frame(stage_duration, lambda_approx, lambda_sim, upper_CI, lower_CI)

lambda_plot <- ggplot(lambda_df) +
  geom_col( aes(x = stage_duration, y = lambda_sim)) +
  geom_errorbar(aes(x = stage_duration, ymin = lower_CI, ymax = upper_CI)) +
  theme_bw()

lambda_plot


#Elasticity analysis plots

image2(D1_elasticity$elasticities)
image2(D2_elasticity$elasticities)
image2(D3_elasticity$elasticities)
image2(D4_elasticity$elasticities)
image2(D5_elasticity$elasticities)

