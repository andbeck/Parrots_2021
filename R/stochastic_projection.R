## YSA STOCHASTIC MODEL OF POPULATION GROWTH

#### Libraries ----
library(popbio)
library(tidyverse)

#### Functions ----

## matrix model function
source("R/make_projection_matrix.R")

## stochastic population growth function
source("R/stochastic_proj.R")

#### YSA Data ----

## YSA breeding biology data 2006-2014 from Bonaire
source("R/YSA_life_history_data.R")

## YSA vital rates

set.seed(2021)

## Sample from vital rates

#number of simulations
n <- 1000

#fledgling survival
s1 <- sapply(1:n, function(x) betaval(0.71, 0.2))

#juvenile survival
s2 <- sapply(1:n, function(x) betaval(0.875, 0.075))

#adult survival
s3 <- sapply(1:n, function(x) betaval(0.875, 0.075))

#fecundity
m3 <- rlnorm(n = n, log(1.6*total_summary$mean_hatch[1]*total_summary$mean_nestling_surv[1]), log(1.01)) #replaced sd with small value for log
hist(m3)             

## Create lists of survival, fecundity, and duration

# survival 
survival_df <- data.frame(s1, s2, s3)
colnames(survival_df)<- c()
survival_list <- asplit(survival_df, 1)

# fecundity 
fecundity_df <- data.frame(0, 0, m3)
colnames(fecundity_df)<- c()
fecundity_list <- asplit(fecundity_df, 1)

# duration list
duration_list <- rep(list(c(1, 24/12, 10)), times = n) #repeat so that length is the same as survival and fecundity


#### Simulate list of matrices using the vital rates and make_projection_matrix function ----

matrices <- list()

for(i in 1:n){
  mpm <- make_projection_matrix(survival_list[[i]], fecundity_list[[i]], duration_list[[i]])
  matrices[[i]] <- mpm
}

head(matrices)


#### Life Stage Simulation Analysis for Population in Stochastic Environment ----

#check number of simulations
n

#initial population sizes
n0 <- c(20, 20, 20) #(chosen a bit randomly for now)

##repeat stochastic population growth 

proj_matrices <- list()

for(i in 1:n){
  mp <- stochastic_proj(matrices, n = n0, smax = 70, nmax = 10000, time = 100) #smax initially set to 70 (number of nests in Sam's thesis)
  proj_matrices[[i]] <- mp
}

head(proj_matrices)

#to plot with ggplot, create for loop for pop sizes in each matrix as a data frame

time <- 100

df_plots <- list()

for(i in 1:n){
  mpl <- data.frame(time = 1:time, pop_sizes = proj_matrices[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

#add identifier for each simulation
plot_data <- bind_rows(df_plots, .id = "id")

#plot
ggplot(plot_data, aes(time, pop_sizes, colour=id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")


## Plot mean population size time series with 95% confidence intervals from LSA ----

plot <- plot_data %>% 
  group_by(time) %>% 
  summarise(mean = mean(pop_sizes),
            se_pop_size = sd(pop_sizes)/sqrt(length(pop_sizes)))

# get predictions and 95% CI
plot_pred <- plot %>% mutate(
  pop_size = mean,  
  # lower limit 95% CI
  ll = mean - 1.96 * se_pop_size,
  # upper limit 95% CI
  ul = mean + 1.96 * se_pop_size
)

#plot
ggplot(plot_pred, aes(x= time, y = mean)) +
  geom_line() +
  geom_ribbon(data = plot_pred, aes(ymin = ll, ymax = ul), alpha = 0.2) + 
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Time (years)", y = "Mean Population Size")


#### Calculating final mean population size and standard deviation from LSA -----

pop_sizes <- numeric()
for (i in 1:n) {
  ms <- proj_matrices[[i]]$pop.sizes[100]
  pop_sizes[i] <- ms
}

pop_mean <- mean(pop_sizes)
pop_mean

pop_sd <- sd(pop_sizes)
pop_sd

pop_se <- sd(pop_sizes)/sqrt(length(pop_sizes))
pop_se

#### Calculating stochastic growth rate ----

#### Calculating cumulative distribution function (CDF) ----
#section unfinished, mostly playing with code from popbio package
quasi <- stoch.quasi.ext(matrices, n0= c(20, 20, 20), Nx = 50, tmax = 100, maxruns = 10,
                         nreps = 5000, prob = NULL, sumweight = NULL, verbose = TRUE)


matplot(quasi, xlab="Years", ylab="Quasi-extinction probability",
        type='l', lty=1, col=rainbow(10), las=1,
        main="Time to reach a quasi-extinction threshold of 50 individuals")


stoch.growth.rate(matrices, prob = NULL, maxt = 100,
                  verbose = TRUE)
