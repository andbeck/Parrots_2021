#### Prepare data ----

# 1. YSA breeding biology data
# data/Life History 2006-2014.xlsx - this is YSA breeding biology data from Bonaire
# R/YSA_life_history_data.R - script which read and summarises the YSA breeding biology data

source("R/YSA_life_history_data.R")

# 2. Combine with estimates of survival and stage duration

# Manually creating data frame for now
stage <- c("1", "2", "3")

# Stage classes
class <- c("fledgling", "juvenile", "adult")

# Stage durations
# To age 12 months, Age 13-36 months,
# Age 37 months+ (as 10 years -> needs updating and a reference)
duration <- c(1, (24/12), 10)

# 0.71 from Salinas-Melgoza & Renton 2007 (CHECK), 0.875 from Rodriguez et al 2004
survival <- c(0.71, 0.875, 0.875)

# SD/standard error estimates from YSA_demog_data_master.csv (CHECK), 0.2 from
# Salinas-Melgoza & Renton 2007, 0.075 from Rodriguez et al 2004
survival_SD <- c(0.2, 0.075, 0.075)

# Reproductive output/fecundity
# Clutch size taken as half of 3.2 as sex ratio assumed 1:1 (Sam's thesis), data
# from YSA_demog_data_master_csv
fecundity <- c(0, 0, 1.6*total_summary$mean_hatch[1]*total_summary$mean_nestling_surv[1])

# Reproductive output/fecundity SEs
# SE Sam's thesis (halved) TODO modify to incorporate SEs for hatch and nestling surv:
# total_summary$se_hatch[1], total_summary$se_nestling_surv[1]
fecundity_SD <- c(0, 0, 0.1)

#### Build matrix model ----

source("R/make_projection_matrix.R")

ysa <- make_projection_matrix(survival, fecundity, duration)

lambda <- eigen(ysa)$values[1]

#### Life stage simulation analysis ----

library(popbio)

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

#make survival data frame 
survival_df <- data.frame(s1, s2, s3)
colnames(survival_df)<- c()

#make fecundity data frame
fecundity_df <- data.frame(0, 0, m3)
colnames(fecundity_df)<- c()

# Create a list of vital rates
duration_list <- rep(list(c(1, 24/12, 10)), times = n) #repeat so that length is the same as survival and fecundity

survival_list <- asplit(survival_df, 1)

fecundity_list <- asplit(fecundity_df, 1)

## Create simulated matrices using for loop
matrices <- list()

for(i in 1:n){
  mpm <- make_projection_matrix(survival_list[[i]], fecundity_list[[i]], duration_list[[i]])
  matrices[[i]] <- mpm
}

matrices

## Simulate population

# time 
t <- 100

#for loop for applying pop.projection from popbio package to each of the matrices

proj_matrices <- list()

for(i in 1:n){
  mp <- pop.projection(A = matrices[[i]],
                       n = c(30, 40, 70),
                       iterations = t)
  proj_matrices[[i]] <- mp
}

#to plot with ggplot, create for loop for pop sizes in each matrix as a data frame
df_plots <- list()

for(i in 1:n){
  mpl <- data.frame(time = 1:t, pop_sizes = proj_matrices[[i]]$pop.sizes)
  df_plots[[i]] <- mpl
}

#add identifier for each matrix to plot each simulation
plot_data <- bind_rows(df_plots, .id = "id")
glimpse(plot_data)

#plot
ggplot(plot_data, aes(time, pop_sizes, colour=id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

## Calculate lambdas ----

ysa_lambdas <- sapply(matrices, lambda)
hist(ysa_lambdas)

#mean lambda
ysa_lambda_mean <- mean(ysa_lambdas)
ysa_lambda_mean

#standard error of lambda
ysa_lambda_se <- sd(ysa_lambdas)/sqrt(length(ysa_lambdas))

#elasticities
ysa_elas <- lapply(matrices, elasticity)

#mean elasticity
ysa_elas_mean <- Reduce('+', ysa_elas)/length(ysa_elas)

image2(ysa_elas_mean)

## Needs nest limit and some way to account for large non-breeding population ----
