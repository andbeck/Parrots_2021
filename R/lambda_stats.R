#New lambda values and statistical analysis

#libraries
library(tidyverse)
library(ggfortify)
library(ggsignif)
library(car)
library(gmodels)

#matrices
source("R/duration_10.R")
source("R/duration_20.R")
source("R/duration_30.R")
source("R/poached10_A20.R")
source("R/poached20_A20.R")
source("R/poached30_A20.R")
source("R/poached40_A20.R")
source("R/poached50_A20.R")


#Stochastic lambda function ----
stoch.lambda <- function(maxt = 1000, matrices){
  matrices <- matrix(unlist(matrices), ncol = length(matrices))    
  
  s <- sqrt(dim(matrices)[1])
  n <- dim(matrices)[2]
  
  prob <- rep(1/n, n)
  
  Abar <- numeric(s^2)
  Exy <- numeric(s^4)
  for (i in 1:n) {
    A <- matrices[, i]
    Exy <- Exy + prob[i] * kronecker(A, A)
    Abar <- Abar + prob[i] * A
  }
  
  C <- (Exy - kronecker(Abar, Abar)) * n/(n - 1)
  C <- matrix(C, nrow = s^2)
  Abar <- matrix(Abar, nrow = s)
  ev <- eigen(Abar)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  
  n0 <- w
  
  r <- numeric(maxt)
  for (t in 1:maxt) {
    if (t == 1 || t%%1 == 0) {
      message("Calculating stochastic growth at time ", 
              t)
    }
    col <- sample(1:n, 1, prob = prob)
    A <- matrix(matrices[, col], nrow = s)
    n0 <- A %*% n0
    N <- sum(n0)
    r[t] <- log(N)
    n0 <- n0/N
  }
  loglsim <- mean(r)
  dse <- 1.96 * sqrt(var(r)/maxt)
  CI <- c(loglsim - dse, loglsim + dse)
  stoch <- list(lambda = exp(r), sim = exp(loglsim), 
                sim.CI = exp(CI))
  return(stoch)
}

#number of time intervals
maxt <- 1000


##DURATION 10 ----

##D1A8 
D1A8_stoch <- stoch.lambda(matrices = D1A8_matrices)

##D2A7 
D2A7_stoch <- stoch.lambda(matrices = D2A7_matrices)

##D3A6 
D3A6_stoch <- stoch.lambda(matrices = D3A6_matrices)

##D4A5 
D4A5_stoch <- stoch.lambda(matrices = D4A5_matrices)

##D5A4
D5A4_stoch <- stoch.lambda(matrices = D5A4_matrices)


##DURATION 20 ----

##D1A18 
D1A18_stoch <- stoch.lambda(matrices = D1A18_matrices)

##D2A17 
D2A17_stoch <- stoch.lambda(matrices = D2A17_matrices)

##D3A16 
D3A16_stoch <- stoch.lambda(matrices = D3A16_matrices)

##D4A15 
D4A15_stoch <- stoch.lambda(matrices = D4A15_matrices)

##D5A14 
D5A14_stoch <- stoch.lambda(matrices = D5A14_matrices)



##DURATION 30 ----

##D1A28 
D1A28_stoch <- stoch.lambda(matrices = D1A28_matrices)

##D2A27 
D2A27_stoch <- stoch.lambda(matrices = D2A27_matrices)

##D3A26 
D3A26_stoch <- stoch.lambda(matrices = D3A26_matrices)

##D4A25 
D4A25_stoch <- stoch.lambda(matrices = D4A25_matrices)

##D5A24 
D5A24_stoch <- stoch.lambda(matrices = D5A24_matrices)



#lambda dataframe
lambda_df <- data.frame(lifespan = c(rep(10, 5000), rep(20, 5000), rep(30, 5000)), imm_dur = c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000)), lambda = c(D1A8_stoch$lambda, D2A7_stoch$lambda, D3A6_stoch$lambda, D4A5_stoch$lambda, D5A4_stoch$lambda, D1A18_stoch$lambda, D2A17_stoch$lambda, D3A16_stoch$lambda, D4A15_stoch$lambda, D5A14_stoch$lambda, D1A28_stoch$lambda, D2A27_stoch$lambda, D3A26_stoch$lambda, D4A25_stoch$lambda, D5A24_stoch$lambda))


lambda_df$lifespan <- as.factor(lambda_df$lifespan)
lambda_df$imm_dur <- as.factor(lambda_df$imm_dur)
lambda_df$lambda <- as.numeric(lambda_df$lambda)
is.factor(lambda_df$lifespan)
is.factor(lambda_df$imm_dur)
is.numeric(lambda_df$lambda)

lambda_df$lifespan <- factor(lambda_df$lifespan, levels = c("10", "20", "30"), labels = c("10 year lifespan", "20 year lifespan", "30 year lifespan"))


#plot 
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


#model
mod_lambda <- lm(lambda ~ imm_dur*lifespan, data = lambda_df)

#ANOVA
Anova(mod_lambda)

summary(mod_lambda)

# planned contrasts

contrast1 <- rbind("20 l vs 30 l" = c(0, -1, 1))

fit.contrast(mod_lambda, 'lifespan', contrast1) #not significant

contrast2 <- rbind("10 l vs 20 l" = c(-1, 1, 0))

fit.contrast(mod_lambda, 'lifespan', contrast2) #significant

contrast3 <- rbind("1 dur vs 2 dur" = c(-1, 1, 0, 0, 0))

fit.contrast(mod_lambda, 'imm_dur', contrast3) #significant

contrast4 <- rbind("2 dur vs 3 dur" = c(0, -1, 1, 0, 0))

fit.contrast(mod_lambda, 'imm_dur', contrast4) #significant

contrast5 <- rbind("3 dur vs 4 dur" = c(0, 0, -1, 1, 0))

fit.contrast(mod_lambda, 'imm_dur', contrast5) #significant

contrast6 <- rbind("4 dur vs 5 dur" = c(0, 0, 0, -1, 1))

fit.contrast(mod_lambda, 'imm_dur', contrast6) #significant

contrast7 <- rbind("avg 20-30 l vs 10 l" = c(-1, 1/2, 1/2))

fit.contrast(mod_lambda, 'lifespan', contrast7) #significant


##Poaching ----

##POACHING 10% ----

##P10D1A18 
P10D1A18_stoch <- stoch.lambda(matrices = P10D1A18_matrices)

##P10D2A17 
P10D2A17_stoch <- stoch.lambda(matrices = P10D2A17_matrices)

##P10D3A16 
P10D3A16_stoch <- stoch.lambda(matrices = P10D3A16_matrices)

##P10D4A15 
P10D4A15_stoch <- stoch.lambda(matrices = P10D4A15_matrices)

##D5A14 
P10D5A14_stoch <- stoch.lambda(matrices = P10D5A14_matrices)

##POACHING 20% ----

##P20D1A18 
P20D1A18_stoch <- stoch.lambda(matrices = P20D1A18_matrices)

##P20D2A17 
P20D2A17_stoch <- stoch.lambda(matrices = P20D2A17_matrices)

##P20D3A16 
P20D3A16_stoch <- stoch.lambda(matrices = P20D3A16_matrices)

##P20D4A15 
P20D4A15_stoch <- stoch.lambda(matrices = P20D4A15_matrices)

##D5A14 
P20D5A14_stoch <- stoch.lambda(matrices = P20D5A14_matrices)

##POACHING 30% ----

##P30D1A18 
P30D1A18_stoch <- stoch.lambda(matrices = P30D1A18_matrices)

##P30D2A17 
P30D2A17_stoch <- stoch.lambda(matrices = P30D2A17_matrices)

##P30D3A16 
P30D3A16_stoch <- stoch.lambda(matrices = P30D3A16_matrices)

##P30D4A15 
P30D4A15_stoch <- stoch.lambda(matrices = P30D4A15_matrices)

##D5A14 
P30D5A14_stoch <- stoch.lambda(matrices = P30D5A14_matrices)

##POACHING 40% ----

##P40D1A18 
P40D1A18_stoch <- stoch.lambda(matrices = P40D1A18_matrices)

##P40D2A17 
P40D2A17_stoch <- stoch.lambda(matrices = P40D2A17_matrices)

##P40D3A16 
P40D3A16_stoch <- stoch.lambda(matrices = P40D3A16_matrices)

##P40D4A15 
P40D4A15_stoch <- stoch.lambda(matrices = P40D4A15_matrices)

##D5A14 
P40D5A14_stoch <- stoch.lambda(matrices = P40D5A14_matrices)

##POACHING 50% ----

##P50D1A18 
P50D1A18_stoch <- stoch.lambda(matrices = P50D1A18_matrices)

##P50D2A17 
P50D2A17_stoch <- stoch.lambda(matrices = P50D2A17_matrices)

##P50D3A16 
P50D3A16_stoch <- stoch.lambda(matrices = P50D3A16_matrices)

##P50D4A15 
P50D4A15_stoch <- stoch.lambda(matrices = P50D4A15_matrices)

##D5A14 
P50D5A14_stoch <- stoch.lambda(matrices = P50D5A14_matrices)



#lambda dataframe

poach_lambda_df <- data.frame(poaching = c(rep("0", 5000), rep("10", 5000), rep("20", 5000), rep("30", 5000), rep("40", 5000), rep("50", 5000)), imm_dur = c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000)), lambda = c(D1A18_stoch$lambda, D2A17_stoch$lambda, D3A16_stoch$lambda, D4A15_stoch$lambda, D5A14_stoch$lambda, P10D1A18_stoch$lambda, P10D2A17_stoch$lambda, P10D3A16_stoch$lambda, P10D4A15_stoch$lambda, P10D5A14_stoch$lambda, P20D1A18_stoch$lambda, P20D2A17_stoch$lambda, P20D3A16_stoch$lambda, P20D4A15_stoch$lambda, P20D5A14_stoch$lambda, P30D1A18_stoch$lambda, P30D2A17_stoch$lambda, P30D3A16_stoch$lambda, P30D4A15_stoch$lambda, P30D5A14_stoch$lambda, P40D1A18_stoch$lambda, P40D2A17_stoch$lambda, P40D3A16_stoch$lambda, P40D4A15_stoch$lambda, P40D5A14_stoch$lambda, P50D1A18_stoch$lambda, P50D2A17_stoch$lambda, P50D3A16_stoch$lambda, P50D4A15_stoch$lambda, P50D5A14_stoch$lambda))


poach_lambda_df$poaching <- as.numeric(poach_lambda_df$poaching)
poach_lambda_df$imm_dur <- factor(poach_lambda_df$imm_dur, levels = c("1", "2", "3", "4", "5"), labels = c("1 year", "2 years", "3 years", "4 years", "5 years"))
poach_lambda_df$lambda <- as.numeric(poach_lambda_df$lambda)
is.numeric(poach_lambda_df$poaching)
is.factor(poach_lambda_df$imm_dur)
is.numeric(poach_lambda_df$lambda)


#plot 
ggplot(poach_lambda_df, aes(x = poaching, y = lambda, col = imm_dur)) +
  geom_point()+ geom_smooth(method = lm, se = TRUE)

#model
mod_poach_lambda <- lm(lambda ~ poaching*imm_dur, data = poach_lambda_df)

autoplot(mod_poach_lambda)


#ANCOVA
anova(mod_poach_lambda)

summary(mod_poach_lambda)


