# Tanith's stochastic population growth function, stochastic_proj. 
# Creates a stochastic population projection from a list of simulated matrices, with the option for density dependence of adults and/or whole population
# Inputs are the list of matrices, the initial population vector, the maximum population size, the maximum number of individuals in the adult stage, and the time to project to.
# Output is a time series of population sizes that can be plotted or used to calculate a cumulative distribution function (CDF) for quasi-extinction probabilities.
#Function can be used in a for loop for life stage simulation analysis (LSA).

stochastic_proj <- function (matrices, n, nmax = NULL, smax = NULL, time = 50) { #matrices are a list of simulated matrices, n is the vector of initial stage sizes, nmax is carrying capacity, smax is maximum adults, time is the number of years to project 
  t <- time
  A <- matrix()
  for (i in 1:t) {
    A[i] <- sample(matrices, 1, replace = TRUE) #samples from the list of matrices for each year
  }
  
  pop <- numeric(t)
  for (i in 1:t) {
    n <- A[[i]] %*% n #matrix multiplication, selecting each randomly sampled matrix and multiplying by the stage vector
    if (!is.null(smax)) { #only applies density dependence if specified
      if (n[3] > smax) #if the number of adults exceeds the specified maximum
        n[2] <- n[2] + (n[3]- (n[3] * (smax/n[3]))) #adds the "extra" adults back into the juvenile (non-reproductive) stage 
      n[3] <- n[3] * (smax/n[3]) #density dependence on number of adults (the reproductive stage) 
    }
    if (!is.null(nmax)) { 
      if (sum(n) > nmax) 
        n <- n * (nmax/sum(n)) #density dependence for whole population, if specified
    }
    pop[i] <- sum(n) #the sum of the number of individuals in each stage to give total population size
  }
  pop.project <- list(pop.sizes = pop) #creates a list of population sizes by year
  pop.project
}

