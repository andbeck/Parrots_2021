#' Create a stochastic population projection
#' 
#' Creates a stochastic projection of population growth from a list of simulated
#'  matrices, with the option for simple density dependence (based on a carrying
#'  capacity) of adults and/or whole population. 
#' 
#' @param matrices List of matrices to sample from.
#' @param n Vector of initial population sizes.
#' @param nmax Maximum population size (Default: `NULL`).
#' @param smax Maximum number of adults (Default: `NULL`).
#' @param time Number of time steps to project population for (Default: `50`).
#' 
#' @return  List `pop.project` representing time series of projected total 
#' population sizes.

stochastic_proj <- function (matrices, n, nmax = NULL, smax = NULL, time = 50) { 
  t <- time
  ## Sample from the list of matrices for each time step
  A <- matrix()
  for (i in 1:t) {
    A[i] <- sample(matrices, 1, replace = TRUE) 
  }
  
  ## Matrix multiplication for each randomly sampled matrix and the population 
  ## size vector
  pop <- numeric(t)
  for (i in 1:t) {
    n <- A[[i]] %*% n 
    
  ## If the number of adults exceeds the specified maximum, adds the "extra" 
  ## adults back into the juvenile (non-reproductive) stage and applies density
  ## dependence on number of adults (the reproductive stage) 
    
    if (!is.null(smax)) { 
      if (n[3] > smax) {
        n[2] <- n[2] + (n[3]- (n[3] * (smax/n[3]))) 
      n[3] <- n[3] * (smax/n[3]) 
    }}
    
  ## Density dependence for the whole population, if specified
    if (!is.null(nmax)) { 
      if (sum(n) > nmax)
        n <- n * (nmax/sum(n)) 
    }
    
  ## Sum the number of individuals in each stage to give total population size
    pop[i] <- sum(n)
  }
  
 ## Creates a list of population sizes by year
  pop.project <- list(pop.sizes = pop) 
  return(pop.project)
}

