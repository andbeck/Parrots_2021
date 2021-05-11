#' Create a stochastic population projection with poaching events
#' 
#' Creates a stochastic projection of population growth from two lists of 
#' simulated matrices
#' 
#' @param matrices List of matrices to sample from.
#' @param poached List of matrices to sample from in poaching years.
#' @param n Vector of initial population sizes.
#' @param time Number of time steps to project population for (Default: `50`).
#' @param interval Time step for each poaching event (Default: `10`).
#' 
#' @return  List `pop.project` representing time series of projected total 
#' population sizes.

poach_function <- function (matrices, poached, n, time = 50, interval = 10) { 
  t <- time
  ## Sample from the list of matrices for each time step
  A <- matrix()

  for (i in 1:t){
    if (i %% interval) {
      A[i] <- sample(matrices, 1, replace = TRUE) 
    } else 
      { A[i] <- sample(poached, 1, replace = TRUE)
    }
  }
   
  ## Matrix multiplication for each randomly sampled matrix and the population 
  ## size vector
  
  pop <- numeric(t)
  for (i in 1:t) {
    n <- A[[i]] %*% n
    
    ## Sum the number of individuals in each stage to give total population size
    pop[i] <- sum(n)
  }
  
  ## Creates a list of population sizes by year
  pop.project <- list(pop.sizes = pop) 
  return(pop.project)
}

