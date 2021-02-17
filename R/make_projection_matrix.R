#' Create a stage-based population projection matrix
#'
#' Create a stage-based population projection matrix parameterised according to
#' the fixed stage duration method of Caswell, 2001 (section 6.4.2).  Inputs are
#' vectors of survival, fecundity and stage durations.  The calculation assumes
#' a value of lambda, which may be supplied by the user, to derive stage-based
#' rates of maturation. Optionally, the maturation rate estimates can be refined
#' by iterating the procedure to allow lambda to converge to a fixed value
#' (section 6.4.4, Caswell, 2001).
#'
#' We assume a post-breeding census, so the matrix entries for transitions to
#' the newborn stage (i.e. reproduction) are constructed by combining
#' pre-breeder and breeder survival/maturation with fecundity estimates (see
#' Kendall et al., 2019; doi:10.1016/j.ecolmodel.2019.03.011).
#'
#' @param surv Stage-based survival estimates
#' @param fec Stage-based fecundity (number of offspring) estimates
#' @param dur Stage durations. If `NULL` the matrix reduces to an age-structured
#'   projection matrix (Leslie matrix) (Default: `NULL`)
#' @param lambda Estimate of lambda to use in iterations (Default: `1`)
#' @param iterate Whether to calculate maturation iteratively (Default: `FALSE`)
#' @param tol Tolerance for convergence of lambda value in iteration (Default:
#'   `1e-3`)
#'
#' @return Projection matrix `A` representing transitions in a stage-based
#'   lifecycle with post-breeding census.
make_projection_matrix <- function (surv, fec, dur = NULL, lambda = 1,
                                    iterate = FALSE, tol = 1e-3) {
  dimA <- length(surv)
  ## Expect vectors of survival and fecundity to have same length
  if (length(fec) != dimA) {
    stop("Expecting equal length survival and fecundity vectors")
  }
  ## Expect vectors of survival and stage length to have same length (or set
  ## stage length vector if not specified)
  if (is.null(dur)) {
    dur <- rep(1, dimA)
  } else if (length(dur) != dimA) {
    stop("Expecting equal length survival and stage length vectors")
  }

  ## Calculate probability of growth conditional on survival (requires
  ## surv/lambda < 1)
  surv1 <- surv/lambda
  if (any(surv1 >= 1 + tol)) {
    stop("Expecting survival/lambda < 1")
  }
  gamma <- (surv1^dur - surv1^(dur-1))/(surv1^dur - 1) # Caswell 2001, eq 6.103

  # Create matrix of survival and growth transitions (see Caswell 2001, eqs
  # 6.89, 6.97 and 6.98)
  A <- diag(surv*(1 - gamma))  # P_i
  if (is.matrix(A[-1,-dimA])) {
    diag(A[-1,-dimA]) <- (surv*gamma)[-dimA]  # G_i
  } else {
    ## If matrix is 2x2 the subsetted value is returned as a numeric vector
    ## rather than as a matrix subset and we can't use diag()<-, so replace
    ## the value directly instead
    A[-1,-dimA] <- (surv*gamma)[-dimA]  # G_i
  }

  # For post-breeding census, we incorporate survival/maturation of
  # breeders/pre-breeders into the reproductive transition entries
  # N.B. this formulation allows individuals maturing to breeding class
  # to reproduce
  Fmat <- diag(dimA)
  Fmat[1,] <- Fmat[1,] + fec
  A <- Fmat %*% A

  if (iterate) {
    ## If lambda value produced by the model doesn't match the value used for
    ## estimates of maturation rates, iterate to improve the estimates
    lambda_new <- Re(eigen(A, only.values = TRUE)$values[1])
    if (abs(lambda_new - lambda) < tol) {
      return(A)
    } else {
      return(makeProjectionMatrix(surv, fec, dur, lambda_new, iterate = TRUE, tol = tol))
    }
  }
  return(A)
}