############ Probability distributions ############

#'@useDynLib hahmmr
#'@import Rcpp
NULL

## refer to https://github.com/evanbiederstedt/poilogcpp

#' Returns the density for the Poisson lognormal distribution with parameters mu and sig
#' 
#' @param x vector of integers, the observations
#' @param mu mean of lognormal distribution
#' @param sig standard deviation of lognormal distribution
#' @param log boolean Return the log density if TRUE (default=FALSE)
#' @return NULL
#' @examples
#' p = dpoilog(1, 1, 1)
#' @export
dpoilog <- function(x, mu, sig, log=FALSE){
  if (!(length(x) == length(mu) & length(x) == length(sig))) stop('dpoilog: All parameters must be same length') 
  if (any((x[x!=0]/trunc(x[x!=0]))!=1)) stop('dpoilog: all x must be integers')
  if (any(x<0)) stop('dpoilog: one or several values of x are negative')
  if (!all(is.finite(c(mu,sig)))) {
    stop('dpoilog: all parameters should be finite')
  }
  if (any(is.na(c(x,mu,sig)))) stop('dpoilog: Parameters cannot be NA')
  if (any(sig<=0)) {
      stop(c('dpoilog: sig is not larger than 0', unique(sig)))
  }

  p = poilog1(x=as.integer(x), my=as.double(mu), sig=as.double(sig^2))

  p[p == 0] = 1e-15

  if (log) {
    return(log(p))
  } else {
    return(p)
  }
}

#' calculate joint likelihood of a PLN model
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels
#' @param d numeric Total library size
#' @param mu numeric Global mean expression
#' @param sig numeric Global standard deviation of expression
#' @param phi numeric Fold change of expression
#' @return numeric Joint log likelihood
#' @examples
#' l_lnpois(c(1, 2), c(1, 2), 1, 1, 1)
#' @export
l_lnpois = function(Y_obs, lambda_ref, d, mu, sig, phi = 1) {
    if (any(sig <= 0)) {stop(glue('negative sigma. value: {sig}'))}
    if (length(sig) == 1) {sig = rep(sig, length(Y_obs))}
    sum(log(dpoilog(Y_obs, mu + log(phi * d * lambda_ref), sig)))
}

#' Beta-binomial distribution density function
#' A distribution is beta-binomial if p, the probability of success, 
#' in a binomial distribution has a beta distribution with shape 
#' parameters alpha > 0 and beta > 0
#' For more details, see extraDistr::dbbinom
#'
#' @param x vector of quantiles
#' @param size number of trials (zero or more)
#' @param alpha numeric (default=1)
#' @param beta numeric (default=1)
#' @param log boolean (default=FALSE)
#' @return density values returned as numeric vector
#' @examples
#' dbbinom(1, 1, 1, 1)
#' @export
dbbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
    cppdbbinom(x, size, alpha, beta, log[1L])
}

#' calculate joint likelihood of allele data
#' @param AD numeric vector Variant allele depth
#' @param DP numeric vector Total allele depth
#' @param alpha numeric Alpha parameter of Beta-Binomial distribution
#' @param beta numeric Beta parameter of Beta-Binomial distribution
#' @return numeric Joint log likelihood
#' @examples
#' l_bbinom(c(1, 2), c(1, 2), 1, 1)
#' @export
l_bbinom = function(AD, DP, alpha, beta) {
    sum(dbbinom(AD, DP, alpha, beta, log = TRUE))
}