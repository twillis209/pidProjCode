#' @title Genome-wide pairwise association signal sharing test statistic
#'
#' @description Implementation of GPS test statistic from Li et al. (2015), 'Meta-analysis of shared genetic architecture across ten pediatric autoimmune diseases'
#' 
#' @param u Vector of p-values
#' @param v Vector of p-values
#'
#' @return Numeric value of statistic
#' @export
gps <- function(u, v) {
  if(length(u) != length(v)) {
    stop("Lengths of u and v differ")
  }

  n <- length(u)

  cdf_u <- ecdf_cpp(u, u)
  cdf_v <- ecdf_cpp(v, v)
  cdf_u_v <- bivariate_ecdf_par_cpp(u, v)

  denom <- sqrt(cdf_u*cdf_v - (cdf_u^2)*(cdf_v^2))

  # TODO I believe this only happens when max(u) and max(v) have same index
  if(any(denom == 0)) {
    stop("One or more zero values in the denominator")
  }

  max(sqrt(n/log(n))*abs(cdf_u_v - (cdf_u*cdf_v))/denom)
}

#' @title Vector of GPS statistic maximands
#'
#' @param u Vector of p-values
#' @param v Vector of p-values
#'
#' @return Vector of maximands
#' @export
gps_maximands<- function(u, v) {
  if(length(u) != length(v)) {
    stop("Lengths of u and v differ")
  }

  n <- length(u)

  cdf_u <- ecdf_cpp(u, u)
  cdf_v <- ecdf_cpp(v, v)
  cdf_u_v <- bivariate_ecdf_par_cpp(u, v)

  denom <- sqrt(cdf_u*cdf_v - (cdf_u^2)*(cdf_v^2))

  if(any(denom == 0)) {
    stop("One or more zero values in the denominator")
  }

  sqrt(n/log(n))*abs(cdf_u_v - (cdf_u*cdf_v))/denom
}

#' @title Generates samples from two-part mixture of exponentials
#'
#' @param n Sample size
#' @param rates Rate parameter for each exponential component
#' @param alt_weight Mixing weight for first exponential component; second has weight 1-alt_weight
#' @param pval_scale Flag to return sample on [0,1] scale
#'
#' @return Sample from the mixture distribution
#' @export
mix_rexp <- function(n, rates = c(5,1), alt_weight = 0.01, pval_scale = F) {
  n1 <- sum(runif(n) <= alt_weight)

  if(pval_scale) {
    sam <- exp(-c(rexp(n1, rates[1]), rexp(n-n1, rates[2])))
  } else {
    sam <- c(rexp(n1, rates[1]), rexp(n-n1, rates[2]))
  }

  sample(sam, n)
}

#' @title Generates samples from the null distribution of the GPS statistic using samples from the exponential mixture
#'
#' @param n Sample size
#' @param no_snps Number of SNPs/p-values to generate for each simulated GWAS data set
#' @param rates Rate parameter for each exponential component
#' @param alt_weight Mixing weight for first exponential component; second has weight 1-alt_weight
#'
#' @return Sample from the GPS null distribution
#' @export
rgps <- function(n, no_snps, rates = c(5,1), alt_weight = 0.01) {
  gps_results <- numeric(n)

  for(i in 0:(n-1)) {
    j <- 1

    gps_attempt <- NA

    while(j < 6 & is.na(gps_attempt)) {
    sam <- mix_rexp(2*no_snps, rates = rates, alt_weight = alt_weight, pval_scale = T)

    gps_attempt <- try(gps_test_stat(sam[1:no_snps], sam[(no_snps+1):(2*no_snps)]), silent = T)

    j <- j+1
    }

    if(!is.na(gps_attempt)) {
      gps_results[i+1] <- gps_attempt
    } else {
      stop(sprintf("Failed to generate GPS sample realisation %d after 5 attempts", i+1))
    }
  }

  gps_results
}
