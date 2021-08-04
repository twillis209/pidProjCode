#' @title Genome-wide pairwise association signal sharing test statistic
#'
#' @description Implementation of GPS test statistic from Li et al. (2015), 'Meta-analysis of shared genetic architecture across ten pediatric autoimmune diseases'
#' 
#' @param u Vector of p-values
#' @param v Vector of p-values
#'
#' @return Numeric value of statistic
#' @export
gps_test_stat<- function(u, v) {
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

  max(sqrt(n/log(n))*abs(cdf_u_v - (cdf_u*cdf_v))/denom)
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
  sam <- mix_rexp(n*2*no_snps, rates = rates, alt_weight = alt_weight, pval_scale = T)

  gps_results <- numeric(n)
  # [(1, no_snps, no_snps+1, 2*no_snps)], [(1+2*no_snps, 3*no_snps, 3*no_snps+1, 4*no_snps)]
  for(i in 0:(n-1)) {
    #    print("Indices")
    #    indices <- c((2*i*no_snps+1), ((2*i+1)*no_snps), ((2*i+1)*no_snps+1), ((2*i+2)*no_snps))
    #    print(indices)
    #    
    #    print("Vector lengths")
    #    lengths <- c(length(sam[(2*i*no_snps+1):((2*i+1)*no_snps)]), length(sam[((2*i+1)*no_snps+1):((2*i+2)*no_snps)]))
    #    print(lengths)
    
    gps_results[i+1] <- pidProjCode::gps_test_stat(sam[(2*i*no_snps+1):((2*i+1)*no_snps)], sam[((2*i+1)*no_snps+1):((2*i+2)*no_snps)])
  }

  gps_results
}
