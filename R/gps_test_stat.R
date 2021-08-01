#' @title Genome-wide pairwise association signal statistic
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
