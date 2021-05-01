#' @param u
#' @param v
#' @importFrom mltools empirical_cdf
gps_test_stat<- function(u, v) {
  if(length(u) != length(v)) {
    stop("Lengths of u and v differ")
  }

  n <- length(u)
  cdf_u <- empirical_cdf(u, u)$CDF
  cdf_v <- empirical_cdf(v, v)$CDF
  cdf_u_v <- empirical_cdf(data.table(u=u, v=v), data.table(u=u, v=v))$CDF

  # TODO need to check for 0 values in the denominator
  max(sqrt(n/log(n))*abs(cdf_u_v - (cdf_u*cdf_v))/sqrt(cdf_u*cdf_v - (cdf_u^2)*(cdf_v^2)))
}
