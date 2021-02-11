#' @title Sample from two-Gaussian mixture.
#'
#' @description Draws a sample from a mixture of a standard normal and a non-standard normal with specified mean and standard deviation.
#'
#' @details Mixing weights are specified by first element of parameter vector, say pi, as pi and 1-pi. 
#'
#' @param n Size of sample
#' @param pi Mixing weight for standard normal component
#' @param mean Mean for non-standard normal component
#' @param sd Standard deviation for non-standard normal component
#'
#' @return Sample from mixture process
#' 
#' @importFrom stats runif rnorm
#' @export
#'
#' @examples
#'
#' sam<-rmixture(1e3, 0.1, 2, 3)
#' 
rmixture<-function(n, pi, mean, sd) {
  n0<-sum(runif(n)<=pi)

  c(rnorm(n0), (sd*rnorm(n-n0))+mean)
}
