% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ecdf_cpp}
\alias{ecdf_cpp}
\title{Empirical cdf}
\usage{
ecdf_cpp(reference, sample)
}
\arguments{
\item{reference}{NumericVector containing reference points to which to fit empirical cdf}

\item{sample}{NumericVector containing points at which to evaluate empirical cdf}
}
\value{
NumericVector of estimated quantiles of sample values
}
\description{
This function estimates an empirical cdf using the reference parameter and then
evaluates the estimated empirical cdf at the points specified in the sample parameter.
}
\author{
Tom Willis
}
