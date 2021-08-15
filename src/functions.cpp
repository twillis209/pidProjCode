#include <Rcpp.h>
#include <Rmath.h>
#include <RcppParallel.h>
#include <Eigen/Dense>
#include <fastCDFOnSample.h>

using namespace Rcpp;
using namespace RcppParallel;

//' @title Empirical cdf
//'
//' @description This function estimates an empirical cdf using the reference parameter and then 
//' evaluates the estimated empirical cdf at the points specified in the sample parameter. 
//' 
//' @param reference NumericVector containing reference points to which to fit empirical cdf
//' @param sample NumericVector containing points at which to evaluate empirical cdf
//' @return NumericVector of estimated quantiles of sample values
//' 
//' @export
//' @author Tom Willis
// [[Rcpp::export]]
NumericVector ecdf_cpp(NumericVector reference, NumericVector sample) {
  NumericVector sortedRef = clone(reference);

  std::sort(sortedRef.begin(), sortedRef.end());

  NumericVector estimatedQuantiles(sample.size());

  for(int i = 0; i < sample.size(); ++i) {
    estimatedQuantiles[i] = (std::upper_bound(sortedRef.begin(), sortedRef.end(), sample[i])-sortedRef.begin());
  }

  return estimatedQuantiles/((double) sortedRef.size());
}

//' @title Bivariate empirical cdf
//'
//' @description This function fits a bivariate ecdf to the specified pair of variables and then evaluates the ecdf at the same points.
//'
//' @details Assumes u_ref and v_ref specify pairs of points in the correct order.
//' 
//' @param u_ref NumericVector reference points in first dimension
//' @param v_ref NumericVector reference points in second dimension
//'
//' @return NumericVector of estimated quantiles of specified pairs
//' @export
//' @author Tom Willis
// [[Rcpp::export]]
NumericVector bivariate_ecdf_cpp(NumericVector u_ref, NumericVector v_ref) {
  int n = u_ref.size();
  NumericVector emp_quantiles(n);

  // Naive implementation
  for(int i = 0; i < n; ++i) {

    int count = 0;

    for(int j = 0; j < n; ++j) {
      if(u_ref[j] <= u_ref[i] && v_ref[j] <= v_ref[i]) {
        count++;
      }
    }

    emp_quantiles[i] = count;
  }

  return emp_quantiles/((double) n);
}

// [[Rcpp::depends(RcppParallel)]]
struct bivariate_ecdf_worker : public Worker {
  const RVector<double> u_vec;
  const RVector<double> v_vec;

  RVector<double> emp_quantiles;

  bivariate_ecdf_worker(const RVector<double> u_vec,
                        const RVector<double> v_vec,
                        RVector<double> emp_quantiles)
    : u_vec(u_vec), v_vec(v_vec), emp_quantiles(emp_quantiles) {}

  int n = u_vec.length();

  void operator()(std::size_t begin, std::size_t end) {

    for(int i = begin; i < end; ++i) {

      int count = 0;

      for(int j = 0; j < n; ++j) {
        if(u_vec[j] <= u_vec[i] && v_vec[j] <= v_vec[i]) {
          count++;
        }
      }

      emp_quantiles[i] = count;
    }
  }

};

//' @title Parallel implementation of bivariate empirical cdf
//'
//' @description This function fits a bivariate ecdf to the specified pair of variables and then evaluates the ecdf at the same points.
//'
//' @details Assumes u_ref and v_ref specify pairs of points in the correct order.
//' 
//' @param u_ref NumericVector reference points in first dimension
//' @param v_ref NumericVector reference points in second dimension
//'
//' @return NumericVector of estimated quantiles of specified pairs
//' 
//' @export
//' @author Tom Willis
// [[Rcpp::export]]
NumericVector bivariate_ecdf_par_cpp(NumericVector u_ref, NumericVector v_ref) {

  int n = u_ref.size();

  NumericVector emp_quantiles(n);

  const RVector<double> u_ref_safe(u_ref);
  const RVector<double> v_ref_safe(v_ref);
  RVector<double> emp_safe(emp_quantiles);

  bivariate_ecdf_worker ecdf_worker(u_ref_safe, v_ref_safe, emp_safe);

  parallelFor(0, n, ecdf_worker);

  return emp_quantiles/((double) n);
}

//' @title Bivariate empirical cdf using novel divide-and-conquer algorithm
//'
//' 
//' @param u_ref NumericVector reference points in first dimension
//' @param v_ref NumericVector reference points in second dimension
//'
//' @return NumericVector of estimated quantiles of specified pairs
//' 
//' @export
//' @author Tom Willis
// [[Rcpp::export]]
NumericVector bivariate_ecdf_lw_cpp(NumericVector u_ref, NumericVector v_ref) {
  int n = u_ref.size();

  Eigen::ArrayXd toAdd = Eigen::ArrayXd::Constant(n, 1.);

  Eigen::ArrayXXd ptr(2, n);

  for(int i = 0; i < n; i++) {
    ptr(0, i) = u_ref[i];
    ptr(1, i) = v_ref[i];
  }

  // TODO surely it's not necessary to copy about like this?
  Eigen::ArrayXd ecdf_arr = StOpt::fastCDFOnSample(ptr, toAdd);

  // TODO could just skip this step?
  //std::vector<double> ecdf_vec(ecdf_arr.data(), ecdf_arr.data() + ecdf_arr.size());

  NumericVector ecdf_nvec(ecdf_arr.data(), ecdf_arr.data() + ecdf_arr.size());
  // ArrayXd does not provide begin() and end() methods
  //  NumericVector ecdf_nvec(ecdf_vec.begin(), ecdf_vec.end());

  return ecdf_nvec;
}
