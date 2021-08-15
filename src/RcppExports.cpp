// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ecdf_cpp
NumericVector ecdf_cpp(NumericVector reference, NumericVector sample);
RcppExport SEXP _pidProjCode_ecdf_cpp(SEXP referenceSEXP, SEXP sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type reference(referenceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sample(sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_cpp(reference, sample));
    return rcpp_result_gen;
END_RCPP
}
// bivariate_ecdf_cpp
NumericVector bivariate_ecdf_cpp(NumericVector u_ref, NumericVector v_ref);
RcppExport SEXP _pidProjCode_bivariate_ecdf_cpp(SEXP u_refSEXP, SEXP v_refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u_ref(u_refSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v_ref(v_refSEXP);
    rcpp_result_gen = Rcpp::wrap(bivariate_ecdf_cpp(u_ref, v_ref));
    return rcpp_result_gen;
END_RCPP
}
// bivariate_ecdf_par_cpp
NumericVector bivariate_ecdf_par_cpp(NumericVector u_ref, NumericVector v_ref);
RcppExport SEXP _pidProjCode_bivariate_ecdf_par_cpp(SEXP u_refSEXP, SEXP v_refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u_ref(u_refSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v_ref(v_refSEXP);
    rcpp_result_gen = Rcpp::wrap(bivariate_ecdf_par_cpp(u_ref, v_ref));
    return rcpp_result_gen;
END_RCPP
}
// bivariate_ecdf_lw_cpp
NumericVector bivariate_ecdf_lw_cpp(NumericVector u_ref, NumericVector v_ref);
RcppExport SEXP _pidProjCode_bivariate_ecdf_lw_cpp(SEXP u_refSEXP, SEXP v_refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u_ref(u_refSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v_ref(v_refSEXP);
    rcpp_result_gen = Rcpp::wrap(bivariate_ecdf_lw_cpp(u_ref, v_ref));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pidProjCode_ecdf_cpp", (DL_FUNC) &_pidProjCode_ecdf_cpp, 2},
    {"_pidProjCode_bivariate_ecdf_cpp", (DL_FUNC) &_pidProjCode_bivariate_ecdf_cpp, 2},
    {"_pidProjCode_bivariate_ecdf_par_cpp", (DL_FUNC) &_pidProjCode_bivariate_ecdf_par_cpp, 2},
    {"_pidProjCode_bivariate_ecdf_lw_cpp", (DL_FUNC) &_pidProjCode_bivariate_ecdf_lw_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pidProjCode(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
