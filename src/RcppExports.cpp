// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// singleKmeansC
List singleKmeansC(const arma::mat& x, arma::mat& centers);
RcppExport SEXP _packageKmeans_singleKmeansC(SEXP xSEXP, SEXP centersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type centers(centersSEXP);
    rcpp_result_gen = Rcpp::wrap(singleKmeansC(x, centers));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_packageKmeans_singleKmeansC", (DL_FUNC) &_packageKmeans_singleKmeansC, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_packageKmeans(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
