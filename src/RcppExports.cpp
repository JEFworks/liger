// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// gseaRandCore
Rcpp::List gseaRandCore(arma::vec sset, arma::vec eso, int nsamples, int seed);
RcppExport SEXP liger_gseaRandCore(SEXP ssetSEXP, SEXP esoSEXP, SEXP nsamplesSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type sset(ssetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eso(esoSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    __result = Rcpp::wrap(gseaRandCore(sset, eso, nsamples, seed));
    return __result;
END_RCPP
}
// gseaBulkCore
Rcpp::List gseaBulkCore(arma::mat setm, arma::vec eso, int nsamples, int seed);
RcppExport SEXP liger_gseaBulkCore(SEXP setmSEXP, SEXP esoSEXP, SEXP nsamplesSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type setm(setmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eso(esoSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    __result = Rcpp::wrap(gseaBulkCore(setm, eso, nsamples, seed));
    return __result;
END_RCPP
}
