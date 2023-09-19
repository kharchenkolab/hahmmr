// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppdbbinom
NumericVector cppdbbinom(const NumericVector& x, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& log_prob);
RcppExport SEXP _hahmmr_cppdbbinom(SEXP xSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP log_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_prob(log_probSEXP);
    rcpp_result_gen = Rcpp::wrap(cppdbbinom(x, size, alpha, beta, log_prob));
    return rcpp_result_gen;
END_RCPP
}
// cpp_dgpois
NumericVector cpp_dgpois(const NumericVector& x, const NumericVector& alpha, const NumericVector& beta, const bool& log_prob);
RcppExport SEXP _hahmmr_cpp_dgpois(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP log_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_prob(log_probSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dgpois(x, alpha, beta, log_prob));
    return rcpp_result_gen;
END_RCPP
}
// logSumExp
double logSumExp(const arma::vec& x);
RcppExport SEXP _hahmmr_logSumExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logSumExp(x));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_compute
double likelihood_compute(Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, arma::cube logPi, int n, int m);
RcppExport SEXP _hahmmr_likelihood_compute(SEXP logphiSEXP, SEXP logprobSEXP, SEXP logPiSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type logphi(logphiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type logprob(logprobSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_compute(logphi, logprob, logPi, n, m));
    return rcpp_result_gen;
END_RCPP
}
// forward_backward_compute
Rcpp::NumericMatrix forward_backward_compute(Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, arma::cube logPi, int n, int m);
RcppExport SEXP _hahmmr_forward_backward_compute(SEXP logphiSEXP, SEXP logprobSEXP, SEXP logPiSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type logphi(logphiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type logprob(logprobSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_backward_compute(logphi, logprob, logPi, n, m));
    return rcpp_result_gen;
END_RCPP
}
// viterbi_compute
Rcpp::NumericVector viterbi_compute(Rcpp::NumericVector log_delta, Rcpp::NumericMatrix logprob, arma::cube logPi, int n, int m, Rcpp::NumericMatrix nu, Rcpp::NumericVector z);
RcppExport SEXP _hahmmr_viterbi_compute(SEXP log_deltaSEXP, SEXP logprobSEXP, SEXP logPiSEXP, SEXP nSEXP, SEXP mSEXP, SEXP nuSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type log_delta(log_deltaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type logprob(logprobSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi_compute(log_delta, logprob, logPi, n, m, nu, z));
    return rcpp_result_gen;
END_RCPP
}
// roman2int_internal
int roman2int_internal(Rcpp::StringVector letters, int nchar);
RcppExport SEXP _hahmmr_roman2int_internal(SEXP lettersSEXP, SEXP ncharSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type letters(lettersSEXP);
    Rcpp::traits::input_parameter< int >::type nchar(ncharSEXP);
    rcpp_result_gen = Rcpp::wrap(roman2int_internal(letters, nchar));
    return rcpp_result_gen;
END_RCPP
}
// fit_lnpois_cpp
arma::rowvec fit_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d);
RcppExport SEXP _hahmmr_fit_lnpois_cpp(SEXP Y_obsSEXP, SEXP lambda_refSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type Y_obs(Y_obsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambda_ref(lambda_refSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_lnpois_cpp(Y_obs, lambda_ref, d));
    return rcpp_result_gen;
END_RCPP
}
// poilog1
std::vector<double> poilog1(std::vector<int> x, std::vector<double> my, std::vector<double> sig);
RcppExport SEXP _hahmmr_poilog1(SEXP xSEXP, SEXP mySEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type my(mySEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(poilog1(x, my, sig));
    return rcpp_result_gen;
END_RCPP
}
// l_lnpois_cpp
double l_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig, double phi);
RcppExport SEXP _hahmmr_l_lnpois_cpp(SEXP Y_obsSEXP, SEXP lambda_refSEXP, SEXP dSEXP, SEXP muSEXP, SEXP sigSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type Y_obs(Y_obsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambda_ref(lambda_refSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lnpois_cpp(Y_obs, lambda_ref, d, mu, sig, phi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hahmmr_cppdbbinom", (DL_FUNC) &_hahmmr_cppdbbinom, 5},
    {"_hahmmr_cpp_dgpois", (DL_FUNC) &_hahmmr_cpp_dgpois, 4},
    {"_hahmmr_logSumExp", (DL_FUNC) &_hahmmr_logSumExp, 1},
    {"_hahmmr_likelihood_compute", (DL_FUNC) &_hahmmr_likelihood_compute, 5},
    {"_hahmmr_forward_backward_compute", (DL_FUNC) &_hahmmr_forward_backward_compute, 5},
    {"_hahmmr_viterbi_compute", (DL_FUNC) &_hahmmr_viterbi_compute, 7},
    {"_hahmmr_roman2int_internal", (DL_FUNC) &_hahmmr_roman2int_internal, 2},
    {"_hahmmr_fit_lnpois_cpp", (DL_FUNC) &_hahmmr_fit_lnpois_cpp, 3},
    {"_hahmmr_poilog1", (DL_FUNC) &_hahmmr_poilog1, 3},
    {"_hahmmr_l_lnpois_cpp", (DL_FUNC) &_hahmmr_l_lnpois_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_hahmmr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
