// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mvndns
Eigen::VectorXd mvndns(const Eigen::VectorXd& a, const Eigen::VectorXd& b, const Eigen::MatrixXi& NN, const Eigen::MatrixXd& muCoeff, const Eigen::VectorXd& condSd, const Eigen::VectorXd& beta, int NLevel1, int NLevel2);
RcppExport SEXP _VeccTMVN_mvndns(SEXP aSEXP, SEXP bSEXP, SEXP NNSEXP, SEXP muCoeffSEXP, SEXP condSdSEXP, SEXP betaSEXP, SEXP NLevel1SEXP, SEXP NLevel2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXi& >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type muCoeff(muCoeffSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type condSd(condSdSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type NLevel1(NLevel1SEXP);
    Rcpp::traits::input_parameter< int >::type NLevel2(NLevel2SEXP);
    rcpp_result_gen = Rcpp::wrap(mvndns(a, b, NN, muCoeff, condSd, beta, NLevel1, NLevel2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VeccTMVN_mvndns", (DL_FUNC) &_VeccTMVN_mvndns, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_VeccTMVN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
