// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fastNumericUnique
arma::mat fastNumericUnique(arma::mat x, size_t expected);
RcppExport SEXP _Prola2020_Sci_Adv_fastNumericUnique(SEXP xSEXP, SEXP expectedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type expected(expectedSEXP);
    rcpp_result_gen = Rcpp::wrap(fastNumericUnique(x, expected));
    return rcpp_result_gen;
END_RCPP
}
// decayModelExecute
NumericVector decayModelExecute(int outputLength, arma::vec ii, arma::mat allAlphas, arma::mat allBetas, double timeinterval, arma::vec prompts, arma::vec promptStarts, arma::vec promptEnds, List plans, bool grad, bool multithread);
RcppExport SEXP _Prola2020_Sci_Adv_decayModelExecute(SEXP outputLengthSEXP, SEXP iiSEXP, SEXP allAlphasSEXP, SEXP allBetasSEXP, SEXP timeintervalSEXP, SEXP promptsSEXP, SEXP promptStartsSEXP, SEXP promptEndsSEXP, SEXP plansSEXP, SEXP gradSEXP, SEXP multithreadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type outputLength(outputLengthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type allAlphas(allAlphasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type allBetas(allBetasSEXP);
    Rcpp::traits::input_parameter< double >::type timeinterval(timeintervalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prompts(promptsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type promptStarts(promptStartsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type promptEnds(promptEndsSEXP);
    Rcpp::traits::input_parameter< List >::type plans(plansSEXP);
    Rcpp::traits::input_parameter< bool >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< bool >::type multithread(multithreadSEXP);
    rcpp_result_gen = Rcpp::wrap(decayModelExecute(outputLength, ii, allAlphas, allBetas, timeinterval, prompts, promptStarts, promptEnds, plans, grad, multithread));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Prola2020_Sci_Adv_fastNumericUnique", (DL_FUNC) &_Prola2020_Sci_Adv_fastNumericUnique, 2},
    {"_Prola2020_Sci_Adv_decayModelExecute", (DL_FUNC) &_Prola2020_Sci_Adv_decayModelExecute, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_Prola2020_Sci_Adv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}