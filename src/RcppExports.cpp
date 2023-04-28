// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BL
Rcpp::List BL(arma::vec x, arma::vec y, arma:: mat e, arma:: mat c, arma:: mat w, int maxSteps, int n, double hatBeta, arma:: vec hatEta, arma::vec hatAlpha, arma::vec hatb, double hatInvTauSq1, arma:: vec hatInvTauSq2, arma::mat invSigAlpha0, arma:: mat invSigb0, double hatLambdaSqStar1, double hatLambdaSqStar2, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, int progress);
RcppExport SEXP _marble_BL(SEXP xSEXP, SEXP ySEXP, SEXP eSEXP, SEXP cSEXP, SEXP wSEXP, SEXP maxStepsSEXP, SEXP nSEXP, SEXP hatBetaSEXP, SEXP hatEtaSEXP, SEXP hatAlphaSEXP, SEXP hatbSEXP, SEXP hatInvTauSq1SEXP, SEXP hatInvTauSq2SEXP, SEXP invSigAlpha0SEXP, SEXP invSigb0SEXP, SEXP hatLambdaSqStar1SEXP, SEXP hatLambdaSqStar2SEXP, SEXP hatSigmaSqSEXP, SEXP aStarSEXP, SEXP bStarSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type maxSteps(maxStepsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type hatBeta(hatBetaSEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type hatEta(hatEtaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatAlpha(hatAlphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatb(hatbSEXP);
    Rcpp::traits::input_parameter< double >::type hatInvTauSq1(hatInvTauSq1SEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type hatInvTauSq2(hatInvTauSq2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigAlpha0(invSigAlpha0SEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type invSigb0(invSigb0SEXP);
    Rcpp::traits::input_parameter< double >::type hatLambdaSqStar1(hatLambdaSqStar1SEXP);
    Rcpp::traits::input_parameter< double >::type hatLambdaSqStar2(hatLambdaSqStar2SEXP);
    Rcpp::traits::input_parameter< double >::type hatSigmaSq(hatSigmaSqSEXP);
    Rcpp::traits::input_parameter< double >::type aStar(aStarSEXP);
    Rcpp::traits::input_parameter< double >::type bStar(bStarSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(BL(x, y, e, c, w, maxSteps, n, hatBeta, hatEta, hatAlpha, hatb, hatInvTauSq1, hatInvTauSq2, invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2, hatSigmaSq, aStar, bStar, alpha, gamma, progress));
    return rcpp_result_gen;
END_RCPP
}
// BLSS
Rcpp::List BLSS(arma::vec x, arma::vec y, arma:: mat e, arma:: mat c, arma:: mat w, int maxSteps, int n, double hatBeta, arma:: vec hatEta, arma::vec hatAlpha, arma::vec hatb, double hatInvTauSq1, arma:: vec hatInvTauSq2, double sg1, arma::vec sg2, double hatPiEta, double hatPiBeta, arma::mat invSigAlpha0, arma:: mat invSigb0, double hatLambdaSqStar1, double hatLambdaSqStar2, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, double mu0, double nu0, int progress);
RcppExport SEXP _marble_BLSS(SEXP xSEXP, SEXP ySEXP, SEXP eSEXP, SEXP cSEXP, SEXP wSEXP, SEXP maxStepsSEXP, SEXP nSEXP, SEXP hatBetaSEXP, SEXP hatEtaSEXP, SEXP hatAlphaSEXP, SEXP hatbSEXP, SEXP hatInvTauSq1SEXP, SEXP hatInvTauSq2SEXP, SEXP sg1SEXP, SEXP sg2SEXP, SEXP hatPiEtaSEXP, SEXP hatPiBetaSEXP, SEXP invSigAlpha0SEXP, SEXP invSigb0SEXP, SEXP hatLambdaSqStar1SEXP, SEXP hatLambdaSqStar2SEXP, SEXP hatSigmaSqSEXP, SEXP aStarSEXP, SEXP bStarSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP mu0SEXP, SEXP nu0SEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type maxSteps(maxStepsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type hatBeta(hatBetaSEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type hatEta(hatEtaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatAlpha(hatAlphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatb(hatbSEXP);
    Rcpp::traits::input_parameter< double >::type hatInvTauSq1(hatInvTauSq1SEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type hatInvTauSq2(hatInvTauSq2SEXP);
    Rcpp::traits::input_parameter< double >::type sg1(sg1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sg2(sg2SEXP);
    Rcpp::traits::input_parameter< double >::type hatPiEta(hatPiEtaSEXP);
    Rcpp::traits::input_parameter< double >::type hatPiBeta(hatPiBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigAlpha0(invSigAlpha0SEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type invSigb0(invSigb0SEXP);
    Rcpp::traits::input_parameter< double >::type hatLambdaSqStar1(hatLambdaSqStar1SEXP);
    Rcpp::traits::input_parameter< double >::type hatLambdaSqStar2(hatLambdaSqStar2SEXP);
    Rcpp::traits::input_parameter< double >::type hatSigmaSq(hatSigmaSqSEXP);
    Rcpp::traits::input_parameter< double >::type aStar(aStarSEXP);
    Rcpp::traits::input_parameter< double >::type bStar(bStarSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< int >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(BLSS(x, y, e, c, w, maxSteps, n, hatBeta, hatEta, hatAlpha, hatb, hatInvTauSq1, hatInvTauSq2, sg1, sg2, hatPiEta, hatPiBeta, invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2, hatSigmaSq, aStar, bStar, alpha, gamma, mu0, nu0, progress));
    return rcpp_result_gen;
END_RCPP
}
// RBL
Rcpp::List RBL(arma::vec x, arma::vec y, arma:: mat w, arma:: mat c, arma::mat e, int maxSteps, int n, arma::vec hatAlpha, arma::vec hatb, double hatBeta, arma::vec hatEta, double hatTau, arma::vec hatV, double hatSg1, arma:: vec hatSg2, arma::mat invSigAlpha0, arma:: mat invSigb0, double hatEtaSq1, double hatEtaSq2, double theta, double r1, double r, double a, double b, int progress);
RcppExport SEXP _marble_RBL(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP cSEXP, SEXP eSEXP, SEXP maxStepsSEXP, SEXP nSEXP, SEXP hatAlphaSEXP, SEXP hatbSEXP, SEXP hatBetaSEXP, SEXP hatEtaSEXP, SEXP hatTauSEXP, SEXP hatVSEXP, SEXP hatSg1SEXP, SEXP hatSg2SEXP, SEXP invSigAlpha0SEXP, SEXP invSigb0SEXP, SEXP hatEtaSq1SEXP, SEXP hatEtaSq2SEXP, SEXP thetaSEXP, SEXP r1SEXP, SEXP rSEXP, SEXP aSEXP, SEXP bSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type maxSteps(maxStepsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatAlpha(hatAlphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatb(hatbSEXP);
    Rcpp::traits::input_parameter< double >::type hatBeta(hatBetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatEta(hatEtaSEXP);
    Rcpp::traits::input_parameter< double >::type hatTau(hatTauSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatV(hatVSEXP);
    Rcpp::traits::input_parameter< double >::type hatSg1(hatSg1SEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type hatSg2(hatSg2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigAlpha0(invSigAlpha0SEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type invSigb0(invSigb0SEXP);
    Rcpp::traits::input_parameter< double >::type hatEtaSq1(hatEtaSq1SEXP);
    Rcpp::traits::input_parameter< double >::type hatEtaSq2(hatEtaSq2SEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(RBL(x, y, w, c, e, maxSteps, n, hatAlpha, hatb, hatBeta, hatEta, hatTau, hatV, hatSg1, hatSg2, invSigAlpha0, invSigb0, hatEtaSq1, hatEtaSq2, theta, r1, r, a, b, progress));
    return rcpp_result_gen;
END_RCPP
}
// RBLSS
Rcpp::List RBLSS(arma::vec x, arma::vec y, arma:: mat w, arma:: mat c, arma::mat e, int maxSteps, int n, arma::vec hatAlpha, arma::vec hatb, double hatBeta, arma::vec hatEta, double hatTau, arma::vec hatV, double hatSg1, arma:: vec hatSg2, double sg1, arma:: vec sg2, arma::mat invSigAlpha0, arma:: mat invSigb0, double hatEtaSq1, double hatEtaSq2, double theta, double r1, double r, double a, double b, double hatPiBeta, double hatPiEta, double sh0, double sh1, int progress);
RcppExport SEXP _marble_RBLSS(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP cSEXP, SEXP eSEXP, SEXP maxStepsSEXP, SEXP nSEXP, SEXP hatAlphaSEXP, SEXP hatbSEXP, SEXP hatBetaSEXP, SEXP hatEtaSEXP, SEXP hatTauSEXP, SEXP hatVSEXP, SEXP hatSg1SEXP, SEXP hatSg2SEXP, SEXP sg1SEXP, SEXP sg2SEXP, SEXP invSigAlpha0SEXP, SEXP invSigb0SEXP, SEXP hatEtaSq1SEXP, SEXP hatEtaSq2SEXP, SEXP thetaSEXP, SEXP r1SEXP, SEXP rSEXP, SEXP aSEXP, SEXP bSEXP, SEXP hatPiBetaSEXP, SEXP hatPiEtaSEXP, SEXP sh0SEXP, SEXP sh1SEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type maxSteps(maxStepsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatAlpha(hatAlphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatb(hatbSEXP);
    Rcpp::traits::input_parameter< double >::type hatBeta(hatBetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatEta(hatEtaSEXP);
    Rcpp::traits::input_parameter< double >::type hatTau(hatTauSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatV(hatVSEXP);
    Rcpp::traits::input_parameter< double >::type hatSg1(hatSg1SEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type hatSg2(hatSg2SEXP);
    Rcpp::traits::input_parameter< double >::type sg1(sg1SEXP);
    Rcpp::traits::input_parameter< arma:: vec >::type sg2(sg2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invSigAlpha0(invSigAlpha0SEXP);
    Rcpp::traits::input_parameter< arma:: mat >::type invSigb0(invSigb0SEXP);
    Rcpp::traits::input_parameter< double >::type hatEtaSq1(hatEtaSq1SEXP);
    Rcpp::traits::input_parameter< double >::type hatEtaSq2(hatEtaSq2SEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type hatPiBeta(hatPiBetaSEXP);
    Rcpp::traits::input_parameter< double >::type hatPiEta(hatPiEtaSEXP);
    Rcpp::traits::input_parameter< double >::type sh0(sh0SEXP);
    Rcpp::traits::input_parameter< double >::type sh1(sh1SEXP);
    Rcpp::traits::input_parameter< int >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(RBLSS(x, y, w, c, e, maxSteps, n, hatAlpha, hatb, hatBeta, hatEta, hatTau, hatV, hatSg1, hatSg2, sg1, sg2, invSigAlpha0, invSigb0, hatEtaSq1, hatEtaSq2, theta, r1, r, a, b, hatPiBeta, hatPiEta, sh0, sh1, progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_marble_BL", (DL_FUNC) &_marble_BL, 23},
    {"_marble_BLSS", (DL_FUNC) &_marble_BLSS, 29},
    {"_marble_RBL", (DL_FUNC) &_marble_RBL, 25},
    {"_marble_RBLSS", (DL_FUNC) &_marble_RBLSS, 31},
    {NULL, NULL, 0}
};

RcppExport void R_init_marble(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
