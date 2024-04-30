// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_gene_to_geneset_scores
Rcpp::NumericMatrix rcpp_gene_to_geneset_scores(int n_gs, Rcpp::IntegerVector gs_index, Rcpp::IntegerVector gs_geneindex, Rcpp::NumericMatrix gene_score);
RcppExport SEXP _goat_rcpp_gene_to_geneset_scores(SEXP n_gsSEXP, SEXP gs_indexSEXP, SEXP gs_geneindexSEXP, SEXP gene_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_gs(n_gsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type gs_index(gs_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type gs_geneindex(gs_geneindexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type gene_score(gene_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_gene_to_geneset_scores(n_gs, gs_index, gs_geneindex, gene_score));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_geneset_null
Rcpp::NumericVector rcpp_geneset_null(Rcpp::NumericVector gene_scores, Rcpp::IntegerVector geneset_sizes, int max_geneset_size, int niter);
RcppExport SEXP _goat_rcpp_geneset_null(SEXP gene_scoresSEXP, SEXP geneset_sizesSEXP, SEXP max_geneset_sizeSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gene_scores(gene_scoresSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type geneset_sizes(geneset_sizesSEXP);
    Rcpp::traits::input_parameter< int >::type max_geneset_size(max_geneset_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_geneset_null(gene_scores, geneset_sizes, max_geneset_size, niter));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_null_distributions
Rcpp::List rcpp_null_distributions(Rcpp::NumericVector gene_scores, Rcpp::IntegerVector geneset_sizes, int max_geneset_size, int niter);
RcppExport SEXP _goat_rcpp_null_distributions(SEXP gene_scoresSEXP, SEXP geneset_sizesSEXP, SEXP max_geneset_sizeSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gene_scores(gene_scoresSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type geneset_sizes(geneset_sizesSEXP);
    Rcpp::traits::input_parameter< int >::type max_geneset_size(max_geneset_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_null_distributions(gene_scores, geneset_sizes, max_geneset_size, niter));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_dsnorm_logsum
double rcpp_dsnorm_logsum(Rcpp::NumericVector x_min_mean, int N, double xi, double sd);
RcppExport SEXP _goat_rcpp_dsnorm_logsum(SEXP x_min_meanSEXP, SEXP NSEXP, SEXP xiSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x_min_mean(x_min_meanSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_dsnorm_logsum(x_min_mean, N, xi, sd));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_goat_rcpp_gene_to_geneset_scores", (DL_FUNC) &_goat_rcpp_gene_to_geneset_scores, 4},
    {"_goat_rcpp_geneset_null", (DL_FUNC) &_goat_rcpp_geneset_null, 4},
    {"_goat_rcpp_null_distributions", (DL_FUNC) &_goat_rcpp_null_distributions, 4},
    {"_goat_rcpp_dsnorm_logsum", (DL_FUNC) &_goat_rcpp_dsnorm_logsum, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_goat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}