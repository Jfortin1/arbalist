// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// count_fragment_size_distributions
Rcpp::IntegerVector count_fragment_size_distributions(std::string fragment_file);
RcppExport SEXP _arbalist_count_fragment_size_distributions(SEXP fragment_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type fragment_file(fragment_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(count_fragment_size_distributions(fragment_file));
    return rcpp_result_gen;
END_RCPP
}
// create_pseudobulk_file
void create_pseudobulk_file(Rcpp::Nullable<Rcpp::CharacterVector> fragment_files, std::string output_file, Rcpp::Nullable<Rcpp::CharacterVector> cellnames);
RcppExport SEXP _arbalist_create_pseudobulk_file(SEXP fragment_filesSEXP, SEXP output_fileSEXP, SEXP cellnamesSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type fragment_files(fragment_filesSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_file(output_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type cellnames(cellnamesSEXP);
    create_pseudobulk_file(fragment_files, output_file, cellnames);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_arbalist_count_fragment_size_distributions", (DL_FUNC) &_arbalist_count_fragment_size_distributions, 1},
    {"_arbalist_create_pseudobulk_file", (DL_FUNC) &_arbalist_create_pseudobulk_file, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_arbalist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
