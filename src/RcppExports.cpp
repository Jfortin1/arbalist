// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// initialize_from_hdf5
SEXP initialize_from_hdf5(std::string file, std::string name, size_t nrow, size_t ncol, bool byrow, bool forced);
RcppExport SEXP _arbalist_initialize_from_hdf5(SEXP fileSEXP, SEXP nameSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP byrowSEXP, SEXP forcedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type name(nameSEXP);
    Rcpp::traits::input_parameter< size_t >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< size_t >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    Rcpp::traits::input_parameter< bool >::type forced(forcedSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_from_hdf5(file, name, nrow, ncol, byrow, forced));
    return rcpp_result_gen;
END_RCPP
}
// initialize_from_memory
SEXP initialize_from_memory(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol, bool byrow);
RcppExport SEXP _arbalist_initialize_from_memory(SEXP xSEXP, SEXP iSEXP, SEXP pSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_from_memory(x, i, p, nrow, ncol, byrow));
    return rcpp_result_gen;
END_RCPP
}
// irlba_realized
Rcpp::List irlba_realized(SEXP input, int rank, int nthreads, int seed);
RcppExport SEXP _arbalist_irlba_realized(SEXP inputSEXP, SEXP rankSEXP, SEXP nthreadsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< int >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(irlba_realized(input, rank, nthreads, seed));
    return rcpp_result_gen;
END_RCPP
}
// irlba_tatami
Rcpp::List irlba_tatami(SEXP input, int rank, int nthreads, int seed);
RcppExport SEXP _arbalist_irlba_tatami(SEXP inputSEXP, SEXP rankSEXP, SEXP nthreadsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< int >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(irlba_tatami(input, rank, nthreads, seed));
    return rcpp_result_gen;
END_RCPP
}
// lsi_matrix_stats
Rcpp::List lsi_matrix_stats(SEXP mat, int nthreads);
RcppExport SEXP _arbalist_lsi_matrix_stats(SEXP matSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(lsi_matrix_stats(mat, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// tatami_dim
Rcpp::IntegerVector tatami_dim(SEXP input);
RcppExport SEXP _arbalist_tatami_dim(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_dim(input));
    return rcpp_result_gen;
END_RCPP
}
// tatami_column
Rcpp::NumericVector tatami_column(SEXP input, int i);
RcppExport SEXP _arbalist_tatami_column(SEXP inputSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_column(input, i));
    return rcpp_result_gen;
END_RCPP
}
// tatami_row
Rcpp::NumericVector tatami_row(SEXP input, int i);
RcppExport SEXP _arbalist_tatami_row(SEXP inputSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_row(input, i));
    return rcpp_result_gen;
END_RCPP
}
// apply_subset
SEXP apply_subset(SEXP input, Rcpp::IntegerVector subset, bool row);
RcppExport SEXP _arbalist_apply_subset(SEXP inputSEXP, SEXP subsetSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_subset(input, subset, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_transpose
SEXP apply_transpose(SEXP input);
RcppExport SEXP _arbalist_apply_transpose(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_transpose(input));
    return rcpp_result_gen;
END_RCPP
}
// apply_bind
SEXP apply_bind(Rcpp::List input, bool row);
RcppExport SEXP _arbalist_apply_bind(SEXP inputSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type input(inputSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_bind(input, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_addition
SEXP apply_addition(SEXP input, Rcpp::NumericVector val, bool row);
RcppExport SEXP _arbalist_apply_addition(SEXP inputSEXP, SEXP valSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_addition(input, val, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_multiplication
SEXP apply_multiplication(SEXP input, Rcpp::NumericVector val, bool row);
RcppExport SEXP _arbalist_apply_multiplication(SEXP inputSEXP, SEXP valSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_multiplication(input, val, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_subtraction
SEXP apply_subtraction(SEXP input, Rcpp::NumericVector val, bool right, bool row);
RcppExport SEXP _arbalist_apply_subtraction(SEXP inputSEXP, SEXP valSEXP, SEXP rightSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type right(rightSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_subtraction(input, val, right, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_division
SEXP apply_division(SEXP input, Rcpp::NumericVector val, bool right, bool row);
RcppExport SEXP _arbalist_apply_division(SEXP inputSEXP, SEXP valSEXP, SEXP rightSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type right(rightSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_division(input, val, right, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_log
SEXP apply_log(SEXP input, double base);
RcppExport SEXP _arbalist_apply_log(SEXP inputSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< double >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_log(input, base));
    return rcpp_result_gen;
END_RCPP
}
// apply_log1p
SEXP apply_log1p(SEXP input);
RcppExport SEXP _arbalist_apply_log1p(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_log1p(input));
    return rcpp_result_gen;
END_RCPP
}
// apply_abs
SEXP apply_abs(SEXP input);
RcppExport SEXP _arbalist_apply_abs(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_abs(input));
    return rcpp_result_gen;
END_RCPP
}
// apply_sqrt
SEXP apply_sqrt(SEXP input);
RcppExport SEXP _arbalist_apply_sqrt(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_sqrt(input));
    return rcpp_result_gen;
END_RCPP
}
// apply_round
SEXP apply_round(SEXP input);
RcppExport SEXP _arbalist_apply_round(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_round(input));
    return rcpp_result_gen;
END_RCPP
}
// apply_exp
SEXP apply_exp(SEXP input);
RcppExport SEXP _arbalist_apply_exp(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_exp(input));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_arbalist_initialize_from_hdf5", (DL_FUNC) &_arbalist_initialize_from_hdf5, 6},
    {"_arbalist_initialize_from_memory", (DL_FUNC) &_arbalist_initialize_from_memory, 6},
    {"_arbalist_irlba_realized", (DL_FUNC) &_arbalist_irlba_realized, 4},
    {"_arbalist_irlba_tatami", (DL_FUNC) &_arbalist_irlba_tatami, 4},
    {"_arbalist_lsi_matrix_stats", (DL_FUNC) &_arbalist_lsi_matrix_stats, 2},
    {"_arbalist_tatami_dim", (DL_FUNC) &_arbalist_tatami_dim, 1},
    {"_arbalist_tatami_column", (DL_FUNC) &_arbalist_tatami_column, 2},
    {"_arbalist_tatami_row", (DL_FUNC) &_arbalist_tatami_row, 2},
    {"_arbalist_apply_subset", (DL_FUNC) &_arbalist_apply_subset, 3},
    {"_arbalist_apply_transpose", (DL_FUNC) &_arbalist_apply_transpose, 1},
    {"_arbalist_apply_bind", (DL_FUNC) &_arbalist_apply_bind, 2},
    {"_arbalist_apply_addition", (DL_FUNC) &_arbalist_apply_addition, 3},
    {"_arbalist_apply_multiplication", (DL_FUNC) &_arbalist_apply_multiplication, 3},
    {"_arbalist_apply_subtraction", (DL_FUNC) &_arbalist_apply_subtraction, 4},
    {"_arbalist_apply_division", (DL_FUNC) &_arbalist_apply_division, 4},
    {"_arbalist_apply_log", (DL_FUNC) &_arbalist_apply_log, 2},
    {"_arbalist_apply_log1p", (DL_FUNC) &_arbalist_apply_log1p, 1},
    {"_arbalist_apply_abs", (DL_FUNC) &_arbalist_apply_abs, 1},
    {"_arbalist_apply_sqrt", (DL_FUNC) &_arbalist_apply_sqrt, 1},
    {"_arbalist_apply_round", (DL_FUNC) &_arbalist_apply_round, 1},
    {"_arbalist_apply_exp", (DL_FUNC) &_arbalist_apply_exp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_arbalist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
