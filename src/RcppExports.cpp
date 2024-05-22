// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// aggregate_counts
Rcpp::NumericMatrix aggregate_counts(SEXP input, Rcpp::IntegerVector grouping, int nthreads, bool binarize);
RcppExport SEXP _arbalist_aggregate_counts(SEXP inputSEXP, SEXP groupingSEXP, SEXP nthreadsSEXP, SEXP binarizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type grouping(groupingSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type binarize(binarizeSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregate_counts(input, grouping, nthreads, binarize));
    return rcpp_result_gen;
END_RCPP
}
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
// fragments_to_regions
SEXP fragments_to_regions(std::string fragment_file, std::string output_file, std::string output_group, Rcpp::CharacterVector seqnames, Rcpp::List region_ids, Rcpp::List region_starts, Rcpp::List region_ends, Rcpp::Nullable<Rcpp::CharacterVector> cellnames, int num_regions, int deflate_level, int chunk_dim);
RcppExport SEXP _arbalist_fragments_to_regions(SEXP fragment_fileSEXP, SEXP output_fileSEXP, SEXP output_groupSEXP, SEXP seqnamesSEXP, SEXP region_idsSEXP, SEXP region_startsSEXP, SEXP region_endsSEXP, SEXP cellnamesSEXP, SEXP num_regionsSEXP, SEXP deflate_levelSEXP, SEXP chunk_dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type fragment_file(fragment_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_file(output_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_group(output_groupSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type region_ids(region_idsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type region_starts(region_startsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type region_ends(region_endsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type cellnames(cellnamesSEXP);
    Rcpp::traits::input_parameter< int >::type num_regions(num_regionsSEXP);
    Rcpp::traits::input_parameter< int >::type deflate_level(deflate_levelSEXP);
    Rcpp::traits::input_parameter< int >::type chunk_dim(chunk_dimSEXP);
    rcpp_result_gen = Rcpp::wrap(fragments_to_regions(fragment_file, output_file, output_group, seqnames, region_ids, region_starts, region_ends, cellnames, num_regions, deflate_level, chunk_dim));
    return rcpp_result_gen;
END_RCPP
}
// fragments_to_tiles
SEXP fragments_to_tiles(std::string fragment_file, int tile_size, std::string output_file, std::string output_group, Rcpp::IntegerVector seqlengths, Rcpp::CharacterVector seqnames, Rcpp::Nullable<Rcpp::CharacterVector> cellnames, int deflate_level, int chunk_dim);
RcppExport SEXP _arbalist_fragments_to_tiles(SEXP fragment_fileSEXP, SEXP tile_sizeSEXP, SEXP output_fileSEXP, SEXP output_groupSEXP, SEXP seqlengthsSEXP, SEXP seqnamesSEXP, SEXP cellnamesSEXP, SEXP deflate_levelSEXP, SEXP chunk_dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type fragment_file(fragment_fileSEXP);
    Rcpp::traits::input_parameter< int >::type tile_size(tile_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_file(output_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_group(output_groupSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seqlengths(seqlengthsSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type cellnames(cellnamesSEXP);
    Rcpp::traits::input_parameter< int >::type deflate_level(deflate_levelSEXP);
    Rcpp::traits::input_parameter< int >::type chunk_dim(chunk_dimSEXP);
    rcpp_result_gen = Rcpp::wrap(fragments_to_tiles(fragment_file, tile_size, output_file, output_group, seqlengths, seqnames, cellnames, deflate_level, chunk_dim));
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
    {"_arbalist_aggregate_counts", (DL_FUNC) &_arbalist_aggregate_counts, 4},
    {"_arbalist_count_fragment_size_distributions", (DL_FUNC) &_arbalist_count_fragment_size_distributions, 1},
    {"_arbalist_fragments_to_regions", (DL_FUNC) &_arbalist_fragments_to_regions, 11},
    {"_arbalist_fragments_to_tiles", (DL_FUNC) &_arbalist_fragments_to_tiles, 9},
    {"_arbalist_irlba_realized", (DL_FUNC) &_arbalist_irlba_realized, 4},
    {"_arbalist_irlba_tatami", (DL_FUNC) &_arbalist_irlba_tatami, 4},
    {"_arbalist_lsi_matrix_stats", (DL_FUNC) &_arbalist_lsi_matrix_stats, 2},
    {"_arbalist_create_pseudobulk_file", (DL_FUNC) &_arbalist_create_pseudobulk_file, 3},
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
