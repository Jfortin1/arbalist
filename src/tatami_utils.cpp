#include "config.h"

#include <array>
#include "tatami/tatami.hpp"
#include "tatami/ext/hdf5/load_hdf5_matrix.hpp"
#include "Rcpp.h"
#include "tatamize.h"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector tatami_dim(SEXP input) {
    auto shared = extract_NumericMatrix_shared(input);
    return Rcpp::IntegerVector::create(shared->nrow(), shared->ncol());
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_column(SEXP input, int i) {
    auto shared = extract_NumericMatrix_shared(input);
    Rcpp::NumericVector output(shared->nrow());
    auto wrk = shared->new_column_workspace();
    shared->column_copy(i-1, static_cast<double*>(output.begin()), wrk.get());
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_row(SEXP input, int i) {
    auto shared = extract_NumericMatrix_shared(input);
    Rcpp::NumericVector output(shared->ncol());
    auto wrk = shared->new_row_workspace();
    shared->row_copy(i-1, static_cast<double*>(output.begin()), wrk.get());
    return output;
}
