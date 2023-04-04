#include "tatamize.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

//[[Rcpp::export(rng=false)]]
SEXP apply_subset(SEXP input, Rcpp::IntegerVector subset, bool row) {
    auto shared = extract_NumericMatrix_shared(input);

    // Is this a contiguous block?
    bool consecutive = true;
    for (size_t i = 1, end = subset.size(); i < end; ++i) {
        if (subset[i] - subset[i - 1] != 1) {
            consecutive = false;
            break;
        }
    }

    if (consecutive) {
        int start = (subset.size() ? subset[0] : 0);
        int end = (subset.size() ? subset[subset.size() - 1] + 1 : 0);
        if (row) {
            return new_MatrixChan(tatami::make_DelayedSubsetBlock<0>(shared, start, end));
        } else {
            return new_MatrixChan(tatami::make_DelayedSubsetBlock<1>(shared, start, end));
        }
    }

    // Otherwise, we get to 1-based indices.
    std::vector<int> resub(subset.begin(), subset.end());
    for (auto& x : resub) {
        --x; 
    } 

    if (row) {
        return new_MatrixChan(tatami::make_DelayedSubset<0>(shared, std::move(resub)));
    } else {
        return new_MatrixChan(tatami::make_DelayedSubset<1>(shared, std::move(resub)));
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_transpose(SEXP input) {
    auto shared = extract_NumericMatrix_shared(input);
    return new_MatrixChan(tatami::make_DelayedTranspose(shared));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_bind(Rcpp::List input, bool row) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    collected.reserve(input.size());

    for (size_t i = 0, end = input.size(); i < end; ++i) {
        Rcpp::RObject current = input[i];
        collected.push_back(extract_NumericMatrix_shared(current));
    }

    if (row) {
        return new_MatrixChan(tatami::make_DelayedBind<0>(std::move(collected)));
    } else {
        return new_MatrixChan(tatami::make_DelayedBind<1>(std::move(collected)));
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_addition(SEXP input, Rcpp::NumericVector val, bool row) {
    auto shared = extract_NumericMatrix_shared(input);
    if (val.size() == 1) {
        return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::DelayedAddScalarHelper(val[0])));
    }

    if (row) {
        return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedAddVectorHelper<0>(std::move(val))));
    } else {
        return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedAddVectorHelper<1>(std::move(val))));
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_multiplication(SEXP input, Rcpp::NumericVector val, bool row) {
    auto shared = extract_NumericMatrix_shared(input);
    if (val.size() == 1) {
        return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::DelayedAddScalarHelper(val[0])));
    }

    if (row) {
        return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedAddVectorHelper<0>(std::move(val))));
    } else {
        return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedAddVectorHelper<1>(std::move(val))));
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_subtraction(SEXP input, Rcpp::NumericVector val, bool right, bool row) {
    auto shared = extract_NumericMatrix_shared(input);
    if (val.size() == 1) {
        if (right) {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::DelayedSubtractScalarHelper<true>(val[0])));
        } else {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::DelayedSubtractScalarHelper<false>(val[0])));
        }
    }

    if (right) {
        if (row) {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<true, 0>(std::move(val))));
        } else {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<false, 1>(std::move(val))));
        }
    } else {
        if (row) {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<true, 0>(std::move(val))));
        } else {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<false, 1>(std::move(val))));
        }
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_division(SEXP input, Rcpp::NumericVector val, bool right, bool row) {
    auto shared = extract_NumericMatrix_shared(input);
    if (val.size() == 1) {
        if (right) {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::DelayedDivideScalarHelper<true>(val[0])));
        } else {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::DelayedDivideScalarHelper<false>(val[0])));
        }
    }

    if (right) {
        if (row) {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<true, 0>(std::move(val))));
        } else {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<false, 1>(std::move(val))));
        }
    } else {
        if (row) {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<true, 0>(std::move(val))));
        } else {
            return new_MatrixChan(tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<false, 1>(std::move(val))));
        }
    }
}
