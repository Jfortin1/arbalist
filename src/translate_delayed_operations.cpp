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
    int start = (subset.size() ? subset[0] - 1 : 0);
    int end = (subset.size() ? subset[subset.size() - 1] : 0);
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
    return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricAddScalar<>(val[0])));
  }
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricAddVector<>(std::move(val), row)));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_multiplication(SEXP input, Rcpp::NumericVector val, bool row) {
  auto shared = extract_NumericMatrix_shared(input);
  if (val.size() == 1) {
    return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricMultiplyScalar<>(val[0])));
  }
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricMultiplyVector<>(std::move(val), row)));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_subtraction(SEXP input, Rcpp::NumericVector val, bool right, bool row) {
  auto shared = extract_NumericMatrix_shared(input);
  if (val.size() == 1) {
    if (right) {
      return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricSubtractScalar<true>(val[0])));
    } else {
      return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricSubtractScalar<false>(val[0])));
    }
  }
  
  if (right) {
    return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricSubtractVector<true>(std::move(val), row)));
  } else {
    return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricSubtractVector<false>(std::move(val), row)));
  }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_division(SEXP input, Rcpp::NumericVector val, bool right, bool row) {
  auto shared = extract_NumericMatrix_shared(input);
  if (val.size() == 1) {
    if (right) {
      return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricDivideScalar<true>(val[0])));
    } else {
      return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricDivideScalar<false>(val[0])));
    }
  }
  
  if (right) {
    return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricDivideVector<true>(std::move(val), row)));
  } else {
    return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricDivideVector<false>(std::move(val), row)));
  }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_log(SEXP input, double base) {
  auto shared = extract_NumericMatrix_shared(input);
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::DelayedUnaryIsometricLog<>(base)));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_log1p(SEXP input) {
  auto shared = extract_NumericMatrix_shared(input);
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::DelayedUnaryIsometricLog1p<>()));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_abs(SEXP input) {
  auto shared = extract_NumericMatrix_shared(input);
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::DelayedUnaryIsometricAbs<>()));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_sqrt(SEXP input) {
  auto shared = extract_NumericMatrix_shared(input);
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::DelayedUnaryIsometricSqrt<>()));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_round(SEXP input) {
  auto shared = extract_NumericMatrix_shared(input);
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::DelayedUnaryIsometricRound<>()));
}

//[[Rcpp::export(rng=false)]]
SEXP apply_exp(SEXP input) {
  auto shared = extract_NumericMatrix_shared(input);
  return new_MatrixChan(tatami::make_DelayedUnaryIsometricOperation(shared, tatami::DelayedUnaryIsometricExp<>()));
}