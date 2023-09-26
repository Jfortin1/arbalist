#include "config.h"

#include "Rcpp.h"
#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::List lsi_matrix_stats(SEXP mat, int nthreads) {
  auto converted = Rtatami::BoundNumericPointer(mat);
  const auto& shared = converted->ptr;
  int NR = shared->nrow(), NC = shared->ncol();
  Rcpp::NumericVector sums(NC);
  Rcpp::IntegerVector detected(NR);
  
  std::vector<std::vector<int> > tmp_detected(nthreads);
  
  if (shared->prefer_rows()) {
    ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int length) -> void {
      std::vector<double> vbuffer(length);
      
      std::vector<double> soutput(NC);
      auto& doutput = tmp_detected[thread];
      doutput.resize(NR);
      
      if (shared->sparse()) {
        auto wrk = tatami::consecutive_extractor<true, true>(shared.get(), 0, NR, start, length);
        std::vector<int> ibuffer(length);
        for (int r = 0; r < NR; ++r) {
          auto range = wrk->fetch(r, vbuffer.data(), ibuffer.data());
          int total = 0;
          for (int i = 0; i < range.number; ++i) {
            total += range.value[i] > 0;
            soutput[range.index[i]] += range.value[i];
          }
          doutput[r] = total;
        }
        
      } else {
        auto wrk = tatami::consecutive_extractor<true, false>(shared.get(), 0, NR, start, length);
        for (int r = 0; r < NR; ++r) {
          auto found = wrk->fetch(r, vbuffer.data());
          int total = 0;
          for (int i = 0; i < length; ++i) {
            total += found[i] > 0;
            soutput[i + start] += found[i];
          }
          doutput[r] = total;
        }
      }
      
      std::copy_n(soutput.begin() + start, length, sums.begin() + start);
    }, NC, nthreads);
    
  } else {
    ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int length) -> void {
      std::vector<double> vbuffer(NR);
      auto& doutput = tmp_detected[thread];
      doutput.resize(NR);
      
      if (shared->sparse()) {
        auto wrk = tatami::consecutive_extractor<false, true>(shared.get(), start, length);
        std::vector<int> ibuffer(NR);
        for (int c = start, end = start + length; c < end; ++c) {
          auto range = wrk->fetch(c, vbuffer.data(), ibuffer.data());
          double total = 0;
          for (int i = 0; i < range.number; ++i) {
            total += range.value[i];
            doutput[range.index[i]] += range.value[i] > 0;
          }
          sums[c] = total;
        }
      } else {
        auto wrk = tatami::consecutive_extractor<false, false>(shared.get(), start, length);
        for (int c = start, end = start + length; c < end; ++c) {
          auto found = wrk->fetch(c, vbuffer.data());
          double total = 0;
          for (int i = 0; i < NR; ++i) {
            total += found[i];
            doutput[i + start] += found[i] > 0;
          }
          sums[c] = total;
        }
      }
      
    }, NC, nthreads);
  }
  
  for (const auto& tmp : tmp_detected) {
    if (!tmp.empty()) { // i.e., a thread was assigned and used.
      for (int r = 0; r < NR; ++r) {
        detected[r] += tmp[r];
        
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("sums") = sums,
    Rcpp::Named("frequency") = detected
  );
}