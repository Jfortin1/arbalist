#include "config.h"

#include "Rcpp.h"
#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::List data_distribtion(SEXP mat, int max) {
  auto converted = Rtatami::BoundNumericPointer(mat);
  const auto& shared = converted->ptr;
  int NR = shared->nrow(), NC = shared->ncol();
  std::vector<std::vector<int> > counts(NC);
  
  if (shared->prefer_rows()) {
    ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int length) -> void {
      std::vector<double> vbuffer(length);
      
      auto& doutput = counts[thread];
      doutput.resize(max);
      
      if (shared->sparse()) {
        auto wrk = tatami::consecutive_extractor<true, true>(shared.get(), 0, NR, start, length);
        std::vector<int> ibuffer(length);
        for (int r = 0; r < NR; ++r) {
          auto range = wrk->fetch(r, vbuffer.data(), ibuffer.data());
          for (int i = 0; i < range.number; ++i) {
            doutput[range.value[i]-1] += 1;
          }
        }
        
      } else {
        auto wrk = tatami::consecutive_extractor<true, false>(shared.get(), 0, NR, start, length);
        for (int r = 0; r < NR; ++r) {
          auto found = wrk->fetch(r, vbuffer.data());
          for (int i = 0; i < length; ++i) {
            doutput[found[i] -1] += 1;
          }
        }
      }
    }, NC, NC);
    
  } else {
    ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int length) -> void {
      std::vector<double> vbuffer(NR);
      auto& doutput = counts[thread];
      doutput.resize(max);
      
      if (shared->sparse()) {
        auto wrk = tatami::consecutive_extractor<false, true>(shared.get(), start, length);
        std::vector<int> ibuffer(NR);
        for (int c = start, end = start + length; c < end; ++c) {
          auto range = wrk->fetch(c, vbuffer.data(), ibuffer.data());
          for (int i = 0; i < range.number; ++i) {
            doutput[range.value[i]-1] += 1;
          }
        }
      } else {
        auto wrk = tatami::consecutive_extractor<false, false>(shared.get(), start, length);
        for (int c = start, end = start + length; c < end; ++c) {
          auto found = wrk->fetch(c, vbuffer.data());
          for (int i = 0; i < NR; ++i) {
            doutput[found[i]-1] += 1;
          }
        }
      }
    }, NC, NC);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("distribution") = counts
  );
}