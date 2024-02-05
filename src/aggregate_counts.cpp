#include "config.h"

#include "Rcpp.h"
#include "tatamize.h"
#include <iostream>

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix aggregate_counts(SEXP input, Rcpp::IntegerVector grouping, int nthreads) {
    auto shared = extract_NumericMatrix_shared(input);
    int NR = shared->nrow();
    int NC = shared->ncol();
    int ngroups = *std::max_element(grouping.begin(), grouping.end()) + 1;
    std::printf("NC %i\n",NC);
    std::printf("NR %i\n",NR);
    std::printf("ngroups %i\n",ngroups);
    
    Rcpp::NumericMatrix output(NR, ngroups);

    if (shared->prefer_rows()) {
        std::vector<double> buffer(NR * ngroups);

        ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int length) -> void {
          std::vector<double> vbuffer(length);

          if (shared->sparse()) {
            std::printf("1 %i, %i, %i",thread,start, length);
              auto wrk = tatami::consecutive_extractor<true, true>(shared.get(), 0, NR, start, length);
              std::vector<int> ibuffer(NC);
              for (int r = 0; r < NR; ++r) { 
                    auto ptr = buffer.data() + r * ngroups;
                    auto range = wrk->fetch(r, vbuffer.data(), ibuffer.data());
                    for (int i = 0; i < range.number; i++) {
                        ptr[grouping[range.index[i]]] += range.value[i];
                  }
              }

          } else {
            std::printf("2 %i, %i, %i",thread,start, length);
              auto wrk = tatami::consecutive_extractor<true, false>(shared.get(), 0, NR, start, length);
              for (int r = 0; r < NR; ++r) { 
                  auto ptr = buffer.data() + r * ngroups;
                  auto found = wrk->fetch(r, vbuffer.data());
                  for (int i = 0; i < NC; i++) {
                      ptr[grouping[i]] += found[i];
                  }
              }
          }
      }, NC, nthreads);

      // Applying a transposition.
      for (int r = 0; r < NR; ++r) {
          auto left = r;
          auto right = ngroups * r;
          for (int c = 0; c < ngroups; ++c, left += NR, ++right) {
              output[left] = buffer[right];
          }
        }

    } else {
        std::vector<std::vector<double> > partials(nthreads);
        std::vector<int> starts(nthreads, -1), lengths(nthreads, -1);

        ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int end) -> void {
          int length = end - start;
          starts[thread] = start;
          lengths[thread] = length;

          std::vector<double> vbuffer(length);
          auto& doutput = partials[thread];
          doutput.resize(NR * ngroups);

          if (shared->sparse()) {
            std::printf("3\n");
              auto wrk = tatami::consecutive_extractor<false, true>(shared.get(), start, length);
              std::vector<int> ibuffer(length);
              for (int c = 0; c < NC; ++c) { 
                  auto range = wrk->fetch(c, vbuffer.data(), ibuffer.data());
                  auto ptr = doutput.data() + grouping[c] * NR;
                  for (int i = 0; i < range.number; i++) {
                    std::printf("%i, %i, %i, %i\n",i, range.number, range.index[i],range.value[i]);
                      ptr[range.index[i]] += range.value[i];
                  }
              }
          } else {
            std::printf("4 %i, %i, %i",thread,start, length);
              auto wrk = tatami::consecutive_extractor<false, false>(shared.get(), start, length);
              for (int c = 0; c < NC; ++c) { 
                  auto found = wrk->fetch(c, vbuffer.data());
                  auto ptr = doutput.data() + grouping[c] * NR;
                  for (int i = 0; i < NR; i++) {
                      ptr[start + i] += found[i];
                  }
              }
          }
      }, NR, nthreads);

        // Copying values over.
        for (int c = 0; c < ngroups; ++c) {
            for (int t = 0; t < nthreads; ++t) {
                if (starts[t] == -1) { 
                    continue; 
                }
                auto source = partials[t].begin() + c * lengths[t];
                std::copy(source, source + lengths[t], output.begin() + c * NR + starts[t]);
            }
        }
        
      }
    
      return output;
}