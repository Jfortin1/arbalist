#include "config.h"

#include "Rcpp.h"
#include "tatamize.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix aggregate_counts(SEXP input, Rcpp::IntegerVector grouping, int nthreads) {
    auto shared = extract_NumericMatrix_shared(input);
    int NR = shared->nrow();
    int NC = shared->ncol();
    int ngroups = *std::max_element(grouping.begin(), grouping.end()) + 1;
    
    //Rprintf("NR: %i, NC: %i, ngroups: %i\n", NR, NC, ngroups);

    Rcpp::NumericMatrix output(NR, ngroups);

    if (shared->prefer_rows()) {
        std::vector<double> buffer(NR * ngroups);

        ARBALIST_CUSTOM_PARALLEL([&](int thread, int start, int length) -> void {
          std::vector<double> vbuffer(length);

          if (shared->sparse()) {
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
          //Rprintf("thread: %i, start: %i, end: %i\n", thread, start, end);
          int length = end - start;
          starts[thread] = start;
          lengths[thread] = length;

          std::vector<double> vbuffer(length);
          auto& doutput = partials[thread];
          doutput.resize(NR * ngroups);

          if (shared->sparse()) {
              auto wrk = tatami::consecutive_extractor<false, true>(shared.get(), start, length);
              std::vector<int> ibuffer(length);
              for (int c = 0; c < NC; ++c) { 
                  auto range = wrk->fetch(c, vbuffer.data(), ibuffer.data());
                  auto ptr = doutput.data() + grouping[c] * NR;
                  for (int i = 0; i < range.number; i++) {
                      ptr[range.index[i]] += range.value[i];
                    //Rprintf("c: %i, grouping.c: %i, i: %i, range.number: %i, range.index: %i, range.value: %f, ptr.value: %f\n", c, grouping[c], i, range.number, range.index[i], range.value[i], ptr[range.index[i]]);
                    //Rprintf("%i %i %i %i %i %i %i %i %i %i %i %i\n", int (doutput[0]), int (doutput[1]),int (doutput[2]), int (doutput[3]),int (doutput[4]),int (doutput[5]),int (doutput[6]),int (doutput[7]),int (doutput[8]),int (doutput[9]),int (doutput[10]),int (doutput[11]));
                    //Rprintf("%i %i %i %i %i %i %i %i %i %i %i %i\n", int (doutput[12]), int (doutput[13]),int (doutput[14]), int (doutput[15]),int (doutput[16]),int (doutput[17]),int (doutput[18]),int (doutput[19]),int (doutput[20]),int (doutput[21]),int (doutput[22]),int (doutput[23]));
                    //Rprintf("%i %i %i %i %i %i %i %i %i %i %i %i\n", int (doutput[24]), int (doutput[25]),int (doutput[26]), int (doutput[27]),int (doutput[28]),int (doutput[29]),int (doutput[30]),int (doutput[31]),int (doutput[32]),int (doutput[33]),int (doutput[34]),int (doutput[35]));
                    //Rprintf("%i %i %i %i %i %i %i %i %i %i %i %i\n", int (doutput[36]), int (doutput[37]),int (doutput[38]), int (doutput[39]),int (doutput[40]),int (doutput[41]),int (doutput[42]),int (doutput[43]),int (doutput[44]),int (doutput[45]),int (doutput[46]),int (doutput[47]));
                  }
             }
          } else {
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