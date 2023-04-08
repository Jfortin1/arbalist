#include "config.h"

#include "Rcpp.h"
#include "tatamize.h"

Rcpp::NumericMatrix aggregate_counts(SEXP input, Rcpp::IntegerVector grouping, int nthreads) {
    auto shared = extract_NumericMatrix_shared(input);
    size_t NR = shared->nrow();
    size_t NC = shared->ncol();
    int ngroups = *std::max_element(grouping.begin(), grouping.end()) + 1;

    Rcpp::NumericMatrix output(NR, ngroups);

    if (shared->prefer_rows()) {
        std::vector<double> buffer(NR * ngroups);

        run_parallel3(shared->nrow(), [&](int thread, size_t start, size_t end) -> void {
            std::vector<double> vbuffer(NC);
            auto wrk = shared->new_row_workspace();

            if (shared->sparse()) {
                std::vector<int> ibuffer(NC);
                for (size_t r = 0; r < NR; ++r) { 
                    auto ptr = buffer.data() + r * ngroups;
                    auto range = shared->sparse_row(r, vbuffer.data(), ibuffer.data(), wrk.get());
                    for (size_t i = 0; i < range.number; i++) {
                        ptr[grouping[range.index[i]]] += range.value[i];
                    }
                }

            } else {
                for (size_t r = 0; r < NR; ++r) { 
                    auto ptr = buffer.data() + r * ngroups;
                    auto found = shared->row(r, vbuffer.data(), wrk.get());
                    for (size_t i = 0; i < NC; i++) {
                        ptr[grouping[i]] += found[i];
                    }
                }
            }
        }, nthreads);

        // Applying a transposition.
        for (size_t r = 0; r < NR; ++r) {
            auto left = r;
            auto right = ngroups * r;
            for (int c = 0; c < ngroups; ++c, left += NR, ++right) {
                output[left] = buffer[right];
            }
        }

    } else {
        std::vector<std::vector<double> > partials(nthreads);
        std::vector<size_t> starts(nthreads, -1), lengths(nthreads, -1);

        run_parallel3(shared->nrow(), [&](int thread, size_t start, size_t end) -> void {
            size_t length = end - start;
            auto wrk = shared->new_column_workspace(start, length);
            starts[thread] = start;
            lengths[thread] = length;

            std::vector<double> vbuffer(length);
            auto& doutput = partials[thread];
            doutput.resize(NR * ngroups);

            if (shared->sparse()) {
                std::vector<int> ibuffer(length);
                for (size_t c = 0; c < NC; ++c) { 
                    auto range = shared->sparse_column(c, vbuffer.data(), ibuffer.data(), wrk.get());
                    auto ptr = doutput.data() + grouping[c] * NR;
                    for (size_t i = 0; i < range.number; i++) {
                        ptr[range.index[i]] += range.value[i];
                    }
                }

            } else {
                for (size_t c = 0; c < NC; ++c) { 
                    auto found = shared->column(c, wrk.get());
                    auto ptr = doutput.data() + grouping[c] * NR;
                    for (size_t i = 0; i < NR; i++) {
                        ptr[start + i] += found[i];
                    }
                }
            }
        }, nthreads);

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
