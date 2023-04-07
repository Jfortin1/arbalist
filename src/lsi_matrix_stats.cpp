#include "config.h"

#include "Rcpp.h"
#include "tatamize.h"
#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::List lsi_matrix_stats(SEXP mat, int nthreads) {
    auto shared = extract_NumericMatrix_shared(mat);
    size_t NR = shared->nrow(), NC = shared->ncol();
    Rcpp::NumericVector sums(NC);
    Rcpp::IntegerVector detected(NR);

    std::vector<std::vector<int> > tmp_detected(nthreads);

    if (shared->prefer_rows()) {
        std::vector<std::vector<double> > tmp_sums(nthreads);

        run_parallel3(shared->ncol(), [&](int thread, size_t start, size_t end) -> void {
            size_t length = end - start;
            auto wrk = shared->new_row_workspace(start, length);
            std::vector<double> vbuffer(length);

            auto& soutput = tmp_sums[thread];
            soutput.resize(NR);
            auto& doutput = tmp_detected[thread];
            doutput.resize(NR);

            if (shared->sparse()) {
                std::vector<int> ibuffer(length);
                for (size_t r = 0; r < NR; ++r) {
                    auto range = shared->sparse_row(r, vbuffer.data(), ibuffer.data(), wrk.get()); 

                    int total = 0;
                    for (size_t i = 0; i < range.number; ++i) {
                        total += range.value[i] > 0;
                        soutput[range.index[i]] += range.value[i];
                    }
                    doutput[r] = total;
                }
            } else {
                for (size_t r = 0; r < NR; ++r) {
                    auto found = shared->row(r, vbuffer.data(), wrk.get()); 
                    int total = 0;
                    for (size_t i = 0; i < length; ++i) {
                        total += found[i] > 0;
                        soutput[i + start] += found[i];
                    }
                    doutput[r] = total;
                }
            }

            std::copy(soutput.begin() + start, soutput.begin() + end, sums.begin() + start);
        }, nthreads);

    } else {
        run_parallel3(shared->ncol(), [&](int thread, size_t start, size_t end) -> void {
            auto wrk = shared->new_column_workspace();
            std::vector<double> vbuffer(NR);
            auto& doutput = tmp_detected[thread];
            doutput.resize(NR);

            if (shared->sparse()) {
                std::vector<int> ibuffer(NR);
                for (size_t c = start; c < end; ++c) {
                    auto range = shared->sparse_column(c, vbuffer.data(), ibuffer.data(), wrk.get()); 
                    double total = 0;
                    for (size_t i = 0; i < range.number; ++i) {
                        total += range.value[i];
                        doutput[range.index[i]] += range.value[i] > 0;
                    }
                    sums[c] = total;
                }

            } else {
                for (size_t c = start; c < end; ++c) {
                    auto found = shared->column(c, vbuffer.data(), wrk.get()); 
                    double total = 0;
                    for (size_t i = 0; i < NR; ++i) {
                        total += found[i];
                        doutput[i + start] += found[i] > 0;
                    }
                    sums[c] = total;
                }
            }

        }, nthreads);
    }

    for (const auto& tmp : tmp_detected) {
        if (!tmp.empty()) { // i.e., a thread was assigned and used.
            for (size_t r = 0; r < NR; ++r) {
                detected[r] += tmp[r];

            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums,
        Rcpp::Named("frequency") = detected
    );
}
