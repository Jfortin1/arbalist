#include "config.h"

#include "Rcpp.h"
#include "irlba/irlba.hpp"
#include "irlba/parallel.hpp"
#include "tatamize.h"

#include <vector>
#include <algorithm>
#include <iostream>

// STOLEN FROM LIBSCRAN:

struct SparseComponents {
    std::vector<size_t> ptrs;
    std::vector<double> values;
    std::vector<int> indices;
};

template<typename T, typename IDX>
SparseComponents sparse_by_row(const tatami::Matrix<T, IDX>* mat, int nthreads) {
    SparseComponents output;
    size_t NR = mat->nrow(), NC = mat->ncol();
    auto& ptrs = output.ptrs;
    ptrs.resize(NR + 1);

    /*** First round, to fetch the number of zeros in each row. ***/
    ARBALIST_CUSTOM_PARALLEL(NR, [&](size_t start, size_t end) -> void {
        std::vector<double> xbuffer(NC);
        std::vector<int> ibuffer(NC);
        auto wrk = mat->new_row_workspace();
        for (size_t r = start; r < end; ++r) {
            auto range = mat->sparse_row(r, xbuffer.data(), ibuffer.data(), wrk.get());
            ptrs[r + 1] = range.number;
        }
    }, nthreads);

    /*** Second round, to populate the vectors. ***/
    for (size_t r = 0; r < NR; ++r) {
        ptrs[r + 1] += ptrs[r];
    }
    output.values.resize(ptrs.back());
    output.indices.resize(ptrs.back());

    ARBALIST_CUSTOM_PARALLEL(NR, [&](size_t start, size_t end) -> void {
        auto wrk = mat->new_row_workspace();
        for (size_t r = start; r < end; ++r) {
            auto offset = ptrs[r];
            mat->sparse_row_copy(r, output.values.data() + offset, output.indices.data() + offset, wrk.get(), tatami::SPARSE_COPY_BOTH);
        }
    }, nthreads);

    return output;
}

template<typename T, typename IDX>
SparseComponents sparse_by_column(const tatami::Matrix<T, IDX>* mat, int nthreads) {
    size_t NR = mat->nrow(), NC = mat->ncol();

    /*** First round, to fetch the number of zeros in each row. ***/
    size_t cols_per_thread = std::ceil(static_cast<double>(NC) / nthreads);
    std::vector<std::vector<size_t> > threaded_nonzeros_per_row(nthreads);

    ARBALIST_CUSTOM_PARALLEL(nthreads, [&](int start, int end) -> void { // Trivial allocation of one job per thread.
        for (int t = start; t < end; ++t) {
            size_t startcol = cols_per_thread * t, endcol = std::min(startcol + cols_per_thread, NC);
            if (startcol < endcol) {
                std::vector<size_t>& nonzeros_per_row = threaded_nonzeros_per_row[t];
                nonzeros_per_row.resize(NR);
                std::vector<double> xbuffer(NR);
                std::vector<int> ibuffer(NR);
                auto wrk = mat->new_column_workspace();

                for (size_t c = startcol; c < endcol; ++c) {
                    auto range = mat->sparse_column(c, xbuffer.data(), ibuffer.data(), wrk.get());
                    for (size_t i = 0; i < range.number; ++i) {
                        ++(nonzeros_per_row[range.index[i]]);
                    }
                }
            }
        }
    }, nthreads);

    // There had better be at least one thread!
    std::vector<size_t> nonzeros_per_row = std::move(threaded_nonzeros_per_row[0]);
    for (int t = 1; t < nthreads; ++t) {
        auto it = nonzeros_per_row.begin();
        for (auto x : threaded_nonzeros_per_row[t]) {
            *it += x;
            ++it;
        }
    }

    /*** Second round, to populate the vectors. ***/
    SparseComponents output;
    output.ptrs.resize(NR + 1);
    size_t total_nzeros = 0;
    for (size_t r = 0; r < NR; ++r) {
        total_nzeros += nonzeros_per_row[r];
        output.ptrs[r + 1] = total_nzeros;

    }
    output.values.resize(total_nzeros);
    output.indices.resize(total_nzeros);

    // Splitting by row this time, because columnar extraction can't be done safely.
    nthreads = 1;
    size_t rows_per_thread = std::ceil(static_cast<double>(NR) / nthreads);
    auto ptr_copy = output.ptrs;

    ARBALIST_CUSTOM_PARALLEL(nthreads, [&](int start, int end) -> void { // Trivial allocation of one job per thread.
    for (int t = start; t < end; ++t) {
        size_t startrow = rows_per_thread * t, endrow = std::min(startrow + rows_per_thread, NR);

        if (startrow < endrow) {
            size_t length = endrow - startrow;
            auto wrk = mat->new_column_workspace(startrow, length);
            std::vector<double> xbuffer(length);
            std::vector<int> ibuffer(length);

            for (size_t c = 0; c < NC; ++c) {
                auto range = mat->sparse_column(c, xbuffer.data(), ibuffer.data(), wrk.get());
                for (size_t i = 0; i < range.number; ++i) {
                    auto r = range.index[i];
                    auto& offset = ptr_copy[r];
                    output.values[offset] = range.value[i];
                    output.indices[offset] = c;
                    ++offset;
                }
            }
        }
    }
    }, nthreads);

    return output;
}

template<typename T, typename IDX>
SparseComponents extract_sparse_for_pca(const tatami::Matrix<T, IDX>* mat, int nthreads) {
    if (mat->prefer_rows()) {
        return sparse_by_row(mat, nthreads);
    } else {
        return sparse_by_column(mat, nthreads);
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List irlba_realized(SEXP input, int rank, int nthreads, int seed) {
    auto shared = extract_NumericMatrix_shared(input);
    size_t NR = shared->nrow();
    size_t NC = shared->ncol();

    // Extracting as row-major sparse matrix and transposing it.
    auto info = extract_sparse_for_pca(shared.get(), nthreads);
    irlba::ParallelSparseMatrix<> mat(NC, NR, std::move(info.values), std::move(info.indices), std::move(info.ptrs), nthreads);

    irlba::EigenThreadScope t(nthreads);
    irlba::Irlba runner;
    runner.set_number(rank);
    runner.set_seed(seed);

    Eigen::MatrixXd u;
    Eigen::MatrixXd v;
    Eigen::VectorXd d;
    runner.run(mat, u, v, d);

    Rcpp::NumericMatrix out_u(u.rows(), u.cols());
    std::copy(u.data(), u.data() + u.size(), out_u.begin());
    Rcpp::NumericMatrix out_v(v.rows(), v.cols());
    std::copy(v.data(), v.data() + v.size(), out_v.begin());
    Rcpp::NumericVector out_d(d.begin(), d.end());

    return Rcpp::List::create(
        Rcpp::Named("u") = out_u,
        Rcpp::Named("v") = out_v,
        Rcpp::Named("d") = out_d
    );
}
