#include "config.h"

#include "Rcpp.h"
#include "irlba/irlba.hpp"
#include "irlba/parallel.hpp"

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
    IDX NR = mat->nrow(), NC = mat->ncol();
    auto& ptrs = output.ptrs;
    ptrs.resize(NR + 1);

    /*** First round, to fetch the number of zeros in each row. ***/
    tatami::Options opt;
    opt.sparse_extract_index = false;
    opt.sparse_extract_value = false;

    ARBALIST_CUSTOM_PARALLEL([&](int, IDX start, IDX len) -> void {
        auto wrk = tatami::consecutive_extractor<true, true>(mat, start, len, opt);
        for (IDX r = start, end = start + len; r < end; ++r) {
            auto range = wrk->fetch(r, NULL, NULL);
            ptrs[r + 1] = range.number;
        }
    }, NR, nthreads);

    /*** Second round, to populate the vectors. ***/
    for (IDX r = 0; r < NR; ++r) {
        ptrs[r + 1] += ptrs[r];
    }
    output.values.resize(ptrs.back());
    output.indices.resize(ptrs.back());

    ARBALIST_CUSTOM_PARALLEL([&](int, IDX start, IDX len) -> void {
        auto wrk = tatami::consecutive_extractor<true, true>(mat, start, len);
        for (IDX r = start, end = start + len; r < end; ++r) {
            auto offset = ptrs[r];
            wrk->fetch_copy(r, output.values.data() + offset, output.indices.data() + offset);
        }
    }, NR, nthreads);

    return output;
}

template<typename T, typename IDX>
SparseComponents sparse_by_column(const tatami::Matrix<T, IDX>* mat, int nthreads) {
    IDX NR = mat->nrow(), NC = mat->ncol();

    /*** First round, to fetch the number of zeros in each row. ***/
    tatami::Options opt;
    opt.sparse_extract_index = true;
    opt.sparse_extract_value = false;

    std::vector<std::vector<IDX> > threaded_nonzeros_per_row(nthreads); // create separate vectors to avoid false sharing.

    ARBALIST_CUSTOM_PARALLEL([&](int t, IDX start, IDX len) -> void {
        auto wrk = tatami::consecutive_extractor<false, true>(mat, start, len, opt);
        auto& mine = threaded_nonzeros_per_row[t];
        mine.resize(NR);
        std::vector<IDX> ibuffer(NR);

        for (IDX i = start, end = start + len; i < end; ++i) {
            auto range = wrk->fetch(i, NULL, ibuffer.data());
            for (IDX j = 0; j < range.number; ++j) {
                ++(mine[range.index[j]]);
            }
        }
    }, NC, nthreads);

    std::vector<IDX> nonzeros_per_row = std::move(threaded_nonzeros_per_row[0]);
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
    IDX total_nzeros = 0;
    for (IDX r = 0; r < NR; ++r) {
        total_nzeros += nonzeros_per_row[r];
        output.ptrs[r + 1] = total_nzeros;

    }
    output.values.resize(total_nzeros);
    output.indices.resize(total_nzeros);
    auto ptr_copy = output.ptrs;

    // Splitting by row this time, because columnar extraction can't be done safely.
    ARBALIST_CUSTOM_PARALLEL([&](int t, int start, int len) -> void { 
        auto wrk = tatami::consecutive_extractor<false, true>(mat, 0, NC, start, len);
        std::vector<T> vbuffer(len);
        std::vector<IDX> ibuffer(len);

        for (IDX c = 0; c < NC; ++c) {
            auto range = wrk->fetch(c, vbuffer.data(), ibuffer.data());
            for (IDX i = 0; i < range.number; ++i) {
                auto r = range.index[i];
                auto& offset = ptr_copy[r];
                output.values[offset] = range.value[i];
                output.indices[offset] = c;
                ++offset;
            }
        }
    }, NR, nthreads);

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
    auto converted = Rtatami::BoundNumericPointer(input);
    const auto& shared = converted->ptr;
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
