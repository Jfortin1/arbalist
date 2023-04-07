#include "config.h"

#include "Rcpp.h"
#include "irlba/irlba.hpp"
#include "irlba/parallel.hpp"
#include "tatamize.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>

struct TatamiWrapper {
    TatamiWrapper(const tatami::NumericMatrix* m, int n) : mat(m), nthreads(n) {}

private:
    const tatami::NumericMatrix* mat;
    int nthreads;

public:
    size_t rows() const {
        return mat->nrow();
    }

    size_t cols() const {
        return mat->ncol();
    }

public:
    struct Workspace {
        Workspace(const tatami::NumericMatrix* mat, int nthreads, bool multiply) : 
            values(nthreads), indices(nthreads) // make sure nthreads entries are available in indices, even if empty, as the dense methods still access them.
        { 
            if (mat->prefer_rows()) {
                size_t NC = mat->ncol();

                if (multiply || nthreads == 1) {
                    rowwork.reserve(nthreads);
                    for (int t = 0; t < nthreads; ++t) {
                        rowwork.push_back(mat->new_row_workspace());
                        values[t].resize(NC);
                        if (mat->sparse()) {
                            indices[t].resize(NC);
                        }
                    }

                } else {
                    size_t start = 0;
                    size_t per_thread = std::ceil(static_cast<double>(NC) / nthreads);
                    rowbwork.reserve(nthreads);
                    staging.resize(nthreads);

                    for (int t = 0; t < nthreads; ++t) {
                        size_t end = std::min(NC, start + per_thread);
                        size_t length = end - start;
                        rowbwork.push_back(mat->new_row_workspace(start, length, true));
                        start = end;

                        values[t].resize(length);
                        staging[t].resize(NC);
                        if (mat->sparse()) {
                            indices[t].resize(length);
                        }
                    }
                }

            } else {
                size_t NR = mat->nrow();

                if (!multiply || nthreads == 1) {
                    colwork.reserve(nthreads);
                    for (int t = 0; t < nthreads; ++t) {
                        colwork.push_back(mat->new_column_workspace(true));
                        values[t].resize(NR);
                        if (mat->sparse()) {
                            indices[t].resize(NR);
                        }
                    }

                } else {
                    size_t start = 0;
                    size_t per_thread = std::ceil(static_cast<double>(NR) / nthreads);
                    colbwork.reserve(nthreads);
                    staging.resize(nthreads);

                    for (int t = 0; t < nthreads; ++t) {
                        size_t end = std::min(NR, start + per_thread);
                        size_t length = end - start;
                        colbwork.push_back(mat->new_column_workspace(start, length, true));
                        start = end;

                        values[t].resize(length);
                        staging[t].resize(NR);
                        if (mat->sparse()) {
                            indices[t].resize(length);
                        }
                    }
                }
            }

            if (nthreads > 1) {
                jobs.reserve(nthreads);
            }
        }

        std::vector<std::vector<double> > values;
        std::vector<std::vector<int> > indices;

        std::vector<std::vector<double> > staging; // avoid false sharing during parallelized addition.

        std::vector<std::shared_ptr<tatami::RowWorkspace> > rowwork;
        std::vector<std::shared_ptr<tatami::ColumnWorkspace> > colwork;
        std::vector<std::shared_ptr<tatami::RowBlockWorkspace> > rowbwork;
        std::vector<std::shared_ptr<tatami::ColumnBlockWorkspace> > colbwork;

        std::vector<std::thread> jobs;
    };

    Workspace workspace() const {
        return Workspace(mat, nthreads, true);
    }

    Workspace adjoint_workspace() const {
        return Workspace(mat, nthreads, false);
    }

private:
    template<bool column>
    static auto& get_work(Workspace& work) {
        if constexpr(column) {
            return work.colwork;
        } else {
            return work.rowwork;
        }
    }

    template<bool column>
    static auto& get_block_work(Workspace& work) {
        if constexpr(column) {
            return work.colbwork;
        } else {
            return work.rowbwork;
        }
    }

    template<bool column>
    size_t get_dimension() const {
        if constexpr(column) {
            return mat->ncol();
        } else {
            return mat->nrow();
        }
    }

    template<bool column>
    size_t get_other_dimension() const {
        if constexpr(column) {
            return mat->nrow();
        } else {
            return mat->ncol();
        }
    }

    template<bool column, class TatamiWork>
    auto get_sparse_range(size_t x, double* vbuffer, int* ibuffer, TatamiWork* work) const {
        if constexpr(column) {
            return mat->sparse_column(x, vbuffer, ibuffer, work);
        } else {
            return mat->sparse_row(x, vbuffer, ibuffer, work);
        }
    }

    template<bool column, class TatamiWork>
    auto get_dense(size_t x, double* vbuffer, TatamiWork* work) const {
        if constexpr(column) {
            return mat->column(x, vbuffer, work);
        } else {
            return mat->row(x, vbuffer, work);
        }
    }

private:
    template<bool column, bool sparse, class Right>
    void running_serial(const Right& rhs, Eigen::VectorXd& output, Workspace& work) const {
        auto& current = get_work<column>(work).front();
        auto vbuffer = work.values.front().data();
        auto ibuffer = work.indices.front().data();

        output.setZero();
        size_t end = get_dimension<column>();
        size_t other = get_other_dimension<column>();

        for (size_t x = 0; x < end; ++x) {
            auto mult = rhs.coeff(x);

            if constexpr(sparse) {
                auto found = get_sparse_range<column>(x, vbuffer, ibuffer, current.get());
                for (size_t i = 0; i < found.number; ++i) {
                    output.coeffRef(found.index[i]) += found.value[i] * mult;
                }
            } else {
                auto found = get_dense<column>(x, vbuffer, current.get());
                for (size_t i = 0; i < other; ++i) {
                    output.coeffRef(i) += found[i] * mult;
                }
            }
        }
    }

    template<bool column, bool sparse, class Right>
    void running_parallel(const Right& rhs, Eigen::VectorXd& output, Workspace& work) const {
        auto& jobs = work.jobs;
        auto& currentwork = get_block_work<column>(work);
        output.setZero();
        size_t end = get_dimension<column>();
        jobs.clear();

        for (int t = 0; t < nthreads; ++t) {
            jobs.emplace_back([&](int thread) -> void {
                auto& current = currentwork[thread];
                auto vbuffer = work.values[thread].data();
                auto ibuffer = work.indices[thread].data();

                auto& sums = work.staging[thread]; // do the sums in a staging area to avoid false sharing problems.
                std::fill(sums.begin(), sums.end(), 0);
                const auto& pos = current->block();

                for (size_t x = 0; x < end; ++x) {
                    auto mult = rhs.coeff(x);

                    if constexpr(sparse) {
                        auto found = get_sparse_range<column>(x, vbuffer, ibuffer, current.get());
                        for (size_t i = 0; i < found.number; ++i) {
                            sums[found.index[i]] += found.value[i] * mult;
                        }
                    } else {
                        auto found = get_dense<column>(x, vbuffer, current.get());
                        for (size_t i = 0; i < pos.second; ++i) {
                            sums[pos.first + i] += found[i] * mult;
                        }
                    }
                }

                std::copy_n(sums.begin() + pos.first, pos.second, output.begin() + pos.first);
            }, t);
        }

        for (int t = 0; t < nthreads; ++t) {
            jobs[t].join();
        }
    }

    template<bool column, bool sparse, class Right>
    void running(const Right& rhs, Eigen::VectorXd& output, Workspace& work) const {
        if (nthreads == 1) {
            running_serial<column, sparse>(rhs, output, work);
        } else {
            running_parallel<column, sparse>(rhs, output, work);
        }
    }

private:
    template<bool column, bool sparse, class Right, class TatamiWork>
    double direct_raw(const Right& rhs, size_t x, double* vbuffer, int* ibuffer, TatamiWork* current) const {
        double total = 0;

        if constexpr(sparse) {
            auto found = get_sparse_range<column>(x, vbuffer, ibuffer, current);
            for (size_t i = 0; i < found.number; ++i) {
                total += found.value[i] * rhs.coeff(found.index[i]);
            }
        } else {
            size_t other = rhs.size();
            auto found = get_dense<column>(x, vbuffer, current);
            for (size_t i = 0; i < other; ++i) {
                total += found[i] * rhs.coeff(i);
            }
        }

        return total;
    }

    template<bool column, bool sparse, class Right>
    void direct_serial(const Right& rhs, Eigen::VectorXd& output, Workspace& work) const {
        auto& current = get_work<column>(work).front();
        auto vbuffer = work.values.front().data();
        auto ibuffer = work.indices.front().data();

        output.setZero();
        size_t end = get_dimension<column>();
        for (size_t x = 0; x < end; ++x) {
            output.coeffRef(x) = direct_raw<column, sparse>(rhs, x, vbuffer, ibuffer, current.get());
        }
    }

    template<bool column, bool sparse, class Right>
    void direct_parallel(const Right& rhs, Eigen::VectorXd& output, Workspace& work) const {
        auto& jobs = work.jobs;
        auto& currentwork = get_work<column>(work);

        output.setZero();
        size_t start = 0;
        size_t ntasks = get_dimension<column>();
        size_t per_thread = std::ceil(static_cast<double>(ntasks) / nthreads);
        jobs.clear();

        for (int t = 0; t < nthreads; ++t) {
            auto end = std::min(ntasks, start + per_thread);

            jobs.emplace_back([&](int thread, size_t s, size_t e) -> void {
                auto& current = currentwork[thread];
                auto vbuffer = work.values[thread].data();
                auto ibuffer = work.indices[thread].data();
                for (size_t x = s; x < e; ++x) {
                    output.coeffRef(x) = direct_raw<column, sparse>(rhs, x, vbuffer, ibuffer, current.get());
                }
            }, t, start, end);

            start = end;
        }

        for (int t = 0; t < nthreads; ++t) {
            jobs[t].join();
        }
    }

    template<bool column, bool sparse, class Right>
    void direct(const Right& rhs, Eigen::VectorXd& output, Workspace& work) const {
        if (nthreads == 1) {
            direct_serial<column, sparse>(rhs, output, work);
        } else {
            direct_parallel<column, sparse>(rhs, output, work);
        }
    }

public:
    template<class Right>
    void multiply(const Right& rhs, Workspace& work, Eigen::VectorXd& output) const {
        if (mat->prefer_rows()) {
            if (mat->sparse()) {
                direct<false, true>(rhs, output, work);
            } else {
                direct<false, false>(rhs, output, work);
            }
        } else {
            if (mat->sparse()) {
                running<true, true>(rhs, output, work);
            } else {
                running<true, false>(rhs, output, work);
            }
        }
    }

    template<class Right>
    void adjoint_multiply(const Right& rhs, Workspace& work, Eigen::VectorXd& output) const {
        if (mat->prefer_rows()) {
            if (mat->sparse()) {
                running<false, true>(rhs, output, work);
            } else {
                running<false, false>(rhs, output, work);
            }
        } else {
            if (mat->sparse()) {
                direct<true, true>(rhs, output, work);
            } else {
                direct<true, false>(rhs, output, work);
            }
        }
    }

    Eigen::MatrixXd realize() const {
        Eigen::MatrixXd output(mat->nrow(), mat->ncol());
        tatami::convert_to_dense<false>(mat, output.data());
        return output;
    }
};

//[[Rcpp::export(rng=false)]]
Rcpp::List irlba_tatami(SEXP input, int rank, int nthreads, int seed) {
    auto shared = extract_NumericMatrix_shared(input);
    size_t NR = shared->nrow();
    size_t NC = shared->ncol();

    auto trans = tatami::make_DelayedTranspose(std::move(shared));
    TatamiWrapper mat(trans.get(), nthreads);

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
