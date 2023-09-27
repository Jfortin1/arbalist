#ifndef TATAMI_STATS_UTILS_HPP
#define TATAMI_STATS_UTILS_HPP

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "../utils/Oracles.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

#ifndef TATAMI_CUSTOM_PARALLEL
#ifndef _OPENMP
#include <thread>
#endif
#include <string>
#include <stdexcept>
#endif

/**
 * @file utils.hpp
 *
 * @brief Utilities for computing matrix statistics.
 */

namespace tatami {

/**
 * Apply a function to a set of tasks in parallel.
 * This can be done using:
 *
 * - OpenMP, if available and enabled by the compiler.
 * - Using a custom parallelization scheme, by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro. 
 *   This should accept the `fun`, `tasks` and `threads` arguments as below.
 * - `<thread>`, otherwise.
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current thread.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `thread`, the thread number executing this task range.
 *   This will be passed as a `size_t`.
 * - `task_start`, the start index of the task range.
 *   This will be passed as a `Index_`.
 * - `task_length`, the number of tasks in the task range.
 *   This will be passed as a `Index_`.
 * @tparam Index_ Integer type for the number of tasks.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * @param tasks Number of tasks.
 * @param threads Number of threads.
 */
template<bool parallel_ = true, class Function_, typename Index_>
void parallelize(Function_ fun, Index_ tasks, size_t threads) {
    if constexpr(parallel_) {
        if (threads > 1) {
#ifndef TATAMI_CUSTOM_PARALLEL
            Index_ worker_size = (tasks / threads) + (tasks % threads > 0); // Ceiling of an integer division.
            threads = (tasks / worker_size) + (tasks % worker_size > 0); // Set the actual number of required threads.
            std::vector<std::string> errors(threads);

#if defined(_OPENMP)
            #pragma omp parallel for num_threads(threads)
            for (size_t t = 0; t < threads; ++t) {
                Index_ start = worker_size * t; // Will not overflow due to the above recomputation of 'threads'.
                Index_ remaining = tasks - start; // Must be positive, as otherwise 'tasks % worker_size = 0' and the iteration wouldn't even get here.

                try {
                    fun(t, start, std::min(remaining, worker_size)); // Use 'remaining' to avoid potential overflow from computing 'end = start + worker_size'.
                } catch (std::exception& e) {
                    errors[t] = e.what();
                } catch (...) {
                    errors[t] = "unknown error in thread " + std::to_string(t);
                }
            }

#else
            Index_ first = 0;
            std::vector<std::thread> workers;
            workers.reserve(threads);

            for (size_t t = 0; t < threads && first < tasks; ++t) {
                Index_ remaining = tasks - first;
                Index_ len = std::min(remaining, worker_size);
                workers.emplace_back([&fun,&errors](int t, Index_ first, Index_ len) -> void {
                    try {
                        fun(t, first, len);
                    } catch (std::exception& e) {
                        errors[t] = e.what();
                    } catch (...) {
                        errors[t] = "unknown error in thread " + std::to_string(t);
                    }
                }, t, first, len);
                first += len;
            }

            for (auto& wrk : workers) {
                wrk.join();
            }
#endif

            for (const auto& e : errors) {
                if (!e.empty()) {
                    throw std::runtime_error(e);
                }
            }

#else
            TATAMI_CUSTOM_PARALLEL(std::move(fun), tasks, threads);
#endif
            return;
        }
    }

    fun(0, 0, tasks);
    return;
}

/**
 * @tparam row_ Whether to perform extraction on rows.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column index.
 * @tparam Args_ Types of further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @param mat Matrix to iterate over.
 * @param iter_start Index of the first row/column of the iteration range.
 * @param iter_length Number of rows/columns in the iteration range.
 * @param args Further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @return An `Extractor` object for iteration over consecutive rows/columns in `[iter_start, iter_start + iter_length)`.
 *
 * This function is equivalent to `new_extractor()` but additionally calls `Extractor::set_oracle()` with a `ConsecutiveOracle` instance.
 * `Matrix` implementations that are oracle-aware can then perform pre-fetching of future accesses for greater performance.
 * Of course, this assumes that the iteration over the target dimension does actually involve consecutive elements from `iter_start` to `iter_start + iter_length`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>* mat, Index_ iter_start, Index_ iter_length, Args_&&... args) {
    auto ext = new_extractor<row_, sparse_>(mat, std::forward<Args_>(args)...);
    if (mat->uses_oracle(row_)) {
        ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(iter_start, iter_length));
    }
    return ext;
}

namespace stats {

/**
 * Count the total number of groups, typically for per-group memory allocations.
 *
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Index_ Integer type for the number of observations.
 *
 * @param[in] group Pointer to an array of group assignments per observation.
 * Each assignment should be an integer in `[0, G)` where `G` is the total number of groups.
 * @param n Number of observations, i.e., the length of the array referenced by `group`.
 *
 * @return Total number of groups, i.e., `G`.
 * Note that not all groups may actually have non-zero occurrences in `group`.
 */
template<typename Group_, typename Size_>
size_t total_groups(const Group_* group, Size_ n) {
    if (n) {
        return static_cast<size_t>(*std::max_element(group, group + n)) + 1;
    } else {
        return 0;
    }
}

/**
 * Count the occurrences of each group.
 *
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Index_ Integer type for the number of observations.
 *
 * @param[in] group Pointer to an array of group assignments per observation.
 * Each assignment should be an integer in `[0, G)` where `G` is the total number of groups.
 * @param n Number of observations, i.e., the length of the array referenced by `group`.
 *
 * @return Vector of length equal to `G`, containing the number of occurrences of each group.
 */
template<typename Group_, typename Size_>
std::vector<Size_> tabulate_groups(const Group_* group, Size_ n) {
    auto ngroups = total_groups(group, n);
    std::vector<Size_> group_sizes(ngroups);
    for (Size_ r = 0; r < n; ++r) {
        ++(group_sizes[group[r]]);
    }
    return group_sizes;
}

}

}

#endif
