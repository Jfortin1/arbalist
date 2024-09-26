#ifndef TATAMI_R_PARALLELIZE_HPP
#define TATAMI_R_PARALLELIZE_HPP

/**
 * @cond
 */
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN
/**
 * @endcond
 */

#include "manticore/manticore.hpp"
#include <thread>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

/**
 * @file parallelize.hpp
 *
 * @brief Safely parallelize for unknown matrices.
 */

namespace tatami_r {

/**
 * Retrieve a global `manticore::Executor` object for all **tatami_r** uses.
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 *
 * @return Reference to a global `manticore::Executor`.
 */
inline manticore::Executor& executor() {
    // This should end up resolving to a single instance, even across dynamically linked libraries:
    // https://stackoverflow.com/questions/52851239/local-static-variable-linkage-in-a-template-class-static-member-function
    static manticore::Executor mexec;
    return mexec;
}

/**
 * @tparam Function_ Function to be executed.
 * @tparam Index_ Integer type for the task indices.
 *
 * @param fun Function to run in each thread.
 * This is a lambda that should accept three arguments:
 * - Integer containing the thread ID.
 * - Integer specifying the index of the first task to be executed in a thread.
 * - Integer specifying the number of tasks to be executed in a thread.
 * @param ntasks Number of tasks to be executed.
 * @param nthreads Number of threads to parallelize over.
 *
 * This function is a drop-in replacement for `tatami::parallelize()`.
 * The series of integers from `[0, ntasks)` is split into `nthreads` contiguous ranges.
 * Each range is used as input to a call to `fun` within a thread created by the standard `<thread>` library. 
 * Serialization can be achieved via `<mutex>` in most cases, or `manticore::Executor::run()` if the task must be performed on the main thread (see `executor()`).
 *
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 */ 
template<class Function_, class Index_>
void parallelize(Function_ fun, Index_ ntasks, int nthreads) {
    if (ntasks == 0) {
        return;
    }

    if (nthreads <= 1 || ntasks == 1) {
        fun(0, 0, ntasks);
        return;
    }

    Index_ tasks_per_worker = ntasks / nthreads;
    int remainder = ntasks % nthreads;
    if (tasks_per_worker == 0) {
        tasks_per_worker = 1; 
        remainder = 0;
        nthreads = ntasks;
    }

    auto& mexec = executor();
    mexec.initialize(nthreads, "failed to execute R command");

    std::vector<std::thread> runners;
    runners.reserve(nthreads);
    std::vector<std::exception_ptr> errors(nthreads);

    Index_ start = 0;
    for (int w = 0; w < nthreads; ++w) {
        Index_ length = tasks_per_worker + (w < remainder);

        runners.emplace_back([&](int id, Index_ s, Index_ l) {
            try {
                fun(id, s, l);
            } catch (...) {
                errors[id] = std::current_exception();
            }
            mexec.finish_thread();
        }, w, start, length);

        start += length;
    }

    mexec.listen();
    for (auto& x : runners) {
        x.join();
    }

    for (const auto& err : errors) {
        if (err) {
            std::rethrow_exception(err);
        }
    }
}

}

/**
 * @cond
 */
#endif
/**
 * @endcond
 */

#endif
