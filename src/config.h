#ifndef CONFIG_H
#define CONFIG_H

#include "Rcpp.h"

#include <thread>
#include <cmath>
#include <vector>

// Must be included before any knncolle Annoy include,
// so as to get rid of stderr usage in R libraries.
#define __ERROR_PRINTER_OVERRIDE__  REprintf 

template<class Function>
void run_parallel(int total, Function fun, int nthreads) {
    if (nthreads == 1) {
        fun(0, total);
        return;
    }

    int jobs_per_worker = std::ceil(static_cast<double>(total)/nthreads);
    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    int first = 0;

    for (int w = 0; w < nthreads && first < total; ++w, first += jobs_per_worker) {
        int last = std::min(first + jobs_per_worker, total);
        workers.emplace_back(fun, first, last);
    }

    for (auto& wrk : workers) {
        wrk.join();
    }
}

template<class Function>
void run_parallel2(int nthreads, Function fun) {
    if (nthreads == 1) {
        fun(0);
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(nthreads);

    for (int w = 0; w < nthreads; ++w) {
        workers.emplace_back(fun, w);
    }

    for (auto& wrk : workers) {
        wrk.join();
    }
}

#define ARBALIST_CUSTOM_PARALLEL run_parallel
#define IRLBA_CUSTOM_PARALLEL run_parallel2

#endif
