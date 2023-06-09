#ifndef CONFIG_H
#define CONFIG_H

#include "Rtatami.h"

// Must be included before any knncolle Annoy include,
// so as to get rid of stderr usage in R libraries.
#define __ERROR_PRINTER_OVERRIDE__  REprintf 

// Correct handling of the main thread.
#define ARBALIST_CUSTOM_PARALLEL tatami_r::parallelize

// Special case when the number of threads == number of jobs.
template<class Function_>
void arbalist_run_parallel(int nthreads, Function_ fun) {
    if (nthreads == 1) {
        fun(0);
    } else {
        tatami_r::parallelize([&](int, size_t start, size_t len) {
            for (size_t i = 0; i < len; ++i) {
                fun(start + i);
            }
        }, nthreads, nthreads);
    }
}

#define IRLBA_CUSTOM_PARALLEL arbalist_run_parallel

#endif
