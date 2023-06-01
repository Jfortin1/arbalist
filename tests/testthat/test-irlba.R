# This tests the various IRLBA utilities.
# library(testthat); library(arbalist); source("test-irlba.R")

library(testthat); library(arbalist); 

set.seed(1000)
x <- Matrix::rsparsematrix(1000, 500, 0.1)
y <- round(abs(x)*10)
ptr <- beachmat::initializeCpp(y)

yr <- Matrix::t(y)
r <- new("dgRMatrix", x=yr@x, j=yr@i, p=yr@p, Dim=dim(y))
rptr <- beachmat::initializeCpp(r)

test_that("IRLBA works for a realized matrix", {
    serial <- arbalist:::irlba_realized(ptr, rank=10, nthreads=1, seed=42)
    parallel <- arbalist:::irlba_realized(ptr, rank=10, nthreads=2, seed=42)
    expect_identical(serial, parallel)
    
    rserial <- arbalist:::irlba_realized(rptr, rank=10, nthreads=1, seed=42)
    expect_identical(rserial, serial)

    rparallel <- arbalist:::irlba_realized(rptr, rank=10, nthreads=2, seed=42)
    expect_identical(rserial, rparallel)
})

test_that("IRLBA works for a tatami matrix", {
    ref <- arbalist:::irlba_realized(ptr, rank=10, nthreads=1, seed=42)

    serial <- arbalist:::irlba_tatami(ptr, rank=10, nthreads=2, seed=42)
    parallel <- arbalist:::irlba_tatami(ptr, rank=10, nthreads=2, seed=42)
    expect_identical(serial, parallel)
    
    rserial <- arbalist:::irlba_tatami(rptr, rank=10, nthreads=1, seed=42)
    expect_identical(rserial, serial)

    rparallel <- arbalist:::irlba_tatami(rptr, rank=10, nthreads=2, seed=42)
    expect_identical(rserial, rparallel)
})

