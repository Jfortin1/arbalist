# This tests the various LSI-related utilities.
# library(testthat); library(arbalist); source("test-iterativeLSI.R")

test_that("LSI-related utilities work as expected", {
    library(testthat); library(arbalist); 

    set.seed(1000)
    x <- Matrix::rsparsematrix(1000, 500, 0.1)
    y <- round(abs(x)*10)
    ptr <- initializeCpp(y)

    info <- arbalist:::lsi_matrix_stats(ptr, nthreads=1)
    expect_equal(info$sums, Matrix::colSums(y))
    expect_equal(info$frequency, Matrix::rowSums(y > 0))
    expect_identical(info, arbalist:::lsi_matrix_stats(ptr, nthreads=2))

    yr <- Matrix::t(y)
    r <- new("dgRMatrix", x=yr@x, j=yr@i, p=yr@p, Dim=dim(y))
    rptr <- initializeCpp(r)
    
    rinfo <- arbalist:::lsi_matrix_stats(rptr, nthreads=1)
    expect_identical(info, rinfo)
    expect_identical(rinfo, arbalist:::lsi_matrix_stats(rptr, nthreads=2))
})

