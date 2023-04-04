# This checks the initialization procedure.
# library(testthat); library(arbalist); source("test-initializeCpp.R")

am_i_ok <- function(ref, ptr, mutate=identity) {
    expect_identical(dim(ref), arbalist:::tatami_dim(ptr))
    for (i in seq_len(ncol(ref))) {
        expected <- mutate(ref[,i])
        expect_identical(expected, arbalist:::tatami_column(ptr, i))
    }
}

set.seed(1000)
x <- Matrix::rsparsematrix(1000, 100, 0.1)
y <- round(abs(x)*10)

test_that("initialization works correctly with a dgCMatrix", {
    ptr <- initializeCpp(y)
    am_i_ok(y, ptr)
})

test_that("initialization works correctly with a dgRMatrix", {
    z <- new("dgRMatrix", x=y@x, j=y@i, p=y@p, Dim=rev(y@Dim))
    ptr <- initializeCpp(z)
    am_i_ok(z, ptr)
})

library(HDF5Array)
library(rhdf5)
dump_hdf5 <- function(y, xtype="H5T_NATIVE_DOUBLE", itype="H5T_NATIVE_INT") {
    tmp <- tempfile(fileext=".h5")
    h5createFile(tmp)
    h5createGroup(tmp, "YAY")
    h5write(y@p, tmp, "YAY/indptr")
    h5write(dim(y), tmp, "YAY/shape")

    h5createDataset(tmp, "YAY/indices", dims=length(y@i), H5type=itype)
    h5write(y@i, tmp, "YAY/indices")

    h5createDataset(tmp, "YAY/data", dims=length(y@x), H5type=xtype)
    h5write(y@x, tmp, "YAY/data")

    TENxMatrix(filepath=tmp, group="YAY")
}

test_that("initialization works correctly with H5SparseMatrices", {
    A <- dump_hdf5(y)
    z <- initializeCpp(A)
    am_i_ok(A, z)

    # Auto-upgrates to a larger integer size.
    A <- dump_hdf5(y * 1e5)
    z <- initializeCpp(A)
    am_i_ok(A, z)

    # Respects doubles.
    A <- dump_hdf5(y * 1.5)
    z <- initializeCpp(A)
    am_i_ok(floor(A), z)

    z <- initializeCpp(A, force.integer=FALSE)
    am_i_ok(A, z)

    # Respects integer data types.
    A <- dump_hdf5(y, xtype="H5T_NATIVE_UINT16")
    z <- initializeCpp(A)
    am_i_ok(A * 1.0, z)

    # Force the use of a larger integer type for 'i'.
    A <- dump_hdf5(round(abs(Matrix::rsparsematrix(1e5, 10, 0.01) * 10)))
    z <- initializeCpp(A, no.sparse.copy=FALSE)
    am_i_ok(A, z)
})

test_that("initialization works correctly with DelayedArray", {
    z <- DelayedArray(y)
    ptr <- initializeCpp(z)
    am_i_ok(y, ptr)
})

test_that("initialization works correctly with DelayedArray transposition", {
    z0 <- DelayedArray(y)
    z <- t(z0)
    ptr <- initializeCpp(z)
    am_i_ok(t(y), ptr)
})

test_that("initialization works correctly with DelayedArray subsetting", {
    z0 <- DelayedArray(y)

    rkeep <- sample(nrow(y), 100)
    z <- z0[rkeep,]
    ptr <- initializeCpp(z)
    am_i_ok(y[rkeep,], ptr)

    ckeep <- sample(ncol(y), 10)
    z <- z0[,ckeep]
    ptr <- initializeCpp(z)
    am_i_ok(y[,ckeep], ptr)

    z <- z0[rkeep,ckeep]
    ptr <- initializeCpp(z)
    am_i_ok(y[rkeep,ckeep], ptr)

    rkeep <- 100:200
    ckeep <- 5:30
    z <- z0[rkeep,ckeep]
    ptr <- initializeCpp(z)
    am_i_ok(y[rkeep,ckeep], ptr)
})

test_that("initialization works correctly with DelayedArray combining", {
    z0 <- DelayedArray(y)

    x2 <- Matrix::rsparsematrix(99, 100, 0.1)
    y2 <- round(abs(x)*10)
    z <- rbind(z0, DelayedArray(y2))
    ptr <- initializeCpp(z)
    am_i_ok(rbind(y, y2), ptr)

    x2 <- Matrix::rsparsematrix(1000, 50, 0.1)
    y2 <- round(abs(x)*10)
    z <- cbind(z0, DelayedArray(y2))
    ptr <- initializeCpp(z)
    am_i_ok(cbind(y, y2), ptr)
})

test_that("initialization works correctly with scalar arithmetic", {
    z0 <- DelayedArray(y)

    z <- z0 + 1
    ptr <- initializeCpp(z)
    am_i_ok(y + 1, ptr)

    z <- z0 * 10
    ptr <- initializeCpp(z)
    am_i_ok(y * 10, ptr)

    z <- z0 - 1
    ptr <- initializeCpp(z)
    am_i_ok(y - 1, ptr)

    z <- 1 - z0
    ptr <- initializeCpp(z)
    am_i_ok(1 - y, ptr)

    z <- z0 / 10
    ptr <- initializeCpp(z)
    am_i_ok(y / 10, ptr)

    z <- 10 / z0
    ptr <- initializeCpp(z)
    am_i_ok(10 / y, ptr)
})

test_that("initialization works correctly with unary operations", {
    z0 <- DelayedArray(y)

    z <- +z0 
    ptr <- initializeCpp(z)
    am_i_ok(y, ptr)

    z <- -z0 
    ptr <- initializeCpp(z)
    am_i_ok(-y, ptr)
})

test_that("initialization works correctly with vector arithmetic", {
    z0 <- DelayedArray(y)

    {
        vr <- runif(nrow(y))

        z <- z0 + vr
        ptr <- initializeCpp(z)
        am_i_ok(y + vr, ptr)

        z <- z0 * vr
        ptr <- initializeCpp(z)
        am_i_ok(y * vr, ptr)

        z <- z0 - vr
        ptr <- initializeCpp(z)
        am_i_ok(y - vr, ptr)

        z <- vr - z0
        ptr <- initializeCpp(z)
        am_i_ok(vr - y, ptr)

        z <- z0 / vr
        ptr <- initializeCpp(z)
        am_i_ok(y / vr, ptr)

        z <- vr / z0
        ptr <- initializeCpp(z)
        am_i_ok(vr / y, ptr)
    }

    {
        vc <- runif(ncol(y))

        z <- sweep(z0, 2, vc, "+")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) + vc), ptr)

        z <- sweep(z0, 2, vc, "*")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) * vc), ptr)

        z <- sweep(z0, 2, vc, "-")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) - vc), ptr)

        z <- sweep(z0, 2, vc, "/")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) / vc), ptr)
    }
})


