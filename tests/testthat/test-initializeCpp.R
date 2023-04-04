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


