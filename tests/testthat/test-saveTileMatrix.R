# Tests the saveTileMatrix function.
# library(testthat); library(arbalist); source("test-saveTileMatrix.R")

library(HDF5Array)
reference_counter <- function(fname, seq.lengths, allowed.cells, tile.size) {
    info <- read.delim(fname, header=FALSE, sep="\t")
    if (!is.null(allowed.cells)) {
        info <- info[info[,4] %in% allowed.cells,,drop=FALSE]
    }

    positions <- ifelse(seq_len(nrow(info)) %% 2 == 1, info[,2], info[,3]) # alternate between start and end.

    nbins <- ceiling(seq.lengths / tile.size)
    offsets <- cumsum(c(0L, nbins))
    m <- match(info[,1], names(seq.lengths))
    bin.id <- offsets[m] + floor(positions / tile.size)

    if (is.null(allowed.cells)) {
        allowed.cells <- info[,4]
        allowed.cells <- allowed.cells[!duplicated(allowed.cells)]
    }

    cid <- match(info[,4], allowed.cells)
    out <- Matrix::sparseMatrix(i=bin.id + 1L, j=cid, x=rep(1, length(cid)), dims=c(sum(nbins), length(allowed.cells)))
    colnames(out) <- allowed.cells

    return(out)
}

test_that("saveTileMatrix compares correctly to the reference", {
    seq.lengths <- c(chrA = 10000, chrB = 100000, chrC = 1000)
    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS)

    ref <- reference_counter(temp, seq.lengths, allowed.cells = LETTERS, tile.size = 500)

    temp.h5 <- tempfile(fileext = ".h5")
    saveTileMatrix(temp, seq.lengths=seq.lengths, output.file=temp.h5, output.name="WHEE", barcodes = LETTERS)
    obs <- H5SparseMatrix(temp.h5, "WHEE")
    colnames(obs) <- LETTERS

    expect_identical(as(obs, "dgCMatrix"), ref)
})
