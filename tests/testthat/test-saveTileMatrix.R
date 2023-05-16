# Tests the saveTileMatrix function.
# library(testthat); library(arbalist); source("test-saveTileMatrix.R")

library(HDF5Array)
reference_counter <- function(fname, seq.lengths, allowed.cells, tile.size) {
    info <- read.delim(fname, header=FALSE, sep="\t", comment.char="#")
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

library(GenomeInfoDb)
test_that("saveTileMatrix compares correctly to the reference", {
    seq.lengths <- c(chrA = 10000, chrB = 100000, chrC = 1000)
    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS)
    ref <- reference_counter(temp, seq.lengths, allowed.cells = LETTERS, tile.size = 500)

    temp.h5 <- tempfile(fileext = ".h5")
    obs <- saveTileMatrix(temp, seq.lengths=seq.lengths, output.file=temp.h5, output.name="WHEE", barcodes = LETTERS)

    expect_identical(nrow(obs$counts), length(obs$tiles))
    expect_identical(as.character(runValue(seqnames(obs$tiles))), names(seq.lengths))
    expect_identical(as(obs$counts, "dgCMatrix"), ref)
})

test_that("saveTileMatrix handles comments correctly", {
    seq.lengths <- c(chrA = 501, chrB = 5001, chr1 = 50001)
    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS, comments=c("Hi", "I am a comment", "YAY"))
    ref <- reference_counter(temp, seq.lengths, allowed.cells = LETTERS, tile.size = 500)

    temp.h5 <- tempfile(fileext = ".h5")
    obs <- saveTileMatrix(temp, seq.lengths=seq.lengths, output.file=temp.h5, output.name="WHEE", barcodes = LETTERS)

    expect_identical(nrow(obs$counts), length(obs$tiles))
    expect_identical(as.character(runValue(seqnames(obs$tiles))), names(seq.lengths))
    expect_identical(as(obs$counts, "dgCMatrix"), ref)
})

test_that("saveTileMatrix works correctly with all cells", {
    seq.lengths <- c(chrA = 1501, chrD = 25001, chrC = 3501) # Scrambling the order to check that everyone respects non-alphanumeric sorting.
    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS)
    ref <- reference_counter(temp, seq.lengths, allowed.cells = NULL, tile.size = 500)

    temp.h5 <- tempfile(fileext = ".h5")
    obs <- saveTileMatrix(temp, seq.lengths=seq.lengths, output.file=temp.h5, output.name="WHEE", barcodes = NULL)

    expect_identical(nrow(obs$counts), length(obs$tiles))
    expect_identical(as.character(runValue(seqnames(obs$tiles))), names(seq.lengths))
    expect_identical(as(obs$counts, "dgCMatrix"), ref)
})

test_that("saveTileMatrix works correctly with restricted cells", {
    seq.lengths <- c(chrC = 1999, chrB = 999, chrA = 29999)
    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS)
    ref <- reference_counter(temp, seq.lengths, allowed.cells = LETTERS[1:5], tile.size = 500)

    temp.h5 <- tempfile(fileext = ".h5")
    obs <- saveTileMatrix(temp, seq.lengths=seq.lengths, output.file=temp.h5, output.name="WHEE", barcodes = LETTERS[1:5])

    expect_identical(nrow(obs$counts), length(obs$tiles))
    expect_identical(as.character(runValue(seqnames(obs$tiles))), names(seq.lengths))
    expect_identical(as(obs$counts, "dgCMatrix"), ref)
})

test_that("saveTileMatrix looks up the sequence lengths", {
    seq.lengths <- c(chrC = 1999, chrB = 999, chrA = 29999)
    temp.fai <- tempfile(fileext=".fai")
    write.table(data.frame(names(seq.lengths), seq.lengths), file=temp.fai, col.names=FALSE, row.names=FALSE)

    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS, comments=paste0("reference_path=", temp.fai))
    ref <- reference_counter(temp, seq.lengths, allowed.cells = NULL, tile.size = 500)

    temp.h5 <- tempfile(fileext = ".h5")
    obs <- saveTileMatrix(temp, seq.lengths=seq.lengths, output.file=temp.h5, output.name="WHEE", barcodes = NULL)

    expect_identical(nrow(obs$counts), length(obs$tiles))
    expect_identical(as.character(runValue(seqnames(obs$tiles))), names(seq.lengths))
    expect_identical(as(obs$counts, "dgCMatrix"), ref)
})
