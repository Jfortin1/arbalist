# Tests the saveRegionMatrix function.
# library(testthat); library(arbalist); source("test-saveRegionMatrix.R")

library(GenomicRanges)
test_that("region sanitization works as expected for GRanges", {
    ref <- GRanges(c("chrA:1-100:+", "chrA:50-200:-", "chrB:20-50:*", "chrB:51-60:+", "chrB:45-55:-", "chrC:100-200:*", "chrC:1-50:+"))
    sanitized <- arbalist:::sanitize_regions(ref)

    expect_identical(sanitized$ids$chrA, c(0L, 1L))
    expect_identical(sanitized$starts$chrA, c(0L, 100L))
    expect_identical(sanitized$ends$chrA, c(49L, 200L))

    expect_identical(sanitized$ids$chrB, c(2L, 3L))
    expect_identical(sanitized$starts$chrB, c(19L, 55L))
    expect_identical(sanitized$ends$chrB, c(44L, 60L))

    expect_identical(sanitized$ids$chrC, c(6L, 5L))
    expect_identical(sanitized$starts$chrC, c(0L, 99L))
    expect_identical(sanitized$ends$chrC, c(50L, 200L))
})

library(GenomicRanges)
test_that("region sanitization works as expected for GRangesLists", {
    ref <- GRanges(c("chrA:1-100:+", "chrA:50-200:-", "chrB:20-50:*", "chrB:51-60:+", "chrB:45-55:-", "chrC:100-200:*", "chrC:1-50:+"))

    # simple partitioning.
    {
        gr1 <- splitAsList(ref, c(1,1,2,2,2,3,3))
        sanitized <- arbalist:::sanitize_regions(gr1)

        expect_identical(sanitized$ids$chrA, 0L)
        expect_identical(sanitized$starts$chrA, 0L)
        expect_identical(sanitized$ends$chrA, 200L)

        expect_identical(sanitized$ids$chrB, 1L)
        expect_identical(sanitized$starts$chrB, 19L)
        expect_identical(sanitized$ends$chrB, 60L)

        expect_identical(sanitized$ids$chrC, c(2L, 2L))
        expect_identical(sanitized$starts$chrC, c(0L, 99L))
        expect_identical(sanitized$ends$chrC, c(50L, 200L))
    }

    # Complex partitioning.
    {
        gr2 <- splitAsList(ref, c(1,3,2,3,2,3,1))
        sanitized <- arbalist:::sanitize_regions(gr2)

        expect_identical(sanitized$ids$chrA, c(0L, 2L))
        expect_identical(sanitized$starts$chrA, c(0L, 100L))
        expect_identical(sanitized$ends$chrA, c(49L, 200L))

        expect_identical(sanitized$ids$chrB, c(1L, 2L))
        expect_identical(sanitized$starts$chrB, c(19L, 55L))
        expect_identical(sanitized$ends$chrB, c(50L, 60L))

        expect_identical(sanitized$ids$chrC, c(0L, 2L))
        expect_identical(sanitized$starts$chrC, c(0L, 99L))
        expect_identical(sanitized$ends$chrC, c(50L, 200L))
    }
})

reference_counter <- function(fname, regions, allowed.cells) {
    info <- read.delim(fname, header=FALSE, sep="\t", comment.char="#")
    if (!is.null(allowed.cells)) {
        info <- info[info[,4] %in% allowed.cells,,drop=FALSE]
    } else {
        allowed.cells <- info[,4]
        allowed.cells <- allowed.cells[!duplicated(allowed.cells)]
    }

    sanitized <- arbalist:::sanitize_regions(regions, decompose = FALSE)
    starts <- GRanges(info[,1], IRanges(info[,2], info[,2]))
    start.overlap <- findOverlaps(starts, sanitized$regions, select="first")
    ends <- GRanges(info[,1], IRanges(info[,3], info[,3]))
    end.overlap <- findOverlaps(ends, sanitized$regions, select="first")

    keep <- is.na(end.overlap) != is.na(start.overlap) | end.overlap == start.overlap
    keep[is.na(keep)] <- FALSE
    start.overlap <- start.overlap[keep]
    end.overlap <- end.overlap[keep]
    fid <- ifelse(is.na(start.overlap), end.overlap, start.overlap)

    cid <- match(info[keep,4], allowed.cells)
    out <- Matrix::sparseMatrix(i=fid, j=cid, x=rep(1, length(cid)), dims=c(length(regions), length(allowed.cells)))
    colnames(out) <- allowed.cells

    return(out)
}

mockRegions <- function(num.fragments, seq.lengths, width.range, num.groups = NULL) {
    seq <- sample(names(seq.lengths), num.fragments, replace=TRUE)
    limits <- (seq.lengths - 1L)[seq] # avoid overlapping the end position.
    starts <- floor(runif(num.fragments) * limits)
    ends <- pmin(starts + floor(runif(num.fragments, width.range[1], width.range[2])), limits)

    regions <- sort(GRanges(seq, IRanges(starts, ends)))
    if (!is.null(num.groups)) {
        g <- sample(num.groups, num.fragments, replace=TRUE)
        regions <- splitAsList(regions, g)
    }

    regions
}

library(HDF5Array)
test_that("saveTileMatrix compares correctly to the reference with GRanges and no restriction", {
    seq.lengths <- c(chrA = 10000, chrB = 100000, chrC = 1000)
    temp <- tempfile(fileext = ".gz")
    mockFragmentFile(temp, seq.lengths, 1e3, cell.names = LETTERS)

    regions <- mockRegions(200, seq.lengths, c(10, 50))
    ref <- reference_counter(temp, regions, allowed.cells = LETTERS)

    temp.h5 <- tempfile(fileext = ".h5")
    saveRegionMatrix(temp, regions, output.file=temp.h5, output.name="WHEE", barcodes = LETTERS)
    obs <- H5SparseMatrix(temp.h5, "WHEE")
    colnames(obs) <- LETTERS

    expect_identical(as(obs, "dgCMatrix"), ref)
})
