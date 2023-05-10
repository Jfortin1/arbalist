#' Create a mock fragment file
#'
#' Mock up a fragment file for examples and testing.
#'
#' @param output.file String containing a path to an output file.
#' @param seq.lengths Named integer vector containing the lengths of the reference sequences used for alignment.
#' Vector names should correspond to the names of the sequences. 
#' @param num.fragments Integer scalar, the average number of fragments per cell.
#' @param cell.names Character vector containing the cell names.
#' The length of this vector is used as the total number of cells.
#' @param width.range Integer vector of length 2, containing the range of possible fragment widths in base pairs.
#' @param read.range Integer vector of length 2, containing the range of the read count supporting each fragment.
#' @param comments Character vector of comments to be added to the start of the file.
#' @param compressed Logical scalar indicating whether the output file should be compressed. 
#'
#' @return A fragment file is created at \code{output.file}.
#' \code{NULL} is invisibly returned.
#'
#' @examples
#' temp <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp, c(chrA=1000, chrB=200000, chrC=200), 
#'     num.fragments=100, cell.names=LETTERS)
#'
#' X <- read.table(temp)
#' head(X)
#'
#' @author Aaron Lun
#'
#' @export
mockFragmentFile <- function(output.file, seq.lengths, num.fragments, cell.names, width.range = c(10, 1000), read.range = c(1, 10), comments = NULL, compressed = TRUE) {
    number <- num.fragments * length(cell.names)
    seq <- sample(names(seq.lengths), number, replace=TRUE)
    limits <- (seq.lengths - 1L)[seq] # avoid overlapping the end position.
    starts <- floor(runif(number) * limits)
    ends <- pmin(starts + floor(runif(number, width.range[1], width.range[2])), limits)

    o <- order(factor(seq, names(seq.lengths)), starts)
    df <- data.frame(seq, starts, ends, name=sample(cell.names, number, replace=TRUE), count=floor(runif(number, read.range[1], read.range[2])))
    df <- df[o,]

    handle <- if (compressed) gzfile(output.file, open="wb") else file(output.file, open="wb")
    on.exit(close(handle))

    if (length(comments)) {
        writeLines(con=handle, paste("# ", comments), sep="\n")
    }
    write.table(file=handle, df, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, eol="\n")

    invisible(NULL)
}
