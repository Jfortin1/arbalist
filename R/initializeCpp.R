#' Initialize the matrix in C++ memory space
#'
#' Initialize a tatami matrix object in C++ memory space from an in-memory or HDF5-backed sparse matrix, or any DelayedArray derivatives thereof. 
#' This object simply references the R memory space and avoids making any copies of its own, so it can be cheaply re-created when needed inside each function.
#'
#' @param x A sparse matrix-like object, typically from the \pkg{Matrix} or \pkg{HDF5Array} packages.
#' Alternatively, a \link{DelayedArray} wrapper around such objects.
#' @param ... Further arguments used by specific methods.
#' This includes:
#' \itemize{
#' \item \code{force.integer}: logical scalar indicating whether double-precision \code{x} should be forced into integers.
#' Only used for loading a HDF5-backed matrix into memory.
#' }
#'
#' @return An external pointer to a C++ object containing a tatami matrix.
#'
#' @details
#' Do not attempt to serialize this object; it contains a pointer to external memory, and will not be valid after a save/load cycle.
#' Users should not be exposed to the returned pointers; rather, each \pkg{arbalest} function should call \code{initializeCpp} at the start to obtain a C++ object for further processing.
#' As mentioned before, this initialization process is very cheap so there is no downside from just recreating the object within each function body.
#'
#' The HDF5-backed sparse matrix is fully loaded into C++ memory upon initialization, to avoid the performance hit from constantly querying the file system.
#' Obviously, this is a rather time-consuming process, so to maintain the cheapness of this function, we use a global cache to store the pointer returned by the first initialization of any HDF5-backed \code{x}.
#' Subsequent calls to this function will then re-use the cached pointer.
#' If the process is using too much memory, users can flush the cache by calling \code{flushFileBackedCache()}.
#' 
#' @examples
#' # Mocking up a count matrix:
#' x <- Matrix::rsparsematrix(1000, 100, 0.1)
#' y <- round(abs(x))
#'
#' stuff <- initializeCpp(y)
#' str(stuff)
#' 
#' @export
#' @aliases
#' initializeCpp,dgCMatrix-method
#' initializeCpp,dgRMatrix-method
#' initializeCpp,H5SparseMatrixSeed-method
#' flushFileBasedCache
#' @import methods
setGeneric("initializeCpp", function(x, ...) standardGeneric("initializeCpp"))

####################################################################################
####################################################################################

#' @export
#' @import Matrix
setMethod("initializeCpp", "dgCMatrix", function(x, ...) initialize_from_memory(x@x, x@i, x@p, nrow(x), ncol(x), byrow=FALSE))

#' @export
setMethod("initializeCpp", "dgRMatrix", function(x, ...) initialize_from_memory(x@x, x@j, x@p, nrow(x), ncol(x), byrow=TRUE))

####################################################################################
####################################################################################

#' @export
#' @importFrom HDF5Array H5SparseMatrixSeed
setMethod("initializeCpp", "H5SparseMatrixSeed", function(x, ..., force.integer=TRUE) {
    ikey <- as.character(force.integer)

    if (x@filepath %in% names(file.based.cache$cache)) {
        fcache <- file.based.cache$cache[[x@filepath]]
        if (x@group %in% names(fcache)) {
            gcache <- fcache[[x@group]]
            if (ikey %in% names(gcache)) {
                return(gcache[[ikey]])
            }
        }
    }

    ptr <- initialize_from_hdf5(
        x@filepath, 
        x@group, 
        nrow(x), 
        ncol(x), 
        byrow=!is(x, "CSC_H5SparseMatrixSeed"), 
        forced=force.integer
    )

    if (!(x@filepath %in% names(file.based.cache$cache))) {
        fcache <- file.based.cache$cache[[x@filepath]]
    } else {
        fache <- list()
    }

    if (x@group %in% names(fcache)) {
        gcache <- fcache[[x@group]]
    } else {
        gcache <- list()
    }

    gcache[[ikey]] <- ptr
    fcache[[x@group]] <- gcache
    file.based.cache$cache[[x@filepath]] <- fcache

    ptr
})

file.based.cache <- new.env()
file.based.cache$cache <- list()

#' @export
flushFileBasedCache <- function() {
    file.based.cache$cache <- list()
    invisible(NULL)
}

####################################################################################
####################################################################################

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedMatrix", function(x, ...) {
    initializeCpp(x@seed, ...)
})

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedAbind", function(x, ...) {
    collected <- lapply(x@seeds, initializeCpp, ...)
    apply_combine(collected, x@along == 1L)
})

#' @export
setMethod("initializeCpp", "DelayedAperm", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)
    if (x@perm[1] == 2L && x@parm[2] == 1L) {
        apply_transpose(seed)
    } else {
        seed
    }
})

#' @export
setMethod("initializeCpp", "DelayedSubset", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)

    for (i in seq_along(x@index)) {
        idx <- x@index[[i]]
        if (!is.null(idx)) {
            seed <- apply_subset(seed, idx, i == 1L)
        }
    }

    seed
})

#' @export
setMethod("initializeCpp", "DelayedSetDimnames", function(x, ...) {
    initializeCpp(x@seed, ...)
})

####################################################################################
####################################################################################

supported.Ops <- c("+", "-", "*", "/")

.apply_arithmetic <- function(seed, op, val, right, row) {
    if (op == "+") {
        return(apply_addition(seed, val, row))
    } else if (op == "*") {
        return(apply_multiplication(seed, val, row))
    } else if (op == "/") {
        return(apply_division(seed, val, right, row))
    } else if (op == "-") {
        return(apply_subtraction(seed, val, right, row))
    }
}

#' @export
setMethod("initializeCpp", "DelayedUnaryIsoOpWithArgs", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)

    # Figuring out the identity of the operation.
    chosen <- NULL
    for (p in supported.Ops) {
        if (identical(x@OP, get(p, envir=baseenv()))) {
            chosen <- p
            break
        }
    }
    if (is.null(chosen)) {
        stop("unknown operation in ", class(x))
    }

    # Saving the left and right args. There should only be one or the other.
    # as the presence of both is not commutative.
    if (length(x@Rargs) + length(x@Largs) !=1) {
        stop("'DelayedUnaryIsoApWithArgs' should operate on exactly one argument")
    }

    right <- length(x@Rargs) > 0
    if (right) {
        args <- x@Rargs[[1]]
        along <- x@Ralong[1]
    } else {
        args <- x@Largs[[1]]
        along <- x@Lalong[1]
    }
    row <- along == 1L

    .apply_arithmetic(seed, op, args, right, row)
})

####################################################################################
####################################################################################

.unary_Ops <- function(seed, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (generic %in% supported.Ops) {
        delayed_op <- NULL
        e1 <- envir$e1
        e2 <- envir$e2

        if (missing(e2)) {
            if (generic == "+") {
                return(seed)
            } else if (generic == "-") {
                return(apply_multiplication(seed, -1, TRUE))
            } else {
                stop("second argument can only be missing for unary '+' or '-'")
            }
        } else {
            right <- is(e1, "DelayedArray") # i.e., is the operation applied to the right of the seed?
            left <- is(e2, "DelayedArray") # i.e., is the operation applied to the left of the seed?

            write_string_scalar(file, name, "side", if (left) "left" else "right")
            val <- if (left) e1 else e2

            # Just guessing that it's applied along the rows, if it's not a scalar.
            return(.apply_arithmetic(seed, generic, val, right, TRUE))
        }
    }
}

#' @export
setMethod("initializeCpp", "DelayedUnaryIsoOpStack", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)
    for (i in rev(seq_along(x@OPS))) { # reverse order, as first operation is first applied (and thus needs to be closer to the leaf of the delayed tree).
        OP <- x@OPS[[i]]
        status <- FALSE 

        if (!status) {
            info <- .unary_Ops(seed, OP)
            if (status <- !is.null(info)) {
                seed <- info
            }
        } 

        if (!status) {
            stop("unknown OPS[[", i, "]] function in ", class(x))
        }
    }

    seed
})
