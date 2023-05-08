#' @export
#' @importFrom utils read.delim
#' @importFrom rhdf5 H5Fopen H5Fclose h5createFile h5createGroup h5createDataset H5Sunlimited
createTileMatrix <- function(fragment.files, output.path, output.name, seq.lengths=NULL, cell.names=NULL, tile.size = 500, max.count = 255, chunk.dim = 20000) {
    if (is.null(seq.lengths)) {
        # Obtaining the sequence lengths, if they weren't already available.
        out <- lapply(fragment.files, .process_fragment_header)
        everything <- unique(lapply(out, function(x) x$reference_path))
        if (length(everything) != 1) {
            stop("all fragment files should be generated with the same 'reference_path' if 'seq.lengths' is not supplied")
        }

        fai <- file.path(everything[[1]], "fasta", "genome.fa.fai")
        fai.info <- read.delim(fai, sep="\t", header=FALSE)[,1:2]
        seq.lengths <- fai.info[,2]
        names(seq.lengths) <- fai.info[,1]
    }

    if (!is.null(cell.names)) {
        if (length(cell.names) != length(fragment.files)) {
            stop("'cell.names' should have the same length as 'fragment.files'")
        }
    }

    if (max.count < 2^8) {
        dtype <- "H5T_NATIVE_UINT8"
    } else if (max.count < 2^16) {
        dtype <- "H5T_NATIVE_UINT16"
    } else {
        dtype <- "H5T_NATIVE_UINT32"
    }

    if (!file.exists(output.path)) {
        h5createFile(output.path)
    }

    (function() {
        handle <- H5Fopen(output.path)
        on.exit(H5Fclose(handle))
        h5createGroup(handle, output.name)

        h5createDataset(
            handle,
            paste0(output.name, "/data"),
            dims = 0,
            maxdims = H5Sunlimited(),
            H5type = dtype,
            chunk = chunk.dim 
        )

        h5createDataset(
            handle,
            paste0(output.name, "/indices"),
            dims    = 0,
            maxdims = H5Sunlimited(),
            H5type  = "H5T_NATIVE_UINT32",
            chunk   = chunk.dim
        )

        h5createDataset(
            handle,
            paste0(output.name, "/indptr"),
            dims    = 1,
            maxdims = H5Sunlimited(),
            H5type = "H5T_NATIVE_UINT64",
            fillValue = 0,
            chunk = chunk.dim 
        )
    })()

    collected <- vector("list", length(fragment.files))
    last <- 0
    for (i in seq_along(fragment.files)) {
        output <- tryCatch(
            dump_fragments_to_files(
                fragment_file = fragment.files[i], 
                tile_size = tile.size, 
                output_file = output.path, 
                output_group = output.name, 
                seqlengths = seq.lengths, 
                seqnames = names(seq.lengths), 
                cellnames = cell.names[[i]], 
                previous_nonzero = last
            ),
            error = function(e) stop("failed to parse fragment file '", fragment.files[i], "'\n- ", err$message)
        )
        last <- output[[1]]
        collected[[i]] <- output[[2]]
    }

    output <- NULL
    if (is.null(cell.names)) {
        output <- collected
        names(output) <- fragment.files
    }
    invisible(output)
}

.process_fragment_header <- function(file) {
    handle <- gzfile(file, open="rb")
    on.exit(close(handle))
    all.headers <- character(0)

    chunk <- 100
    repeat {
        lines <- readLines(handle, n = chunk)
        header <- startsWith(lines, "#")
        all.headers <- c(all.headers, sub("^# ", lines[header]))
        if (length(lines) < n || !all(header)) {
            break
        }
    }

    field <- sub("=.*", "", all.headers)
    value <- sub("[^=]+=", all.headers)
    split(value, field)
}
