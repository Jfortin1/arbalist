#' @export
createTileMatrix <- function(fragment.files, output.path, output.name, seq.lengths=NULL, tile.size = 500, max.count = 255, chunk.dim = 20000) {
    # Obtaining the sequence lengths, if they weren't already available.
    if (is.null(seq.lengths)) {
        out <- lapply(fragment.files, .process_fragment_header)
        everything <- unique(lapply(out, function(x) x$reference_path))
        if (length(everything) != 1) {
            stop("all fragment files should be generated with the same 'reference_path' if 'seq.lengths' is not supplied")
        }

        fai <- file.path(everything[[1]], "genome.fa.fai")
        fai.info <- read.delim(fai, sep="\t", header=FALSE)[,1:2]
        seq.lengths <- fai.info[,2]
        names(seq.lengths) <- fai.info[,1]
    }

    # Determining the maximimum count.
    if (max.index < 2^8) {
        dtype <- "H5T_NATIVE_UINT8"
    } else if (max.index < 2^16) {
        dtype <- "H5T_NATIVE_UINT16"
    } else {
        dtype <- "H5T_NATIVE_UINT32"
    }
    
    # Setting up the file.
    if (!file.exists(file)) {
        h5createFile(file)
    }

    (function() {
        handle <- H5Fopen(file)
        on.exit(H5Fclose(handle))
        h5createGroup(handle, name)

        h5createDataset(
            handle,
            file.path(name, "data"),
            dims = 0,
            maxdims = H5Sunlimited(),
            H5type = dtype,
            chunk = chunk.dim 
        )

        h5createDataset(
            handle,
            file.path(name, "indices"),
            dims    = 0,
            maxdims = H5Sunlimited(),
            H5type  = itype,
            chunk   = chunk.dim
        )
    })()

    for (x in fragment.files) {



    }
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

.fragments_to_tiles <- function(file) {
    contents <- read.delim(file, header=FALSE, sep="\t", colClasses=list("string", "integer", "integer", "string", NULL))
}


