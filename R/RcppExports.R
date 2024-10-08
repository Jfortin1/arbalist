# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

count_fragment_size_distributions <- function(fragment_file) {
    .Call('_arbalist_count_fragment_size_distributions', PACKAGE = 'arbalist', fragment_file)
}

fragments_to_regions <- function(fragment_file, output_file, output_group, seqnames, region_ids, region_starts, region_ends, cellnames, num_regions, deflate_level, chunk_dim) {
    .Call('_arbalist_fragments_to_regions', PACKAGE = 'arbalist', fragment_file, output_file, output_group, seqnames, region_ids, region_starts, region_ends, cellnames, num_regions, deflate_level, chunk_dim)
}

fragments_to_tiles <- function(fragment_file, tile_size, output_file, output_group, seqlengths, seqnames, cellnames, deflate_level, chunk_dim) {
    .Call('_arbalist_fragments_to_tiles', PACKAGE = 'arbalist', fragment_file, tile_size, output_file, output_group, seqlengths, seqnames, cellnames, deflate_level, chunk_dim)
}

lsi_matrix_stats <- function(mat, nthreads) {
    .Call('_arbalist_lsi_matrix_stats', PACKAGE = 'arbalist', mat, nthreads)
}

create_pseudobulk_file <- function(fragment_files, output_file, cellnames) {
    invisible(.Call('_arbalist_create_pseudobulk_file', PACKAGE = 'arbalist', fragment_files, output_file, cellnames))
}

