#' Retrieves genome reference file name 
#' 
#' Helper function to extract genome reference file name from fragment file header.
#' 
#' @param fragment.file String specifying fragment file name
#' @return Named string vector
#' @author Natalie Fox
#' @examples
#' \dontrun{
#' download.file(
#' paste0("https://cf.10xgenomics.com/samples/cell-arc/1.0.0/",
#' "pbmc_granulocyte_sorted_10k/",
#' "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"), 
#' file.path(tempdir(), "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
#' )
#' fragment.file <- file.path(tempdir(), "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
#' extractGenomeRefFilenameFromFragmentFile(fragment.file)
#' }
#'
#' @export
extractGenomeRefFilenameFromFragmentFile <- function(fragment.file){
  handle <- gzfile(fragment.file, open = "rb")
  on.exit(close(handle))
  all.headers <- character(0)
  chunk <- 100
  repeat {
    lines <- readLines(handle, n = chunk)
    header <- startsWith(lines, "#")
    all.headers <- c(all.headers, sub("^# ", "", lines[header]))
    if (length(lines) < chunk || !all(header)) {
      break
    }
  }
  field <- sub("=.*", "", all.headers)
  value <- sub("[^=]+=", "", all.headers)
  info <- split(value, field)
  fai <- file.path(info$reference_path, "fasta", "genome.fa.fai")
  return(fai)
}