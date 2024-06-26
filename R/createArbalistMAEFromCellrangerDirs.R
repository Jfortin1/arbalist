#' Import a single cell ATAC-seq MultiAssayExperiment
#' 
#' Import results from scATAC-seq or multiome cellranger directories into a MultiAssayExperiment.
#' 
#' @param cellranger.res.dirs Vector of strings specifying a sunrise scATAC-seq or multiome cellranger results directory. Vector names need to be sample names.
#' @param output.dir String containing the directory where files should be output while creating the \linkS4class{MultiAssayExperiment}.
#' @param min.frags Number specifying the minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses. Cells containing greater than or equal to min.frags total fragments will be retained.
#' @param max.frags Number specifying the maximum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses. Cells containing less than or equal to max.frags total fragments will be retained.
#' @param tile.size Number specifying size of the tiles across the genome in base pairs.
#' @param seq.lengths Named integer vector containing the lengths of the reference sequences used for alignment.
#' @param gene.grs Genomic Ranges specifying gene coordinates for creating the gene score matrix. If NULL, then the geneset will be selected based on the genome version.
#' @param use.alt.exp Logical for selecting the MultiAssayExperiment structure. TRUE means that there will only be one experiment in the MultiAssayExperiment and all other experiments will be in alternative experiments. This option is only available if the columns are the same for all Matrices. FALSE means that each Matrix will be a separate experiment in the MAE.
#' @param main.exp.name String containing the name of the experiment that will be the main experiment when use.alt.exp is TRUE.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @return A \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#' @export
createArbalistMAEFromCellrangerDirs <- function(
    cellranger.res.dirs,
    output.dir = tempdir(),
    min.frags = 1000,
    max.frags = Inf,
    tile.size = 500,
    seq.lengths = NULL,
    gene.grs = NULL,
    use.alt.exp = FALSE,
    main.exp.name = 'TileMatrix500',
    BPPARAM = bpparam()
) {
  
  # Update the cellranger paths in case they don't point directly at the results
  cellranger.res.dirs <- .updateCellrangerPath(cellranger.res.dirs)

  # Create the MultiAssayExperiment from the list of Experiments
  mae <- createArbalistMAE(
    sample.names = names(cellranger.res.dirs),
    fragment.files = .getFilesFromResDirs(cellranger.res.dirs,c('atac_fragments.tsv.gz','fragments.tsv.gz')),
    filtered.feature.matrix.files = .getFilesFromResDirs(cellranger.res.dirs,c('filtered_feature_bc_matrix.h5','filtered_tf_bc_matrix.h5')),
    barcode.annotation.files = .getFilesFromResDirs(cellranger.res.dirs,c('per_barcode_metrics.csv','singlecell.csv')),
    sample.annotation.files = .getFilesFromResDirs(cellranger.res.dirs,'summary.csv'),
    output.dir = output.dir,
    min.frags = min.frags,
    max.frags = max.frags,
    tile.size = tile.size,
    seq.lengths = seq.lengths,
    gene.grs = gene.grs,
    use.alt.exp = use.alt.exp ,
    main.exp.name = main.exp.name,
    BPPARAM = BPPARAM
  )
  
  return(mae)
}

.getFilesFromResDirs <- function(res.dirs,file.name.options) {
  selected.files <- as.character(sapply(
    res.dirs,
    function(x,file.name.options){
      for(i in file.name.options) {
        potential.file <- paste0(x,'/',i)
        if(file.exists(potential.file)) {
          return(potential.file)
        } else {
          potential.files <- paste0(x,'/',list.files(x),'/',i)
          potential.file <- potential.files[file.exists(potential.files)]
          if(length(potential.file) == 1) {
            return(potential.file)
          }
        }
      }
    },
    file.name.options
  ))
  return(selected.files)
}

.updateCellrangerPath <- function(res.dirs) {
  # find the directory which directly contains the fragment file
  cleaned.res.dir <- as.character(sapply(
    res.dirs,
    function(x){
      if(file.exists(paste0(x,'/atac_fragments.tsv.gz'))) {
        x
      }
      else if(file.exists(paste0(x,'/fragments.tsv.gz'))) {
        x
      }
      else {
        potential.files <- paste0(x,'/',list.files(x),'/fragments.tsv.gz')
        paste0(x,'/',list.files(x))[file.exists(potential.files)]
      }
    }
  ))
  names(cleaned.res.dir) <- names(res.dirs)
  return(cleaned.res.dir)
}

.process_fragment_header <- function(file) {
  handle <- gzfile(file, open="rb")
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
  split(value, field)
}