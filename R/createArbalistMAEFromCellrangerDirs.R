#' Import a single cell ATAC-seq MultiAssayExperiment
#' 
#' Import results from scATAC-seq or multiome Cell Ranger results directories into a MultiAssayExperiment.
#' 
#' @param cellranger.res.dirs Vector of strings specifying a Cell Ranger scATAC-seq or multiome results directory. Vector names need to be sample names.
#' @inheritParams createArbalistMAE
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
    filter.rna.features.without.intervals = TRUE,
    BPPARAM = bpparam()
) {
  
  # Update the cellranger paths in case they don't point directly at the results
  cellranger.res.dirs <- .updateCellRangerPath(cellranger.res.dirs)

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
    use.alt.exp = use.alt.exp,
    main.exp.name = main.exp.name,
    filter.rna.features.without.intervals = filter.rna.features.without.intervals,
    BPPARAM = BPPARAM
  )
  
  return(mae)
}

.getFilesFromResDirs <- function(res.dirs, file.name.options) {
  selected.files <- unlist(sapply(
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

.updateCellRangerPath <- function(res.dirs) {
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