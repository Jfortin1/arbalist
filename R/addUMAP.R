#' Add UMAP transformation from MAE
#' 
#' Creating UMAP embedding from iterative LSI reduced dimensions contained with a MAE and adding the UMAP result back to the MAE
#' 
#' @param mae MultiAssayExperiment
#' @param name.iterative.lsi A string specifying the name of the reduced dimensions to create the UMAP from
#' @param name.umap A string specifying the name of the new UMAP reduced dimensions
#' @param num.neighbors An integer describing the number of nearest neighbors to compute a UMAP. This argument is passed to n_neighbors in uwot::umap().
#' @param min.dist A number specifying effective minimum distance between embedded points. See [uwot::umap]
#' @param metric A string specifying metric for calculating distance. ex. "cosine", "euclidean", etc. Multiple metrics can be specified. See [uwot::umap] for details.
#' @param seed A number used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param sample.cells A integer specifying the number of cells to subsample for UMAP creation. The cells not subsampled will be re-projected using [uwot::umap_transform] to the Embedding. Only recommended for extremely large number of cells.
#' @param num.threads.transform A integer specifying the number of threads for [uwot::umap_transform]
#' @param outlier.quantile A numeric (0 to 1) describing the distance quantile in the subsampled cels (see `sample.cells`) to use to filter poor quality re-projections.
#' @param ... Additional parameters to pass to [uwot::umap()]
#'
#' @return A \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#'
#' @export
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SingleCellExperiment reducedDims altExp reducedDim<- reducedDim altExpNames altExp<- reducedDimNames
addUMAP <- function(
  mae,
  name.iterative.lsi = 'iterativeLSI',
  name.umap = paste0(name.iterative.lsi,'_UMAP'),
  num.neighbors = 40,
  min.dist = 0.4,
  metric = "cosine",
  seed = 1,
  sample.cells = NULL,
  num.threads.transform = 1,
  outlier.quantile = 0.9,
  ...
) {
  
  # Find the IterativeLSI result
  lsi.matrix <- NULL;
  lsi.exp.idx <- NULL;
  lsi.alt.exp.name <- NULL;
  for(exp.idx in seq_len(length(mae))) {
    if(name.iterative.lsi %in% reducedDimNames(mae[[exp.idx]])) {
      lsi.matrix <- reducedDim(mae[[exp.idx]], name.iterative.lsi)
      # record where the iterative LSI was saved so we can put the UMAP in the same place
      lsi.exp.idx <- exp.idx
      break
    }
    if(is.null(lsi.matrix)) {
      for(name.alt.exp in altExpNames(mae[[1]])) {
        if(name.iterative.lsi %in% reducedDimNames(altExp(mae[[exp.idx]], name.alt.exp))) {
          lsi.matrix <- reducedDim(altExp(mae[[exp.idx]], name.alt.exp), name.iterative.lsi)
          # record where the iterative LSI was saved so we can put the UMAP in the same place
          lsi.exp.idx <- exp.idx
          lsi.alt.exp.name <- name.alt.exp
          break
        }
      }
    }
  }
  if(is.null(lsi.matrix)) {
    stop(paste0(name.iterative.lsi,' is not found in mae'))
  }
  
  # Create the UMAP embedding
  umap.res <- getUMAPFromIterativeLSI(
    lsi.matrix,
    num.neighbors = num.neighbors,
    min.dist = min.dist,
    metric = metric,
    seed = seed,
    sample.cells = sample.cells,
    num.threads.transform = num.threads.transform,
    outlier.quantile = outlier.quantile,
    ...
  )
  
  # save the new reduced dimensions to the MAE
  if(is.null(lsi.alt.exp.name)) {
    reducedDim(mae[[exp.idx]], name.umap) <- umap.res
  } else {
    reducedDim(altExp(mae[[exp.idx]], lsi.alt.exp.name), name.umap) <- umap.res
  }
  
  return(mae)
}