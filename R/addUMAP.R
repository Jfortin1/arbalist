#' Add UMAP transformation from MAE
#' 
#' Creating UMAP embedding from iterative LSI reduced dimensions contained with a MAE and adding the UMAP result back to the MAE
#' 
#' @param mae MultiAssayExperiment
#' @param name.iterative.lsi A string specifying the name of the reduced dimensions to create the UMAP from. If there are multiple reduced dimensions in the MultiAssayExperiment with the same name, then specify the correct experiment name in the vector's names. This can also be a vector or strings if you want to create a UMAP based on the combination of reduced dimensions. The UMAP matrix will be added to the same location as we find the first reduced dimension.
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
  if(length(name.umap) > 1) {
    stop('Change name.umap to only be a single string.')
  }
  
  # retrieve the reduced dimension (probably iterative lsi) matrix from the MAE
  reduced.dim.list <- list()
  for(i in name.iterative.lsi) {
    reduced.dim.list[[i]] <- findReducedDimRes(mae,i)
  }
  
  reduced.dim.matrix <- reduced.dim.list[[1]]$matrix;
  if(length(reduced.dim.list) > 1) {
    for(i in 2:length(reduced.dim.list)) {
      if(nrow(reduced.dim.matrix) != nrow(reduced.dim.list[[i]]$matrix)) {
        stop('The reduced dimension matrices need to have the same number of rows.')
      }
      reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix)
    }
  }
  
  # Create the UMAP embedding
  umap.res <- getUMAPFromIterativeLSI(
    reduced.dim.matrix,
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
  if(is.null(reduced.dim.list[[1]]$alt.exp.name)) {
    reducedDim(mae[[reduced.dim.list[[1]]$exp.idx]], name.umap) <- umap.res
  } else {
    reducedDim(altExp(mae[[reduced.dim.list[[1]]$exp.idx]], reduced.dim.list[[1]]$alt.exp.name), name.umap) <- umap.res
  }
  
  return(mae)
}