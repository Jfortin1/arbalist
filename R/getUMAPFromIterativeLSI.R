#' Perform UMAP transformation
#' 
#' Creating UMAP matrix from the iterative LSI matrix
#' 
#' @param x a matrix containing iterativeLSI features by cells
#' @param num.neighbors An number describing the number of nearest neighbors to compute a UMAP. This argument is passed to n_neighbors in uwot::umap().
#' @param min.dist A number specifying effective minimum distance between embedded points. See [uwot::umap]
#' @param metric A string specifying metric for calculating distance. ex. "cosine", "euclidean", etc. Multiple metrics can be specified. See [uwot::umap] for details.
#' @param seed A number used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param sample.cells A integer specifying the number of cells to subsample for UMAP creation. The cells not subsampled will be re-projected using [uwot::umap_transform] to the Embedding. Only recommended for extremely large number of cells.
#' @param num.threads.transform A integer specifying the number of threads for [uwot::umap_transform]
#' @param outlier.quantile A numeric (0 to 1) describing the distance quantile in the subsampled cels (see `sample.cells`) to use to filter poor quality re-projections.
#' @param ... Additional parameters to pass to [uwot::umap()]
#' 
#' @return A data frame containing UMAP features by cells
#' 
#' @author Natalie Fox
#' @importFrom uwot umap umap_transform
#' @importFrom nabor knn
#' @importFrom stats quantile
#' @export
getUMAPFromIterativeLSI <- function(
  x,
  num.neighbors = 40,
  min.dist = 0.4,
  metric = "cosine",
  seed = 1,
  sample.cells = NULL,
  num.threads.transform = 1,
  outlier.quantile = 0.9,
  ...
) {

  # Provide an option for subsampling cells for creating the embedding.
  # This decreases run time and memory requirements but can lower the overall quality of the UMAP Embedding.
  estimate.umap <- FALSE
  if(!is.null(sample.cells) && sample.cells < nrow(x)) {
    selected.cells <- sample(seq_len(nrow(x)), sample.cells)
    cell.names <- rownames(x)
    x.cells.not.selected <- x[-selected.cells, , drop=FALSE]
    x <- x[selected.cells, , drop=FALSE]
    save.model <- TRUE
  }
  
  # Create UMAP embedding
  set.seed(seed)
  umap.res <- uwot::umap(
    x,
    n_neighbors = num.neighbors,
    min_dist = min.dist,
    metric = metric,
    ret_nn = estimate.umap,
    ret_model = estimate.umap,
    ...
  )

  umap.df <- data.frame(umap.res)
  colnames(umap.df) <- paste0("UMAP_Dimension_",seq_len(ncol(umap.df)))
  rownames(umap.df) <- rownames(x)
  
  if(estimate.umap) {
    # Project cells not in the subsampling to the UMAP embedding
    umap.transform.res <- uwot::umap_transform(
      X = x.cells.not.selected,
      model = umap.res,
      n_threads = as.integer(num.threads.transform)
    )
    # Check the distances
    knn.ref <- as.vector(nabor::knn(data=umap.res[[1]], query=umap.res[[1]], k=2)$nn.dists[,-1])
    knn.proj <- as.vector(nabor::knn(data=umap.res[[1]], query=umap.transform.res, k=1)$nn.dists)
    idx.cells.to.remove <- which(knn.proj >= quantile(knn.ref, outlier.quantile))
    umap.transform.res[idx.cells.to.remove, ] <- NA
    
    umap.transform.df <- data.frame(umap.transform.res)
    colnames(umap.transform.df) <- paste0("UMAP_Dimension_",seq_len(ncol(umap.transform.df)))
    rownames(umap.transform.df) <- rownames(x.cells.not.selected)
    rm(umap.transform.res)
    umap.df <- rbind(umap.df, umap.transform.df)
    umap.df <- umap.df[cell.names, , drop=FALSE]
  }
  
  return(umap.df)
}