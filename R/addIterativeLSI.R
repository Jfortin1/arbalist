#' Add new iterative LSI embeddings to MultiAssayExperiment
#' 
#' Calculate iterative LSI dimensionality reduction and add the embedding to the reduceDim field of the specified experiment.
#'
#' @param mae \linkS4class{MultiAssayExperiment} 
#' @param experiment.name String containing the name of the experiment to create the embedding from and add reduced dimensions to.
#' @param embedding.name String containing the name of the new iterativeLSI embedding.
#' @inheritParams iterativeLSI
#'
#' @return \linkS4class{MultiAssayExperiment} with iterative LSI results added to the \linkS4class{SingleCellExperiment} reduced dimensions.
#' 
#' @author Natalie Fox
#' @export
#' @importFrom MultiAssayExperiment assay
#' @importFrom SingleCellExperiment reducedDim<- reducedDims reducedDimNames<- reducedDimNames reducedDims altExp altExpNames
addIterativeLSI <- function(
  mae,
  experiment.name = 'TileMatrix500',
  embedding.name = 'iterativeLSI',
  col.subset = NULL,
  rank = 30,
  iterations = 2,
  first.selection = "Top", # or "Var"
  num.features = 25000,
  lsi.method = 1,
  cluster.method = "Seurat",
  correlation.cutoff = 0.75,
  scale.to = 10000,
  num.threads = 4,
  seed = 5,
  total.features = 500000,
  filter.quantile = 0.995,
  outlier.quantiles = c(0.02, 0.98),
  binarize = TRUE
) {

  # Find the experiment result
  sce.list <- findSCE(mae,experiment.name)
  if(is.null(sce.list)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }
  se <- sce.list$sce
  main.exp.name <- names(mae)[sce.list$sce.idx]
  alt.exp.name <- sce.list$alt.exp.name

  if(embedding.name %in% reducedDimNames(se)) {
    stop(paste0('There is already a reduced dimension called ',embedding.name))
  }
  
  x <- assay(se)
  if(is(x,"DelayedArray")) {
    res <- iterativeLSI(
      x = x,
      col.subset = col.subset,
      rank = rank,
      iterations = iterations,
      first.selection = first.selection,
      num.features = num.features,
      lsi.method = lsi.method,
      cluster.method = cluster.method,
      correlation.cutoff = correlation.cutoff,
      scale.to = scale.to,
      num.threads = num.threads,
      seed = seed,
      total.features = total.features,
      filter.quantile =  filter.quantile,
      outlier.quantiles = outlier.quantiles,
      binarize = binarize
    )
  } else {
    res <- iterativeLSI.sparse.in.mem(
      x = x,
      rank = rank,
      iterations = iterations,
      first.selection = first.selection,
      num.features = num.features,
      lsi.method = lsi.method,
      cluster.method = cluster.method,
      correlation.cutoff = correlation.cutoff,
      scale.to = scale.to,
      num.threads = num.threads,
      seed = seed,
      total.features = total.features,
      filter.quantile =  filter.quantile,
      outlier.quantiles = outlier.quantiles,
      binarize = binarize
    )
  }
  
  if(is.null(alt.exp.name)) {
    if(any(!colnames(mae[[experiment.name]]) %in% rownames(res$embedding))) {
      mae[[experiment.name]] <- mae[[experiment.name]][,rownames(res$embedding)]
    }
    reducedDim(mae[[experiment.name]], embedding.name) <- res$embedding
  } else {
    if(any(!colnames(altExp(mae[[main.exp.name]],alt.exp.name)) %in% rownames(res$embedding))) {
      mae[[main.exp.name]] <- mae[[main.exp.name]][,rownames(res$embedding)]
    }
    reducedDim(altExp(mae[[main.exp.name]],alt.exp.name), embedding.name) <- res$embedding
  }
  
  return(mae)
}
