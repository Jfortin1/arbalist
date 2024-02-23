#' Add new iterative LSI embeddings to MultiAssayExperiment
#'
#' @param mae \linkS4class{MultiAssayExperiment} 
#' @param experiment.name String containing the name of the experiment to create the embedding from and add reduced dimensions to.
#' @param embedding.name String containing the name of the new iterativeLSI embedding.
#' @param rank Integer scalar specifying the rank for irlba_realized.
#' @param iterations Integer scalar specifying number of LSI iterations to perform.
#' @param num.features Integer scalar specifying the number of accessible features to select when selecting the most accessible or most variable features.
#' @param cluster.method String containing cluster method. Currently "kmeans" is the only supported option.
#' @param cluster.k Integer scalar specifying how many clusters to use.
#' @param correlation.cutoff 	Numeric scalar specifying the cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded from analysis.
#' @param scale.to Numeric scalar specifying the center for TF-IDF normalization.
#' @param num.threads Integer scalar specifying the number of threads to be used for parallel computing.
#' @param seed Numeric scalar to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param total.features Integer scalar specifying how many features to consider for use in LSI after ranking the features by the total number of insertions. These features are the only ones used throughout the variance identification and LSI.
#' @param filter.quantile Numeric scalar between 0 and 1 inclusive that indicates the quantile above which features should be removed based on insertion counts prior.
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
  rank = 30,
  iterations = 2,
  num.features = 25000,
  cluster.method = "kmeans",
  cluster.k = 20,
  correlation.cutoff = 0.75,
  scale.to = 10000,
  num.threads = 1,
  seed = 5,
  total.features = 500000,
  filter.quantile = 0.995 
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
  
  res <- iterativeLSI(
    x = assay(se),
    rank = rank,
    iterations = iterations,
    num.features = num.features,
    cluster.method = cluster.method,
    cluster.k = cluster.k,
    correlation.cutoff = correlation.cutoff,
    scale.to = scale.to,
    num.threads = num.threads,
    seed = seed,
    total.features = total.features,
    filter.quantile = filter.quantile 
  )
  
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
