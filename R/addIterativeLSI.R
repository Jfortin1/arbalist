#' Add new iterative LSI embeddings to MultiAssayExperiment
#'
#' @param mae \linkS4class{MultiAssayExperiment} 
#' @param experiment.name String specifying the name of the experiment to create the embedding from and add reduced dimensions to
#' @param embedding.name String sepcifying the name of the new iterativeLSI embedding
#' @param rank Number specifying the rank for irlba_realized
#' @param iterations Number of LSI iterations to perform.
#' @param num.features Number of accessible features to select when selecting the most accessible or most variable features.
#' @param cluster.method String specifying cluster method. Currently "kmeans" is the only supported option.
#' @param cluster.k Number of clusters to use
#' @param correlation.cutoff 	Numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded from analysis.
#' @param scale.to Number specifying the center for TF-IDF normalization.
#' @param num.threads Number of threads to be used for parallel computing.
#' @param seed Number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param total.features Number of features to consider for use in LSI after ranking the features by the total number of insertions. These features are the only ones used throughout the variance identification and LSI.
#' @param filter.quantile Number 0,1 that indicates the quantile above which features should be removed based on insertion counts prior
#'
#' @return \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#' @export
#' @importFrom MultiAssayExperiment assay
#' @importFrom SingleCellExperiment reducedDim<- reducedDims reducedDimNames<- reducedDimNames reducedDims
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
  
  if(embedding.name %in% reducedDimNames(mae[[experiment.name]])) {
    stop(paste0('There is already a reduced dimension called ',embedding.name))
  }
  
  res <- iterativeLSI (
    x = assay(mae[[experiment.name]]),
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
  
  reducedDim(mae[[experiment.name]], embedding.name) <- res$embedding
  
  return(mae)
}