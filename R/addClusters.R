#' Cluster cells
#' 
#' Creating cell clusters from iterative LSI reduced dimensions contained with a MAE and adding the cluster result back to the MAE
#' 
#' @param mae MultiAssayExperiment
#' @param name.iterative.lsi String specifying the name of the reduced dimensions to create clusters from. If there are multiple reduced dimensions in the MultiAssayExperiment with the same name, then specify the correct experiment name in the vector's names. This can also be a vector or strings if you want to cluster based on the combination of reduced dimensions.
#' @param clusters.colname String specifying the column name to save the clusters as in the experiment column data
#' @param cluster.prefix String prefixed to the cluster number for saving in colData
#' @param method String specifying the method for creating clusters. Valid options are "Seurat" or "scran".
#' @param dims.to.use Numeric vector or list of numeric vector specifying which of the columns to use from the reduced dimensions
#' @param force Logical whether to overwrite existing columns with clusters.colname column name
#' @param ... Additional arguments to be provided to Seurat::FindClusters or scran::buildSNNGraph (for example, knn = 50, jaccard = TRUE)
#' 
#' @return A \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#' @export
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject FindNeighbors FindClusters
#' @importFrom SingleCellExperiment altExp reducedDim<- reducedDim altExpNames altExp<- reducedDimNames
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats rnorm
addClusters <- function(
  mae,
  name.iterative.lsi = "iterativeLSI",
  clusters.colname = "Clusters",
  cluster.prefix = "C",
  method = "Seurat",
  dims.to.use = NULL,
  #scaleDims = NULL, 
  #corCutOff = 0.75,
  #nOutlier = 5, 
  #maxClusters = 25,
  #testBias = TRUE,
  #filterBias = FALSE,
  #biasClusters = 0.01,
  #biasCol = "nFrags",
  #biasVals = NULL,
  #biasQuantiles = c(0.05, 0.95),
  #biasEnrich = 10,
  #biasProportion = 0.5,
  #biasPval = 0.05,
  #nPerm = 500,
  force = FALSE,
  ... 
) {
  
  # retrieve the reduced dimension (probably iterative lsi) matrix from the MAE
  reduced.dim.list <- list()
  for(i in name.iterative.lsi) {
    reduced.dim.list[[i]] <- findReducedDimRes(mae,i)
  }

  if(is.null(dims.to.use)) {
    reduced.dim.matrix <- reduced.dim.list[[1]]$matrix;
  } else if(is(dims.to.use,'list')) {
    if(length(dims.to.use) != length(reduced.dim.list)) {
      stop('if including multiple dim.to.use vectors then there must be the same number of entries in the list as values in name.iterative.lsi')
    }
    reduced.dim.matrix <- reduced.dim.list[[1]]$matrix[,dims.to.use[[1]]];
  } else {
    reduced.dim.matrix <- reduced.dim.list[[1]]$matrix[,dims.to.use];
  }
  
  if(length(reduced.dim.list) > 1) {
    for(i in 2:length(reduced.dim.list)) {
      if(nrow(reduced.dim.matrix) != nrow(reduced.dim.list[[i]]$matrix)) {
        stop('The reduced dimension matrices need to have the same number of rows.')
      }
      if(is.null(dims.to.use)) {
        reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix)
      } else if(is(dims.to.use,'list')) {
        reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix[,dims.to.use[[i]]]);
      } else {
        reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix[,dims.to.use]);
      }
    }
  }
  
  if(grep('seurat',tolower(method))) {
    
    # create a fake matrix to fill the input data slot to represent the cell dimension
    mat.fake <- matrix(rnorm(nrow(reduced.dim.matrix) * 3, 10), ncol = nrow(reduced.dim.matrix), nrow = 3)
    colnames(mat.fake) <- rownames(reduced.dim.matrix)
    rownames(mat.fake) <- paste0("t",seq_len(nrow(mat.fake)))
    
    seurat.obj <- Seurat::CreateSeuratObject(mat.fake, project='scATAC', min.cells=0, min.features=0)
    seurat.obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=reduced.dim.matrix, key='PC_', assay='RNA')
    
    seurat.obj <- Seurat::FindNeighbors(seurat.obj, reduction = 'pca', dims = seq_len(ncol(reduced.dim.matrix)))
    seurat.obj <- Seurat::FindClusters(seurat.obj, reduction = 'pca', dims = seq_len(ncol(reduced.dim.matrix)), ...)
    
    # get clusters form Seurat Object
    clust <- seurat.obj@meta.data[,ncol(seurat.obj@meta.data)]
    clust <- paste0(cluster.prefix,match(clust, unique(clust)))
    names(clust) <- rownames(reduced.dim.matrix)
    clust <- clust[rownames(reduced.dim.matrix)]

  } else if(grep('seurat',tolower(method))) {
    g <- scran::buildSNNGraph(x = t(reduced.dim.matrix), d = ncol(reduced.dim.matrix), k = ifelse(exists("...$k"), ...$k, 25))
    clust <- paste0(cluster.prefix,igraph::cluster_walktrap(g)$membership)
  } else {
    stop(paste0(method,' is not one of the current clustering method options. Try "Seurat" or "scran".'))
  }
  
  # save the new reduced dimensions to the MAE
  for(i in 1:length(reduced.dim.list)) {
    if(is.null(reduced.dim.list[[i]]$alt.exp.name)) {
      if(clusters.colname %in% colnames(colData(mae[[reduced.dim.list[[i]]$exp.idx]])) && !force) {
        stop(paste0(clusters.colname,' is already a column name in ',names(mae)[reduced.dim.list[[i]]$exp.idx],' experiment colData. Set force = TRUE if you want to overwrite it.'))
      }
      colData(mae[[reduced.dim.list[[i]]$exp.idx]])[,clusters.colname] <- as.character(clust)
    } else {
      if(clusters.colname %in% colnames(colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))) && !force) {
        stop(paste0(clusters.colname,' is already a column name in ',names(mae)[reduced.dim.list[[i]]$exp.idx],' experiment colData. Set force = TRUE if you want to overwrite it.'))
      }
      colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))[,clusters.colname] <- as.character(clust)
    }
  }
  
  mae
  
}