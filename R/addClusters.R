#' Cluster cells
#' 
#' Creating cell clusters from iterative LSI reduced dimensions contained within a MAE and adding the cluster result back to the MAE.
#' 
#' @param mae \linkS4class{MultiAssayExperiment}
#' @param name.reduced.dim String containing the name of the reduced dimensions to create clusters from. If there are multiple reduced dimensions in the MultiAssayExperiment with the same name, then specify the experiment name in the vector's names. This can also be a vector or strings if you want to cluster based on the combination of reduced dimensions.
#' @param clusters.colname String containing the column name to save the clusters as in the experiment column data.
#' @param cluster.prefix String to prefix to the cluster number for saving in colData result.
#' @param method String containing the method for creating clusters. Valid options are "Seurat" or "scran".
#' @param dims.to.use Numeric vector or list of numeric vectors specifying which of the columns to use from the reduced dimensions. Columns are the reduced dimension features.
#' @param seed Integer scalar for random number generation required for sampling cells to cluster
#' @param num.cells.to.sample Integer scalar specifying the number of cells to sample for subsampling of cells when clustering.
#' @param knn.k Integer scalar specifying the number of nearest neighbors to use during assigning clusters to the cells not subsampled.
#' @param force Logical whether to overwrite existing columns with clusters.colname column name.
#' 
#' @return A \linkS4class{MultiAssayExperiment} with cluster column added to the experiment colData.
#' 
#' @author Natalie Fox
#' @export
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject FindNeighbors FindClusters
#' @importFrom SingleCellExperiment altExp reducedDim<- reducedDim altExpNames altExp<- reducedDimNames
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats rnorm
#' @importFrom nabor knn
addClusters <- function(
  mae,
  name.reduced.dim = "iterativeLSI",
  clusters.colname = "Clusters",
  cluster.prefix = "C",
  method = "Seurat",
  dims.to.use = NULL,
  seed = 17,
  num.cells.to.sample = NULL,
  knn.k = 10,
  force = FALSE
) {
  
  # retrieve the reduced dimension (probably iterative LSI) matrix from the MAE
  reduced.dim.list <- list()
  for(i in name.reduced.dim) {
    reduced.dim.list[[i]] <- findReducedDimRes(mae,i)
  }

  # if dims.to.use is set, filter to the specified dimensions for the first reduced dimension specified
  if(is.null(dims.to.use)) {
    reduced.dim.matrix <- reduced.dim.list[[1]]$matrix
  } else if(is(dims.to.use,'list')) {
    if(length(dims.to.use) != length(reduced.dim.list)) {
      stop('if including multiple dim.to.use vectors then there must be the same number of entries in the list as values in name.reduced.dim')
    }
    reduced.dim.matrix <- reduced.dim.list[[1]]$matrix[,dims.to.use[[1]]]
  } else {
    reduced.dim.matrix <- reduced.dim.list[[1]]$matrix[,dims.to.use]
  }
  
  # if more than one reduced dimension specified then combine them into one matrix (account for dims.to.use)
  if(length(reduced.dim.list) > 1) {
    for(i in 2:length(reduced.dim.list)) {
      if(nrow(reduced.dim.matrix) != nrow(reduced.dim.list[[i]]$matrix)) {
        stop('The reduced dimension matrices need to have the same number of rows.')
      }
      if(is.null(dims.to.use)) {
        reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix)
      } else if(is(dims.to.use,'list')) {
        reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix[,dims.to.use[[i]]])
      } else {
        reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix[,dims.to.use])
      }
    }
  }
  
  # Sample cells for clustering
  set.seed(seed)
  estimating.clusters <- FALSE
  if(!is.null(num.cells.to.sample)) {
    if(num.cells.to.sample < nrow(reduced.dim.matrix)) {
      reduced.dim.matrix.all <- reduced.dim.matrix
      reduced.dim.matrix <- reduced.dim.matrix[sample(seq_len(nrow(reduced.dim.matrix)),num.cells.to.sample),]
      estimating.clusters <- TRUE
    } else {
      stop(paste0('num.cells.sample (',num.cells.to.sample,') needs to be less than the number of cells in the reduced dimentison after dim.to.use filtering (',nrow(reduced.dim.matrix),').'))
    }
  }
  
  if(any(grep('seurat',tolower(method)))) {
    
    # create a fake matrix to fill the input data slot to represent the cell dimension
    mat.fake <- matrix(rnorm(nrow(reduced.dim.matrix) * 3, 10), ncol = nrow(reduced.dim.matrix), nrow = 3)
    colnames(mat.fake) <- rownames(reduced.dim.matrix)
    rownames(mat.fake) <- paste0("t",seq_len(nrow(mat.fake)))
    
    seurat.obj <- Seurat::CreateSeuratObject(mat.fake, project='scATAC', min.cells=0, min.features=0)
    seurat.obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=reduced.dim.matrix, key='PC_', assay='RNA')
    
    seurat.obj <- Seurat::FindNeighbors(seurat.obj, reduction = 'pca', dims = seq_len(ncol(reduced.dim.matrix)))
    seurat.obj <- Seurat::FindClusters(seurat.obj, reduction = 'pca', dims = seq_len(ncol(reduced.dim.matrix)))
    
    # get clusters form Seurat Object
    clust <- seurat.obj@meta.data[,ncol(seurat.obj@meta.data)]
    clust <- paste0(cluster.prefix,match(clust, unique(clust)))
    names(clust) <- rownames(reduced.dim.matrix)

  } else if(any(grep('scran',tolower(method)))) {
    
    g <- scran::buildSNNGraph(x = t(reduced.dim.matrix), d = ncol(reduced.dim.matrix))
    clust <- paste0(cluster.prefix,igraph::cluster_walktrap(g)$membership)
    names(clust) <- rownames(reduced.dim.matrix)
    
  } else {
    stop(paste0(method,' is not one of the current clustering method options. Try "Seurat" or "scran".'))
  }
  
  # if subsampled the cells then estimate the clusters for the remaining cells
  if(estimating.clusters) {
    missing.clusters <- setdiff(rownames(reduced.dim.matrix.all),rownames(reduced.dim.matrix))
    # calculate the k nearest neighbors
    knn.idx <- nabor::knn(data = reduced.dim.matrix, query = reduced.dim.matrix.all[missing.clusters,], k = knn.k)$nn.idx
    # find the most common cluster classificaiton in the nearest neighbors
    clust.estimated <- apply(knn.idx, 1, function(x) { names(rev(sort(table(clust[x]))))[1] })
    names(clust.estimated) <- missing.clusters
  }
  
  # save the new reduced dimensions to the MAE
  for(i in 1:length(reduced.dim.list)) {
    if(is.null(reduced.dim.list[[i]]$alt.exp.name)) {
      if(clusters.colname %in% colnames(colData(mae[[reduced.dim.list[[i]]$exp.idx]])) && !force) {
        stop(paste0(clusters.colname,' is already a column name in ',names(mae)[reduced.dim.list[[i]]$exp.idx],' experiment colData. Set force = TRUE if you want to overwrite it.'))
      }
      colData(mae[[reduced.dim.list[[i]]$exp.idx]])[,clusters.colname] <- rep(NA,nrow(colData(mae[[reduced.dim.list[[i]]$exp.idx]])))
      colData(mae[[reduced.dim.list[[i]]$exp.idx]])[names(clust),clusters.colname] <- as.character(clust)
      if(estimating.clusters) {
        colData(mae[[reduced.dim.list[[i]]$exp.idx]])[names(clust.estimated),clusters.colname] <- as.character(clust.estimated)
      }
    } else {
      if(clusters.colname %in% colnames(colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))) && !force) {
        stop(paste0(clusters.colname,' is already a column name in ',names(mae)[reduced.dim.list[[i]]$exp.idx],' experiment colData. Set force = TRUE if you want to overwrite it.'))
      }
      colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))[,clusters.colname] <- rep(NA,nrow(colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))))
      colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))[names(clust),clusters.colname] <- as.character(clust)
      if(estimating.clusters) {
        colData(altExp(mae[[reduced.dim.list[[i]]$exp.idx]], reduced.dim.list[[i]]$alt.exp.name))[names(clust.estimated),clusters.colname] <- as.character(clust.estimated)
      }
    }
  }
  
  mae

}

cluster.matrix <- function(
    mat,
    method = 'seurat',
    cluster.prefix = "C",
    resolution = 0.8
) {

  if(any(grep('seurat',tolower(method)))) {
    
    # create a fake matrix to fill the input data slot to represent the cell dimension
    mat.fake <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
    colnames(mat.fake) <- rownames(mat)
    rownames(mat.fake) <- paste0("t",seq_len(nrow(mat.fake)))
    
    seurat.obj <- Seurat::CreateSeuratObject(mat.fake, project='scATAC', min.cells=0, min.features=0)
    seurat.obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
    
    seurat.obj <- Seurat::FindNeighbors(seurat.obj, reduction = 'pca', dims = seq_len(ncol(mat)))
    seurat.obj <- Seurat::FindClusters(seurat.obj, reduction = 'pca', dims = seq_len(ncol(mat)), resolution = resolution)
    
    # get clusters form Seurat Object
    clust <- seurat.obj@meta.data[,ncol(seurat.obj@meta.data)]
    clust <- paste0(cluster.prefix,match(clust, unique(clust)))
    names(clust) <- rownames(mat)
    
  } else if(any(grep('scran',tolower(method)))) {
    
    g <- scran::buildSNNGraph(x = t(mat), d = ncol(mat))
    clust <- paste0(cluster.prefix,igraph::cluster_walktrap(g)$membership)
    names(clust) <- rownames(mat)
    
  } else {
    stop(paste0(method,' is not one of the current clustering method options. Try "Seurat" or "scran".'))
  }
  
  return(clust)
}