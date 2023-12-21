#' Clustering cells
#' 
#' Creating cell clusters from iterative LSI reduced dimensions contained with a MAE and adding the cluster result back to the MAE
#' 
#' @param mae MultiAssayExperiment
#' @param name.iterative.lsi String specifying the name of the reduced dimensions to create clusters from
#' @param clusters.colname String specifying the column name to save the clusters as in the experiment column data
#' @param method String specifying the method for creating clsuters. Valid options are "Seurat" or "Scran".
#' @param force Logical whether to overwrite existing columns with clusters.colname column name
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
  #sampleCells = NULL,
  #seed = 1, 
  method = "Seurat",
  #dimsToUse = NULL,
  #scaleDims = NULL, 
  #corCutOff = 0.75,
  #knnAssign = 10, 
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
  #prefix = "C",
  force = FALSE
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
  
  if(grep('seurat',tolower(method))) {
    
    # create a fake matrix to fill the input daya slot to represent the cell dimension
    mat.fake <- matrix(rnorm(nrow(lsi.matrix) * 3, 10), ncol = nrow(lsi.matrix), nrow = 3)
    colnames(mat.fake) <- rownames(lsi.matrix)
    rownames(mat.fake) <- paste0("t",seq_len(nrow(mat.fake)))
    
    seurat.obj <- Seurat::CreateSeuratObject(mat.fake, project='scATAC', min.cells=0, min.features=0)
    seurat.obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=lsi.matrix, key='PC_', assay='RNA')
    
    seurat.obj <- Seurat::FindNeighbors(seurat.obj, reduction = 'pca', dims = seq_len(ncol(lsi.matrix)))
    seurat.obj <- Seurat::FindClusters(seurat.obj, reduction = 'pca', dims = seq_len(ncol(lsi.matrix)))
    
    # get clusters form Seurat Object
    clust <- seurat.obj@meta.data[,ncol(seurat.obj@meta.data)]
    clust <- paste0("Cluster",match(clust, unique(clust)))
    names(clust) <- rownames(lsi.matrix)
    clust <- clust[rownames(lsi.matrix)]

  }
  
  # save the new reduced dimensions to the MAE
  if(is.null(lsi.alt.exp.name)) {
    if(clusters.colname %in% colnames(colData(mae[[exp.idx]])) && !force) {
      stop(paste0(clusters.colname,' is already a column name. Set force = TRUE if you want to overwrite it.'))
    }
    colData(mae[[exp.idx]])[,clusters.colname] <- as.character(clust)
  } else {
    if(clusters.colname %in% colnames(colData(altExp(mae[[exp.idx]], lsi.alt.exp.name))) && !force) {
      stop(paste0(clusters.colname,' is already a column name. Set force = TRUE if you want to overwrite it.'))
    }
    colData(altExp(mae[[exp.idx]], lsi.alt.exp.name))[,clusters.colname] <- as.character(clust)
  }
  
  mae
  
}