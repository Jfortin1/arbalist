#' Differential features
#'
#' Calculate gene score matrix to estimate gene expression and add the matrix to the MultiAssayExperiment.
#' 
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param experiment.name String containing the experiment name for selecting the features to use for differential analysis.
#' @param group.colname String specifying the colname for cell groups.
#' @param group.experiment.name String containing the experiment name for finding the groups in the colData. If NULL, will use experiment.name.
#' @param num.threads Integer scalar specifying the number of threads to be used for parallel computing.
#' @param sampleLabels String containing experiment colData colname for sample labels.
#' 
#' @return \linkS4class{MultiAssayExperiment} with DoubletScore and DoubletEnrichment columns added to the experiment colData.
#'
#' @author Natalie Fox
#' @importFrom presto wilcoxauc
#' @importFrom nabor knn
#' @export
addMarkerFeatures <- function(
    mae,
    experiment.name = 'PeakMatrix',
    group.colname = 'Clusters',
    group.experiment.name = 'TileMatrix500',
    num.threads = 1,
    sampleLabels = "Sample"
) {
  sce.list <- findSCE(mae, experiment.name)
  se <- sce.list$sce
  main.exp.name <- names(mae)[sce.list$sce.idx]
  alt.exp.name <- sce.list$alt.exp.name
  
  if(is.null(group.experiment.name)) {
    group.experiment.name <- experiment.name
  }
  group.sce.list <- findSCE(mae, group.experiment.name)
  
  diff.res <- markerDiff(
    mat = assay(se),
    group.classification = colData(group.sce.list$sce)[,group.colname]
    )
  
  exp.list <- experiments(mae)
  new.exp.name <- paste0(experiment.name,'_marker_features')
  new.se <- SummarizedExperiment(assays=diff.res,rowRanges=rowRanges(se))
  exp.list[[new.exp.name]] <- SummarizedExperiment(assays=diff.res,rowRanges=rowRanges(se))
  colData(exp.list[[new.exp.name]])$ID <- colnames(diff.res[[1]])
  
  el <- ExperimentList(exp.list)
  maplist <- lapply(exp.list, function(se) {
    if(sampleLabels %in% colnames(colData(se))) {
      data.frame(primary = colData(se)[,sampleLabels], colname = colnames(se), stringsAsFactors = FALSE)
    } else {
      data.frame(primary = se$ID, colname = colnames(se), stringsAsFactors = FALSE)
    }
  })
  sampMap <- listToMap(maplist)
  
  # Create and annotate the MultiAssayExperiment
  new.mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
  colData.filler <- matrix(NA,ncol=ncol(colData(mae)),nrow=ncol(new.se))
  colnames(colData.filler) <- colnames(colData(mae))
  rownames(colData.filler) <- colnames(new.se)
  
  colData(new.mae) <- rbind(colData(mae),colData.filler)
  metadata(new.mae) <- metadata(mae)
  
  return(new.mae)
}

#' @importFrom presto wilcoxauc
#' @importFrom nabor knn
markerDiff <- function(mat, group.classification) {
  
  beachmat::flushMemoryCache()
  ptr <- beachmat::initializeCpp(mat, memorize = TRUE)
  stats <- lsi_matrix_stats(ptr, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row

  data.column.hist <- t(as.data.frame(data_distribtion(ptr, max = stats$max)$distribution))
  rownames(data.column.hist) <- NULL

  diff.res <- list()
  for( group.name in unique(group.classification)) {
    match.num <- min(sum(group.classification != group.name), sum(group.classification == group.name))
    group.idxs <- which(group.classification == group.name)
    match.option.idxs <- which(group.classification != group.name)
    match.num <- min(length(match.option.idxs), length(group.idxs))
    knnx <- nabor::knn(data = data.column.hist[match.option.idxs,], query = data.column.hist[group.idxs,], k = match.num)$nn.idx
    bkgd.matched.idxs <- match.option.idxs[unique(as.numeric(knnx[sample(nrow(knnx)),]))[1:match.num]]
    
    mat1 <- as(mat[,group.idxs],"sparseMatrix")# cells in group
    mat2 <- as(mat[,bkgd.matched.idxs],"sparseMatrix")# matched cells with similar column wide distribution
    
    n1 <- ncol(mat1)
    n2 <- ncol(mat2)
    
    df <- wilcoxauc(cbind(mat1,mat2), c(rep(group.name, ncol(mat1)),rep(paste0(group.name,'_matched'), ncol(mat2))))
    df <- df[which(df$group==group.name),]
    
    #Sparse Row Sums
    m1 <- Matrix::rowSums(mat1, na.rm=TRUE)
    m2 <- Matrix::rowSums(mat2, na.rm=TRUE)
    offset <- 1
    log2FC <- log2((m1 + offset) / (m2 + offset))
    log2Mean <- log2(((m1 + offset) + (m2 + offset)) / 2)
    
    diff.res[[group.name]] <- data.frame(
      log2Mean = log2Mean,
      log2FC = log2FC,
      fdr = df$padj, 
      pval = df$pval, 
      mean1 = Matrix::rowMeans(mat1, na.rm=TRUE), 
      mean2 = Matrix::rowMeans(mat2, na.rm=TRUE), 
      n = ncol(mat1),
      auc = df$auc
    )
  }

  return(diff.res)
}