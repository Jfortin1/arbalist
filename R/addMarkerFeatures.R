#' Differential feature analysis
#'
#' For each of the groups specified, perform differential analysis for the specified features by comparing that group of cells to cells with similar distribution across the features from the other cell groups.
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
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList listToMap experiments
#' @export
addMarkerFeatures <- function(
    mae,
    experiment.name = 'PeakMatrix',
    group.colname = 'Clusters',
    group.experiment.name = 'TileMatrix500',
    num.threads = 4,
    sampleLabels = "Sample"
) {
  # Pull the necessary data from the MAE
  sce.list <- findSCE(mae, experiment.name)
  se <- sce.list$sce
  main.exp.name <- names(mae)[sce.list$sce.idx]
  alt.exp.name <- sce.list$alt.exp.name
  if(is.null(group.experiment.name)) {
    group.experiment.name <- experiment.name
  }
  group.sce.list <- findSCE(mae, group.experiment.name)
  
  # Perform differential analysis
  diff.res <- markerDiff(
    mat = assay(se),
    group.classification = colData(group.sce.list$sce)[,group.colname],
    num.threads = num.threads
    )
  
  # Add the new SummarizedExperiment to a list of the existing experiments
  exp.list <- experiments(mae)
  new.exp.name <- paste0(experiment.name,'_marker_features')
  exp.list[[new.exp.name]] <- SummarizedExperiment(assays = diff.res, rowRanges = rowRanges(se))
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
  
  # Create a new MultiAssayExperiment matching the old MAE but with the new experiment
  new.mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
  colData.filler <- matrix(NA, ncol = ncol(colData(mae)), nrow = ncol(exp.list[[new.exp.name]]))
  colnames(colData.filler) <- colnames(colData(mae))
  rownames(colData.filler) <- colnames(exp.list[[new.exp.name]])
  colData(new.mae) <- rbind(colData(mae),colData.filler)
  metadata(new.mae) <- metadata(mae)
  
  return(new.mae)
}

#' @importFrom presto wilcoxauc
#' @importFrom nabor knn
#' @importFrom methods is as
markerDiff <- function(mat, group.classification, num.threads = 4) {
  
  # set up a pointer to the data for profiling the columns for matching cells between groups
  beachmat::flushMemoryCache()
  if(is(mat,"DelayedArray")) {
    ptr <- beachmat.hdf5::initializeCpp(mat, memorize=TRUE)
  } else {
    ptr <- beachmat::initializeCpp(mat, memorize = TRUE)
  }
  # estimate a maximum of the data so that we can get the distribution of the data in one pass and minimizing memory usage.
  row.sums <- tatami.row.sums(ptr, num.threads=num.threads)
  estimated.max <- max(tatami.column.sums(tatami.subset(ptr, which(row.sums == max(row.sums)), by.row = TRUE), num.threads = num.threads))+5
  
  # get the data distribution of each column
  data.column.hist <- matrix(NA, nrow = tatami.dim(ptr)[2], ncol= estimated.max)
  for(i in seq(estimated.max)) {
    data.column.hist[,i] <- tatami.column.sums(tatami.compare(ptr, op = '>=', val = i, by.row = TRUE, right = TRUE), num.threads = num.threads)
  }

  # For each group, compare the cells to matched cells from the other groups to assess feature differences
  diff.res.list <- list()
  for( group.name in unique(group.classification) ) {
    
    # find the nearest neighbours from the cells not in the group to used as a matched group for comparison
    match.num <- min(sum(group.classification != group.name), sum(group.classification == group.name))
    group.idxs <- which(group.classification == group.name)
    match.option.idxs <- which(group.classification != group.name)
    match.num <- min(length(match.option.idxs), length(group.idxs))
    knnx <- nabor::knn(data = data.column.hist[match.option.idxs,], query = data.column.hist[group.idxs,], k = match.num)$nn.idx
    bkgd.matched.idxs <- match.option.idxs[unique(as.numeric(knnx[sample(nrow(knnx)),]))[1:match.num]]
    
    mat1 <- as(mat[,group.idxs],"sparseMatrix") # cells in one group
    mat2 <- as(mat[,bkgd.matched.idxs],"sparseMatrix")# matched cells with similar column wide distribution
    
    # perform wilcoxon test for the features
    diff.res <- wilcoxauc(cbind(mat1,mat2), c(rep(group.name, ncol(mat1)),rep(paste0(group.name,'_matched'), ncol(mat2))))
    diff.res <- diff.res [which(diff.res$group == group.name),]
    
    # pull together differential results
    m1 <- Matrix::rowSums(mat1, na.rm=TRUE)
    m2 <- Matrix::rowSums(mat2, na.rm=TRUE)
    offset <- 1
    log2FC <- log2((m1 + offset) / (m2 + offset))
    log2Mean <- log2(((m1 + offset) + (m2 + offset)) / 2)
    
    diff.res.list[[group.name]] <- data.frame(
      log2Mean = log2Mean,
      log2FC = log2FC,
      fdr = diff.res$padj, 
      pval = diff.res$pval, 
      mean1 = Matrix::rowMeans(mat1, na.rm=TRUE), 
      mean2 = Matrix::rowMeans(mat2, na.rm=TRUE), 
      n = ncol(mat1),
      auc = diff.res$auc
    )
  }

  return(diff.res.list)
}