#' Add doublet scores
#'
#' Calculate doublet scores and add them to the MultiAssayExperiment.
#' 
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param experiment.name String containing thr experiment name to calculate doublet scores.
#' @param seed Numeric scalar to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param num.trials Integer scalar specifying the number of times to simulate the number of cells in the sample for doublets.
#' @param num.threads Integer scalar specifying the number of threads to be used for parallel computing.
#' @param num.dimensions Integer scalar specifying the number of iterative LSI dimensions to use doublet scoring
#' @param plot.out.dir String spcifying the output path for doublet score quality control plots
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how doublet scoring should be parallelized per sample.
#'
#' @return \linkS4class{MultiAssayExperiment} with doublet scores column added to the experiment colData.
#'
#' @author Natalie Fox
#' @importFrom scran buildKNNGraph
#' @importFrom igraph as_adjacency_matrix
#' @importFrom ggplot2 ggplot geom_point aes ggsave theme_classic scale_colour_gradientn
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom BiocGenerics Reduce
addDoubletScores <- function(
    mae,
    experiment.name = 'TileMatrix500',
    seed = 5,
    num.trials = 5,
    num.threads = 1,
    num.dimensions = 30,
    plot.out.dir = '.',
    BPPARAM = bpparam()
) {
  
  # Find the experiment result and get the matrix
  sce.list <- findSCE(mae,experiment.name)
  if(is.null(sce.list)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }
  se <- sce.list$sce
  main.exp.name <- names(mae)[sce.list$sce.idx]
  alt.exp.name <- sce.list$alt.exp.name
  
  per.sample.cell.counts <- sort(table(colData(se)$Sample), decreasing = TRUE)
  if(any(per.sample.cell.counts <= num.dimensions)) {
    warning(paste0(
      paste(names(per.sample.cell.counts)[which(per.sample.cell.counts <= num.dimensions)],collapse=', '),
      ' samples have less than num.dimensions (',num.dimensions,') cells with ',
      paste(as.numeric(per.sample.cell.counts[which(per.sample.cell.counts <= num.dimensions)]),collapse=', '),
      ' cell counts respectively'))
    per.sample.cell.counts <- per.sample.cell.counts[which(per.sample.cell.counts > num.dimensions)]
  }
  if(length(per.sample.cell.counts) > 2) {
    new.order <- NULL
    end.i <- length(per.sample.cell.counts)
    for(i in seq(1,length(per.sample.cell.counts)/2)) {
      if(i >= end.i) {
        if(new.order[length(new.order)] != i) {
          new.order <- c(new.order,i)
        }
        break
      }
      new.order <- c(new.order,i)
      if(i >= end.i) {
        break
      }
      new.order <- c(new.order, end.i)
      end.i <- end.i - 1
      if(i >= end.i) {
        break
      }
      new.order <- c(new.order, end.i)
      end.i <- end.i - 1
    }
    per.sample.cell.counts <- per.sample.cell.counts[new.order]
  }
  
  parallelized.res <- bptry(bplapply(
    names(per.sample.cell.counts)[which(per.sample.cell.counts > num.dimensions)],
    .doublet.calculations.per.sample,
    se = se,
    plot.out.dir = plot.out.dir,
    num.trials = num.trials,
    num.threads = num.threads,
    num.dimensions = num.dimensions,
    BPPARAM = BPPARAM
  ))
  
  doublet.res <- Reduce(rbind, parallelized.res)
  
  colData(mae[[sce.list$sce.idx]])$DoubletScore <- rep(NA,nrow(colData(mae[[sce.list$sce.idx]])))
  colData(mae[[sce.list$sce.idx]])$DoubletEnrichment <- rep(NA,nrow(colData(mae[[sce.list$sce.idx]])))
  
  colData(mae[[sce.list$sce.idx]])[rownames(doublet.res),'DoubletScore'] <- doublet.res$score
  colData(mae[[sce.list$sce.idx]])[rownames(doublet.res),'DoubletEnrichment'] <- doublet.res$enrichment
  
  return(mae)
}

#' @importFrom stats pbinom p.adjust
.doublet.calculations.per.sample <- function(
    sample.name,
    se,
    plot.out.dir = '.',
    num.trials = 5,
    num.threads = 1, 
    num.dimensions = 30,
    num.doublets.per.iteration = 10
){
  
  num.synthetic.doublets <- num.trials * sum(colData(se)$Sample == sample.name)
  mat <- assay(se)[,which(colData(se)$Sample == sample.name)]
  
  # Calculate the iterativeLSI embedding for the sample
  lsi.res <- iterativeLSI(mat, rank = num.dimensions)
  if(is.null(lsi.res)) {
    return(NULL)
  }
  # Simulate Doublets
  beachmat::flushMemoryCache()
  ptr <- beachmat::initializeCpp(mat, memorize=TRUE)
  mat.doublet <- NULL
  num.iterations <- ceiling(num.synthetic.doublets/num.doublets.per.iteration)
  for(i in seq(num.iterations)) {
    if(i == num.iterations) {
      remainder <- num.synthetic.doublets %% num.doublets.per.iteration
      if(remainder != 0) {
        num.doublets.per.iteration <- num.synthetic.doublets %% num.doublets.per.iteration
      }
    }
    random.cells <- rev(sort(sample(1:ncol(mat), 2*num.doublets.per.iteration, replace = FALSE)))
    ptr.subset <- apply_subset(ptr, random.cells, row = FALSE)
    doublet.sparse.column <- as(aggregate_counts(ptr.subset, rep(seq(1,num.doublets.per.iteration),2)-1L, nthreads = num.threads, binarize = FALSE), "sparseMatrix")/2
    if(is.null(mat.doublet)) {
      mat.doublet <- doublet.sparse.column
    } else {
      mat.doublet <- cbind(mat.doublet, doublet.sparse.column)
    }
  }
  colnames(mat.doublet) <- paste0('doublet',seq(ncol(mat.doublet)))
  
  # Project the Doublets into the UMAP embedding
  require(Matrix)
  mat.doublet.normalized <- .apply.tf.idf.normalization(mat.doublet[lsi.res$details$row.subset,], lsi.res$details$ncol, lsi.res$details$row.sums, scale.to = lsi.res$details$scale.to, lsi.method = lsi.res$details$lsi.method) 
  idxNA <- Matrix::which(is.na(mat.doublet.normalized), arr.ind = TRUE)
  if(length(idxNA) > 0){
    mat.doublet.normalized[idxNA] <- 0
  }
  V <- Matrix::t(mat.doublet.normalized) %*% lsi.res$details$svd$u %*% Matrix::diag(1/lsi.res$details$svd$d)
  svdDiag <- matrix(0, nrow = lsi.res$details$num.dimensions, ncol = lsi.res$details$num.dimensions)
  diag(svdDiag) <- lsi.res$details$svd$d
  mat.lsi.doublet <- Matrix::t(svdDiag %*% Matrix::t(V))
  mat.lsi.doublet <- as.matrix(mat.lsi.doublet)
  rownames(mat.lsi.doublet) <- colnames(mat.doublet.normalized)
  colnames(mat.lsi.doublet) <- paste0("LSI",seq_len(ncol(mat.lsi.doublet)))
  
  # Calculate the UMAP embedding
  combined.mat <- rbind(lsi.res$embedding,mat.lsi.doublet)
  umap.res <- scater::calculateUMAP(t(combined.mat))
  umap.res <- as.data.frame(umap.res)
  umap.res$doublet <- c(rep(FALSE,nrow(lsi.res$embedding)), rep(TRUE,num.synthetic.doublets))
  umap.res <- umap.res[rev(seq(nrow(umap.res))),]
  
  if(!is.null(plot.out.dir)) {
    ggplot(umap.res, aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color=doublet)) + theme_classic() + ggtitle(sample.name)
    ggsave(paste0(plot.out.dir,'/',sample.name,'_doublet_umap_plot.pdf'), height = 4, width = 5)
  }
  
  knnDoub <- nabor::knn(lsi.res$embedding, mat.lsi.doublet, k = 10)$nn.idx
  
  countKnn <- rep(0, nrow(lsi.res$details$matSVD))
  names(countKnn) <- rownames(lsi.res$details$matSVD)
  
  tabDoub <- table(as.vector(knnDoub))
  countKnn[as.integer(names(tabDoub))] <-  countKnn[as.integer(names(tabDoub))] + tabDoub
  
  scale.to <- 10000
  scale.by <- scale.to / num.synthetic.doublets
  
  # P-Values
  pvalBinomDoub <- unlist(lapply(seq_along(countKnn), function(x){
    #Round Prediction
    countKnnx <- round(countKnn[x] * scale.by)
    sumKnnx <- round(sum(countKnn) * scale.by)
    pbinom(countKnnx - 1, sumKnnx, 1 / scale.to, lower.tail = FALSE)
  }))
  
  # Adjust
  padjBinomDoub <- p.adjust(pvalBinomDoub, method = "bonferroni")
  
  # Convert To Scores
  doublet.res <- data.frame(
    score = -log10(pmax(padjBinomDoub, 4.940656e-324)),
    enrichment = (countKnn / sum(countKnn)) / (1 / nrow(lsi.res$details$matSVD))
  )
  doublet.res$enrichment <- 10000 * doublet.res$enrichment / length(countKnn) #Enrichment Per 10000 Cells in Data Set
  rownames(doublet.res) <- names(countKnn)
  
  if(!is.null(plot.out.dir)) {
    umap.res <- umap.res[rownames(doublet.res),]
    umap.res$doublet.scores <- doublet.res$score
    umap.res$doublet.enrichment <- doublet.res$enrichment
    ggplot(umap.res, aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color=doublet.res$score)) + theme_classic() + ggtitle(sample.name) + ggplot2::scale_colour_gradientn(colours=c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"))
    ggsave(paste0(plot.out.dir,'/',sample.name,'_doublet_score_umap_plot.pdf'), height = 4, width = 5)
    ggplot(umap.res, aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color=doublet.res$enrichment)) + theme_classic() + ggtitle(sample.name) + ggplot2::scale_colour_gradientn(colours=c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"))
    ggsave(paste0(plot.out.dir,'/',sample.name,'_doublet_enrichment_umap_plot.pdf'), height = 4, width = 5)
  }
  
  return(doublet.res)
}