#' Add doublet scores
#'
#' Calculate doublet scores and add them to the MultiAssayExperiment.
#' 
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param experiment.name String containing the experiment name to calculate doublet scores.
#' @param seed Numeric scalar to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param num.trials Integer scalar specifying the number of times to simulate the number of cells in the sample for doublets.
#' @param num.threads Integer scalar specifying the number of threads to be used for parallel computing.
#' @param num.dimensions Integer scalar specifying the number of iterative LSI dimensions to use doublet scoring.
#' @param num.doublets.per.iteration Integer scalar specifying how many doublets to estimate at one time this allows adjustment of memory usage/run time.
#' @param max.num.synthetic.doublets Integer scalar specifying the maximum number of doublets to run per sample. This is used if num.trials * number of cells in the sample is more than this number.
#' @param plot.out.dir String specifying the output path for doublet score quality control plots.
#' @param selected.samples Vector of strings specifying which samples to calculate doublet scores for. By default (NULL), all that have enough cells will be run.
#' @param cell.depth.column String specifying the column name in the experiment colData that contains the values to use for cell depth. For example for arbalist created MAEs this is probably "fragments". For ArchR created MAEs this might be "nFrags".
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how doublet scoring should be parallelized per iteration of doublets.
#'
#' @return \linkS4class{MultiAssayExperiment} with DoubletScore and DoubletEnrichment columns added to the experiment colData.
#'
#' @author Natalie Fox
#' @importFrom scran buildKNNGraph
#' @importFrom igraph as_adjacency_matrix
#' @importFrom ggplot2 ggplot geom_point aes ggsave theme_classic scale_colour_gradientn
#' @importFrom BiocParallel bplapply bpparam SerialParam bptry
#' @importFrom BiocGenerics Reduce
#' @importFrom stats p.adjust pbinom
#' @importFrom scater calculateUMAP
#' @export
addDoubletScores <- function(
    mae,
    experiment.name = 'TileMatrix500',
    seed = 5,
    num.trials = 5,
    num.threads = 4,
    num.dimensions = 30,
    num.doublets.per.iteration = 200,
    max.num.synthetic.doublets = 2000,
    plot.out.dir = '.',
    selected.samples = NULL,
    cell.depth.column = 'fragments',
    BPPARAM = bpparam()
) {
  
  # Find the experiment result in the MAE so we can retrieve the matrix, this accounts for using the main or alternative experiment
  sce.list <- findSCE(mae, experiment.name)
  if(is.null(sce.list)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }

  # Add DoubletScore and DoubletEnrichment columns to the colData, if not already there
  if(! 'DoubletScore' %in% colnames(colData(mae[[sce.list$sce.idx]]))) {
    colData(mae[[sce.list$sce.idx]])$DoubletScore <- as.numeric(rep(NA,nrow(colData(mae[[sce.list$sce.idx]]))))
  }
  if(! 'DoubletEnrichment' %in% colnames(colData(mae[[sce.list$sce.idx]]))) {
    colData(mae[[sce.list$sce.idx]])$DoubletEnrichment <- as.numeric(rep(NA,nrow(colData(mae[[sce.list$sce.idx]]))))
  }
  
  # Get sample names and check them have enough cells to calculate doublets
  if(is.null(selected.samples)) {
    per.sample.cell.counts <- sort(table(colData(sce.list$sce)$Sample), decreasing = TRUE) # sort largest to smallest number of cells so that if there isn't enough memory we'll error quicker
    if(any(per.sample.cell.counts <= num.dimensions)) {
      warning(paste0(
        paste(names(per.sample.cell.counts)[which(per.sample.cell.counts <= num.dimensions)],collapse=', '),
        ' samples have less than num.dimensions (',num.dimensions,') cells with ',
        paste(as.numeric(per.sample.cell.counts[which(per.sample.cell.counts <= num.dimensions)]),collapse=', '),
        ' cell counts respectively'))
      per.sample.cell.counts <- per.sample.cell.counts[which(per.sample.cell.counts > num.dimensions)]
    }
    selected.samples <- names(per.sample.cell.counts)
  }
  
  # Calculate doublet scores per sample
  for(sample.name in selected.samples) {
    
    # Decide how many doublets to simulate and therefore how many iterations we will need to run
    num.synthetic.doublets <- num.trials * sum(colData(sce.list$sce)$Sample == sample.name)
    if(num.synthetic.doublets > max.num.synthetic.doublets) {
      num.synthetic.doublets <- max.num.synthetic.doublets
    }
    sample.columns <- which(colData(sce.list$sce)$Sample == sample.name)
    if(length(sample.columns) < 2*num.doublets.per.iteration) {
      num.doublets.per.iteration  <- floor(length(sample.columns)/2)
    }
    num.iterations <- ceiling(num.synthetic.doublets/num.doublets.per.iteration)
  
    # Calculate the iterativeLSI embedding for the sample's cells
    cell.depths <- colData(sce.list$sce)[,cell.depth.column]
    cell.names <- colnames(sce.list$sce)
    sample.names <- colData(sce.list$sce)[,'Sample']
    lsi.res <- iterativeLSI(
      assay(sce.list$sce)[,sample.columns],
      cell.names = cell.names[sample.columns],
      sample.names = sample.names[sample.columns],
      cell.depths = cell.depths[sample.columns],
      rank = num.dimensions
      )
    
    gc()
    if(is.null(lsi.res)) {
      next
    }
    # Simulate doublets and project them into the LSI reduced dimensions
    set.seed(seed)
    parallelized.res <- bptry(bplapply(
      seq(num.iterations),
      .doublet.simulation,
      mat = assay(sce.list$sce)[lsi.res$details$row.subset,sample.columns],
      lsi.res = lsi.res,
      num.iterations = num.iterations,
      num.synthetic.doublets = num.synthetic.doublets,
      num.doublets.per.iteration = num.doublets.per.iteration,
      num.threads = num.threads,
      BPPARAM = BPPARAM
    ))
    mat.lsi.doublet.combined <- Reduce(rbind, parallelized.res)
    if(!is(mat.lsi.doublet.combined, "matrix")) {
      stop('something went wrong with the parallelized doublet simulation')
    }
    rownames(mat.lsi.doublet.combined) <- paste0('doublet',seq(nrow(mat.lsi.doublet.combined)))

    # Calculate the UMAP embedding
    umap.res <- scater::calculateUMAP(t(rbind(lsi.res$embedding, mat.lsi.doublet.combined)))
    umap.res <- as.data.frame(umap.res)
    umap.res$doublet <- c(rep(FALSE,nrow(lsi.res$embedding)), rep(TRUE,nrow(mat.lsi.doublet.combined)))
    umap.res <- umap.res[rev(seq(nrow(umap.res))),]
    
    if(!is.null(plot.out.dir)) {
      # Plot doublets relative to the cells
      ggplot(umap.res, aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color=doublet), size = 0.5) + theme_classic() + ggtitle(sample.name)
      ggsave(paste0(plot.out.dir,'/',sample.name,'_doublet_umap_plot.pdf'), height = 4, width = 5)
    }
    
    # Find nearest neighbours to the doublets for calculating a score
    knnDoub <- nabor::knn(lsi.res$embedding, mat.lsi.doublet.combined, k = 10)$nn.idx
    countKnn <- rep(0, nrow(lsi.res$details$matSVD))
    names(countKnn) <- rownames(lsi.res$embedding)
    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <-  countKnn[as.integer(names(tabDoub))] + tabDoub
    
    scale.to <- 10000
    scale.by <- scale.to / num.synthetic.doublets
    
    # P-Values
    pvalBinomDoub <- unlist(lapply(seq_along(countKnn), function(x){
      # Round prediction
      countKnnx <- round(countKnn[x] * scale.by)
      sumKnnx <- round(sum(countKnn) * scale.by)
      pbinom(countKnnx - 1, sumKnnx, 1 / scale.to, lower.tail = FALSE)
    }))
    padjBinomDoub <- p.adjust(pvalBinomDoub, method = "bonferroni")
    
    # Convert adjusted p-values too doublet scores
    doublet.res <- data.frame(
      score = -log10(pmax(padjBinomDoub, 4.940656e-324)),
      enrichment = (countKnn / sum(countKnn)) / (1 / nrow(lsi.res$details$matSVD))
    )
    doublet.res$enrichment <- 10000 * doublet.res$enrichment / length(countKnn) #Enrichment Per 10000 Cells in Data Set
    rownames(doublet.res) <- names(countKnn)
    
    if(!is.null(plot.out.dir)) {
      # Plot UMAP of cells coloured by doublet scores and doublet enrichment
      umap.res <- umap.res[rownames(doublet.res),]
      umap.res$doublet.scores <- doublet.res$score
      umap.res$doublet.enrichment <- doublet.res$enrichment
      
      umap.res <- umap.res[order(umap.res$doublet.scores),]
      ggplot(umap.res, aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color=doublet.scores), size = 0.5) + theme_classic() + ggtitle(sample.name) + ggplot2::scale_colour_gradientn(colours=c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"))
      ggsave(paste0(plot.out.dir,'/',sample.name,'_doublet_score_umap_plot.pdf'), height = 5.5, width = 7)
      
      umap.res <- umap.res[order(umap.res$doublet.enrichment),]
      ggplot(umap.res, aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color=doublet.enrichment), size = 0.5) + theme_classic() + ggtitle(sample.name) + ggplot2::scale_colour_gradientn(colours=c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"))
      ggsave(paste0(plot.out.dir,'/',sample.name,'_doublet_enrichment_umap_plot.pdf'), height = 5.5, width = 7)
    }
    
    # Add doublet scores/enrichment to the colData
    colData(mae[[sce.list$sce.idx]])[rownames(doublet.res),'DoubletScore'] <- doublet.res$score
    colData(mae[[sce.list$sce.idx]])[rownames(doublet.res),'DoubletEnrichment'] <- doublet.res$enrichment
  }

  return(mae)
}

#' @importFrom beachmat tatami.subset
.doublet.simulation <- function(
  i,
  mat,
  lsi.res,
  num.iterations,
  num.synthetic.doublets,
  num.doublets.per.iteration,
  num.threads
){
  
  if(i == num.iterations) {
    # Adjust for the number of doublets not being evenly divisible by the number of iterations
    remainder <- num.synthetic.doublets %% num.doublets.per.iteration
    if(remainder != 0) {
      num.doublets.per.iteration <- remainder
    }
  }
  
  # Randomly select cells for doublets and then combine them in pairs to simulate doublets
  beachmat::flushMemoryCache()
  if(is(mat,"DelayedArray")) {
    ptr <- beachmat.hdf5::initializeCpp(mat, memorize=TRUE)
  } else {
    ptr <- beachmat::initializeCpp(mat, memorize=TRUE)
  }
  random.cells <- sample(seq(1,ncol(mat)), 2, replace = FALSE)
  ptr.subset <- tatami.subset(ptr, subset = random.cells, by.row = FALSE)
  row.sums <- tatami.row.sums(ptr.subset, num.threads = num.threads)
  mat.doublet <- as(matrix(row.sums,ncol=1),'dgCMatrix')/2
  for(j in 2:num.doublets.per.iteration) {
    random.cells <- sample(seq(1,ncol(mat)), 2, replace = FALSE)
    ptr.subset <- tatami.subset(ptr, subset = random.cells, by.row = FALSE)
    row.sums <- tatami.row.sums(ptr.subset, num.threads = num.threads)
    mat.doublet <- cbind(mat.doublet,as(matrix(row.sums,ncol=1),'dgCMatrix')/2)
  }
  gc()
  
  # Project the doublets into the UMAP embedding
  require(Matrix)
  mat.doublet.normalized <- .apply.tf.idf.normalization.in.mem(mat.doublet, lsi.res$details$ncol, lsi.res$details$row.sums, scale.to = lsi.res$details$scale.to, lsi.method = lsi.res$details$lsi.method) 
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

  return(mat.lsi.doublet)
}