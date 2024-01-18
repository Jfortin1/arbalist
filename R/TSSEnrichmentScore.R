
addTSSEnrichmentScores <- function(
    mae,
    experiment.name = 'TileMatrix500',
    gene.grs = NULL,
    force = FALSE,
    BPPARAM = bpparam(),
    ...
) {
  
  # Find the experiment result
  sce.list <- findSCE(mae,'TileMatrix500')
  if(is.null(sce.list)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }
  
  if(is.null(gene.grs)) {
    gene.grs = GRanges(rowData(mae[['GeneExpressionMatrix']])$interval[rowData(mae[['GeneExpressionMatrix']])$interval != 'NA'])
  }
  
  if('TSSEnrichment' %in% colnames(colData(mae[[sce.list$sce.idx]])) & !force) {
    stop('TSSEnrichment is already a colname. set force = TRUE to overwrite the column')
  }
  
  fragment.files <- colData(mae)$fragment_file
  
  tss.scores <- bptry(bplapply(
    names(fragment.files), 
    function(sample.name) {
      sample.tss.scores <- tssEnrichmentScores(
        fragment.files[sample.name],
        as.character(sub('.*#','',rownames(colData(mae[[sce.list$sce.idx]])))[mae[[sce.list$sce.idx]]$Sample == sample.name]),
        res.gene.grs
      )
      names(sample.tss.scores) <- rownames(colData(mae[[sce.list$sce.idx]]))[mae[[sce.list$sce.idx]]$Sample == sample.name]
      return(sample.tss.scores)
    },
    BPPARAM = BPPARAM
  ))
  colData(mae[[sce.list$sce.idx]])$TSSEnrichment <- unlist(tss.scores)[rownames(colData(mae[[sce.list$sce.idx]]))]
  
  #colData(mae[[sce.list$sce.idx]])$TSSEnrichment <- rep(NA,ncol(mae[[sce.list$sce.idx]]))
  #for(i in seq_along(fragment.files)) {
  #  sample.tss.scores <- tssEnrichmentScores(
  #    fragment.files[i],
  #    as.character(sub('.*#','',rownames(colData(mae[[sce.list$sce.idx]])))[mae[[sce.list$sce.idx]]$Sample == names(fragment.files[i])]),
  #    res.gene.grs
  #  )
  #  colData(mae[[sce.list$sce.idx]])[rownames(colData(mae[[sce.list$sce.idx]]))[mae[[sce.list$sce.idx]]$Sample == names(fragment.files[i])],sce.list$sce.idx] <- sample.tss.scores
  #}
  
  mae
}

tssEnrichmentScores <- function(
    fragment.file,
    barcodes,
    gene.grs,
    window = 101,
    norm = 100,
    flank = 2000,
    min.norm = 0.2
) {
  
  TSS <- gene.grs
  
  #Create Window and Flank
  TSS <- GenomicRanges::resize(TSS, 1, fix = "start")
  #strand(TSS) <- "*"
  TSS <- unique(TSS)
  tssWindow <- GenomicRanges::resize(TSS, window, "center")
  tssWindow$type <- "window"
  tssFlank <- c(
    #Positive Flank
    GRanges(seqnames(TSS), IRanges(end(TSS) + flank - norm + 1, end(TSS) + flank)),
    #Negative Flank
    GRanges(seqnames(TSS), IRanges(start(TSS) - flank, start(TSS) - flank + norm - 1))
  )
  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)
  
  #Trim In Case Extending beyond Chromosomes
  tssFeatures <- GenomicRanges::trim(tssFeatures)
  
  countWindow <- saveRegionMatrix(
    fragment.file,
    output.file = tempfile(pattern = "window", tmpdir = tempdir()),
    output.name = 'tss_and_flank',
    regions = tssWindow,
    barcodes = barcodes
  )
  
  countFlank <- saveRegionMatrix(
    fragment.file,
    output.file = tempfile(pattern = "flank", tmpdir = tempdir()),
    output.name = 'tss_and_flank',
    regions = tssFlank,
    barcodes = barcodes
  )
  
  #Normalize per BP
  cWn <- apply(countWindow,2,function(x) {sum(as.integer(x[x != 0]))}) / window
  cFn <- apply(countFlank,2,function(x) {sum(as.integer(x[x != 0]))}) / norm
  
  #Compute scores
  tssScores <- 2 * cWn / (pmax(cFn, min.norm))
  tssScores <- round(tssScores, 3)
  
  file.remove(paste0(tempdir(),'/window.h5'))
  file.remove(paste0(tempdir(),'/flank.h5'))
  
  return(tssScores)
}
