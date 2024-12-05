#' Add gene score matrix
#'
#' Calculate gene score matrix to estimate gene expression and add the matrix to the MultiAssayExperiment.
#' 
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param gene.grs String specifying the experiment name to use rowRanges start as for TSS. If this experiment name is not found in mae, then gene.grs must be specified. If both experiment.name.for.gene.grs and gene.grs are specified, then gene.grs takes precedence.
#' @param experiment.name String containing the experiment name of the tile matrix.
#' @inheritParams calculateGeneScoreMatrix
#'
#' @return \linkS4class{MultiAssayExperiment} with DoubletScore and DoubletEnrichment columns added to the experiment colData.
#'
#' @author Natalie Fox
#' @importFrom SummarizedExperiment rowRanges assay
#' @export
addGeneScoreMatrix <- function(
    mae,
    gene.grs,
    experiment.name = 'TileMatrix500',
    extend.outside.gene = c(1000, 100000),
    gene.model = "exp(-abs(x)/5000) + exp(-1)",
    gene.scale.factor = 5,
    upstream = 5000,
    downstream = 0,
    num.threads = 4
) {

  # retrieve the tile matrix from the MAE
  sce.list <- findSCE(mae,experiment.name)
  if(is.null(sce.list)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }

  # Calculate the gene score matrix
  sce <- calculateGeneScoreMatrix(
    tile.matrix = assay(sce.list$sce),
    gene.grs = gene.grs,
    tile.grs = rowRanges(sce.list$sce),
    extend.outside.gene = extend.outside.gene,
    gene.model = gene.model,
    gene.scale.factor = gene.scale.factor,
    upstream = upstream,
    downstream = downstream
  )
  
  # add the Experiment in the MAE
  new.mae <- c(mae, 'GeneScoreMatrix'=sce)
  
  return(new.mae)
}

#' Create gene score matrix
#'
#' Calculate gene score matrix to estimate gene expression.
#' 
#' @param tile.matrix tile matrix count matrix.
#' @param gene.grs A \linkS4class{GRanges} specifying gene coordinates.
#' @param tile.grs A \linkS4class{GRanges} specifying tile coordinates.
#' @param extend.outside.gene A numeric vector with two values specifying the range to extend outside the gene prior to overlapping genes with tiles.
#' @param gene.model A string giving a "gene model function" used for weighting peaks for gene score calculation. This string
#' should be a function of `x`, where `x` is the stranded distance from the transcription start site of the gene. 
#' @param gene.scale.factor A numeric scaling factor to weight genes based on the inverse of there length i.e. [(Scale Factor)/(Gene Length)]. This
#' is scaled from 1 to the scale factor. Small genes will be the scale factor while extremely large genes will be closer to 1. This scaling helps with
#' the relative gene score value.
#' @param upstream Integer describing the number of bp upstream the gene to extend the gene body. This effectively makes the gene body larger as there are proximal peaks that should be weighted equally to the gene body.
#' @param downstream Integer describing the number of bp downstream the gene to extend the gene body.This effectively makes the gene body larger as there are proximal peaks that should be weighted equally to the gene body.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how doublet scoring should be parallelized per iteration of doublets.
#' @param num.threads Integer scalar specifying the number of threads to use for tatami multiply.
#'
#' @return \linkS4class{SingleCellExperiment}
#'
#' @author Natalie Fox
#' @importFrom SummarizedExperiment SummarizedExperiment colData<- rowRanges<-
#' @importFrom SingleCellExperiment mainExpName<-
#' @importFrom GenomicRanges strand start end findOverlaps resize shift width distance
#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom BiocGenerics which
#' @importFrom BiocParallel bpparam
#' @importFrom beachmat tatami.multiply
#' @export
calculateGeneScoreMatrix <- function(
    tile.matrix,
    gene.grs,
    tile.grs,
    extend.outside.gene = c(1000, 100000),
    gene.model = "exp(-abs(x)/5000) + exp(-1)",
    gene.scale.factor = 5,
    upstream = 5000,
    downstream = 0,
    BPPARAM = BiocParallel::bpparam(),
    num.threads = 4
    ) {
  
  common.seqlevels <- intersect(seqlevels(tile.grs),seqlevels(gene.grs))
  gene.grs <- keepSeqlevels(gene.grs,common.seqlevels, pruning.mode='coarse')
  
  idx.minus <- BiocGenerics::which(strand(gene.grs) == "-")
  idx.other <- BiocGenerics::which(strand(gene.grs) != "-")
  start(gene.grs)[idx.other] <- start(gene.grs)[idx.other] - upstream
  end(gene.grs)[idx.other] <- end(gene.grs)[idx.other] + downstream
  end(gene.grs)[idx.minus] <- end(gene.grs)[idx.minus] + upstream
  start(gene.grs)[idx.minus] <- start(gene.grs)[idx.minus] - downstream
  
  tile.size <- width(tile.grs)[1]
  
  # check which of the gene ranges overlap before extending the genes upstream & downstream
  gene.overlap <- findOverlaps(gene.grs)
  gene.overlap <- gene.overlap[queryHits(gene.overlap) < subjectHits(gene.overlap)]

  # check gene overlap with minimum gene extension
  min.extended.gene.grs <- resize(shift(gene.grs,-min(extend.outside.gene)),width(gene.grs)+min(extend.outside.gene)*2)
  gene.min.extended.overlap <- findOverlaps(gene.grs)
  gene.min.extended.overlap <- gene.min.extended.overlap[queryHits(gene.min.extended.overlap) < subjectHits(gene.min.extended.overlap)]

  # check gene overlap with max gene extension 
  extended.gene.grs <- resize(shift(gene.grs,-max(extend.outside.gene)),width(gene.grs)+max(extend.outside.gene)*2)
  gene.extended.overlap <- findOverlaps(gene.grs)
  gene.extended.overlap <- gene.extended.overlap[queryHits(gene.extended.overlap) < subjectHits(gene.extended.overlap)]

  # find the distance between genes and tiles
  gene.tile.overlap <- findOverlaps(extended.gene.grs,tile.grs)
  gene.tile.dist <- distance(gene.grs[queryHits(gene.tile.overlap)], tile.grs[subjectHits(gene.tile.overlap)])
  gene.tile.dist <- gene.tile.dist * sign(start(tile.grs)[subjectHits(gene.tile.overlap)] - start(gene.grs)[queryHits(gene.tile.overlap)])
  gene.tile.dist[which(strand(gene.grs) == "-")] <- gene.tile.dist[which(strand(gene.grs) == "-")] * -1
  
  # evaluate input model
  x <- gene.tile.dist
  gene.model <- "exp(-abs(x)/5000) + exp(-1)"
  x <- eval(parse(text=gene.model))
  
  # get feature weights related to gene width
  m <- 1 / width(gene.grs)
  gene.grs$weight <- 1 + m * (gene.scale.factor - 1) / (max(m) - min(m))
  x <- x * mcols(gene.grs)$weight[queryHits(gene.tile.overlap)]
  
  # creating weights sparse matrix
  weights.matrix <- Matrix::sparseMatrix(
    i = queryHits(gene.tile.overlap), 
    j = subjectHits(gene.tile.overlap), 
    x = x, 
    dims = c(length(gene.grs), nrow(tile.matrix))
  )
  
  beachmat::flushMemoryCache()
  if(is(tile.matrix,"DelayedArray")) {
    ptr.tiles <- beachmat.hdf5::initializeCpp(tile.matrix, memorize = TRUE)
  } else {
    ptr.tiles <- beachmat::initializeCpp(tile.matrix, memorize = TRUE)
  }
  ptr.weights <- beachmat::initializeCpp(weights.matrix, memorize = TRUE)
  gs.matrix <- beachmat::tatami.multiply(ptr.weights, ptr.tiles, right=TRUE, num.threads=num.threads)

  colnames(gs.matrix) <- colnames(tile.matrix)

  mat.list <- list(counts=gs.matrix)
  se <-  SummarizedExperiment(mat.list)
  rowRanges(se) <- gene.grs
  colnames(se) <- colnames(gs.matrix)
  colData(se)$Sample <- sub("#.*","",colnames(gs.matrix))
  sce <- as(se, 'SingleCellExperiment')
  mainExpName(sce) <- 'GeneScoreMatrix'
  
  return(sce)
}