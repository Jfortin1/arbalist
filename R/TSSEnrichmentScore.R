#' Add TSS Enrichment Scores
#'
#' Calculate and add TSS enrichment scores to the MultiAssayExperiment.
#' 
#' @param mae A \linkS4class{MultiAssayExperiment}
#' @param experiment.name String specifying the experiment name to add the TSS Enrichment to the colData
#' @param experiment.name.for.gene.grs String specifying the experiment name to use rowRanges start as for TSS. If this experiment name is not found in mae, then gene.grs must be specified. If both experiment.name.for.gene.grs and gene.grs are specified, then gene.grs takes precedence.
#' @inheritParams tssEnrichmentScores
#' @param force Logical. If there is already a TSS Enrichment column whether to overwrite it (TRUE) or error (FALSE).
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @author Natalie Fox
#' @importFrom SummarizedExperiment rowData
#' @export
addTSSEnrichmentScores <- function(
    mae,
    experiment.name = 'TileMatrix500',
    experiment.name.for.gene.grs = 'GeneExpressionMatrix',
    gene.grs = NULL,
    window = 101,
    norm = 100,
    flank = 2000,
    min.norm = 0.2,
    force = FALSE,
    BPPARAM = BiocParallel::bpparam()
) {
  
  if(!'fragment_file' %in% colnames(colData(mae))) {
    stop('Please include the paths to the fragment files in the MultiAssayExperiment colData')
  }
  
  # Find the experiment result
  sce.list <- findSCE(mae, experiment.name)
  if(is.null(sce.list)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }
  
  if(is.null(gene.grs)) {
    gene.sce.list <- findSCE(mae,experiment.name.for.gene.grs)
    if(is.null(gene.sce.list)) {
      stop(paste0('Cannot find a ',experiment.name.for.gene.grs,' to get gene coordinates from so please specify gene.grs'))
    }
    gene.grs <- GRanges(rowData(gene.sce.list$sce)$interval[rowData(gene.sce.list$sce)$interval != 'NA'])
    if(length(gene.grs) == 0) {
      stop(paste0('There are no interval specified in the rowData of ',experiment.name.for.gene.grs,' so please specify gene.grs'))
    }
  }
  
  if('TSSEnrichment' %in% colnames(colData(mae[[sce.list$sce.idx]])) & !force) {
    stop('TSSEnrichment is already a colname. set force = TRUE to overwrite the column')
  }
  
  fragment.files <- colData(mae)$fragment_file
  names(fragment.files) <- rownames(colData(mae))
  
  tss.scores <- bptry(bplapply(
    names(fragment.files), 
    function(sample.name, fragment.files, gene.grs, window,norm, flank, min.norm, sce) {
      per.sample.tss.scores <- tssEnrichmentScores(
        fragment.file = fragment.files[sample.name],
        barcodes = as.character(sub('.*#','',rownames(colData(sce)))[colData(sce)$Sample == sample.name]),
        gene.grs = gene.grs,
        window = window,
        norm = norm,
        flank = flank,
        min.norm = min.norm
      );
      as.numeric(per.sample.tss.scores)
    },
    fragment.files = fragment.files,
    gene.grs = gene.grs,
    window = window,
    norm = norm,
    flank = flank,
    min.norm = min.norm,
    sce = sce.list$sce,
    BPPARAM = BPPARAM
  ))
  for(i in seq_along(fragment.files)) {
    names(tss.scores[[i]]) <- paste0(names(fragment.files)[i],'#',as.character(sub('.*#','',rownames(colData(sce.list$sce)))[colData(sce.list$sce)$Sample == names(fragment.files)[i]]))
  }
  colData(mae[[sce.list$sce.idx]])$TSSEnrichment <- unlist(tss.scores)[rownames(colData(sce.list$sce))]
  
  mae
}

#' Calculate TSS Enrichment Scores
#'
#' Calculate TSS enrichment scores from one fragment.files.
#' 
#' @param fragment.file String specifying fragment file.
#' @param barcodes Vector or strings specified the cell barcodes to include from the fragment file.
#' @param gene.grs A \linkS4class{GRanges} specifying genes. The start coordinate will be used as the transcription start site.
#' @param window Number specifying the size in bp of the TSS window.
#' @param norm Number specifying the size in bp of the flanking region window.
#' @param flank Number specifying the bp distance for the flanking region from the TSS.
#' @param min.norm Number specifying the minimum allow for flanking region window. Avoids divide by zero error.
#'
#' @author Natalie Fox
#' @importFrom GenomicRanges resize GRanges trim start end seqnames strand
#' @importFrom IRanges IRanges
tssEnrichmentScores <- function(
    fragment.file,
    barcodes,
    gene.grs,
    window = 101,
    norm = 100,
    flank = 2000,
    min.norm = 0.2
) {
  
  tss.grs <- gene.grs
  
  # create TSS and flanking GRanges
  tss.grs <- GenomicRanges::resize(tss.grs, 1, fix = "start")
  tss.grs <- unique(tss.grs)
  tss.window <- GenomicRanges::resize(tss.grs, window, "center")
  tss.window <- GenomicRanges::trim(tss.window)
  tss.flank <- c(
    # positive Flank
    GRanges(seqnames(tss.grs), IRanges(end(tss.grs) + flank - norm + 1, end(tss.grs) + flank)),
    # negative Flank
    GRanges(seqnames(tss.grs), IRanges(start(tss.grs) - flank, start(tss.grs) - flank + norm - 1))
  )
  tss.flank <- GenomicRanges::trim(tss.flank)

  window.file <- tempfile(pattern = paste0(gsub("window",'\\.','_',gsub('/','_',fragment.file))), tmpdir = tempdir())
  if(file.exists(window.file)) {
    window.file <- tempfile(pattern = paste0(gsub("window",'\\.','_',gsub('/','_',fragment.file))), tmpdir = tempdir())
    if(file.exists(window.file)) {
      stop('file already exists')
    }
  }
  window.count.matrix <- saveRegionMatrix(
    fragment.file,
    output.file = window.file,
    output.name = 'tss_and_flank',
    regions = tss.window,
    barcodes = barcodes
  )
  
  flank.file <- tempfile(pattern = paste0(gsub("flank",'\\.','_',gsub('/','_',fragment.file))), tmpdir = tempdir())
  if(file.exists(window.file)) {
    flank.file <- tempfile(pattern = paste0(gsub("flank",'\\.','_',gsub('/','_',fragment.file))), tmpdir = tempdir())
    if(file.exists(flank.file)) {
      stop('file already exists')
    }
  }
  flank.count.matrix <- saveRegionMatrix(
    fragment.file,
    output.file = flank.file ,
    output.name = 'tss_and_flank',
    regions = tss.flank,
    barcodes = barcodes
  )
  
  # normalize per bp
  window.per.cell.sum <- apply(window.count.matrix,2,function(x) {sum(as.integer(x[x != 0]))}) / window
  flank.per.cell.sum  <- apply(flank.count.matrix,2,function(x) {sum(as.integer(x[x != 0]))}) / norm
  
  # compute scores
  tss.scores <- 2 * window.per.cell.sum / (pmax(flank.per.cell.sum, min.norm))
  tss.scores <- round(tss.scores, 3)
  
  # remove temp window and flank files
  file.remove(window.file)
  file.remove(flank.file)
  
  return(tss.scores)
}
