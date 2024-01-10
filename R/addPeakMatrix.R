#' Create a peak matrix
#'
#' Call peaks using MACSr and then create a peak matrix experiment in the MultiAssayExperiment.
#'
#' @param mae \linkS4class{MultiAssayExperiment}
#' @param genome.size Number specifying the effective genome size or the size of hte genome that is mappable to use in MACSr peak calling. See [MACSr::callpeak] for more details.
#'
#' @return \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocParallel bptry bpparam bplapply
#' @importFrom alabaster.matrix AmalgamatedArray
#' @export
addPeakMatrix <- function(
    mae,
    genome.size,
    pseudobulk.experiment.name = 'TileMatrix500_pseudobulk',
    sc.experiment.name = 'TileMatrix500',
    reproducibility = "2",
    shift = -75,
    extsize = 150,
    method = "q",
    cutOff = 0.1,
    nomodel = TRUE,
    nolambda = TRUE,
    extendSummits = 250,
    output.dir = tempdir(),
    BPPARAM = bpparam(),
    ...
){
  
  ## Calculate peak set ##
  
  # set Significance Threshold Parameters
  qvalue <- NULL
  pvalue <- NULL
  if(tolower(method) == "p"){
    pvalue <- cutOff
  }else{
    qvalue <- cutOff
  }
  peak.files <- mae[[pseudobulk.experiment.name]]$coverage.file
  names(peak.files) <- mae[[pseudobulk.experiment.name]]$ID
  
  groupPeaks <- bptry(bplapply(
    peak.files, 
    function(coverage.file) {
      # run MACSr peak calling
      res.peakcalling <- MACSr::callpeak(
        tfile = coverage.file,
        gsize = genome.size,
        store_bdg = FALSE,
        format = "BED",
        call_summits = TRUE,
        keepduplicates = "all",
        name = basename(coverage.file),
        outdir = dirname(coverage.file),
        nolambda = nolambda,
        nomodel = nomodel,
        shift = shift,
        extsize = extsize,
        qvalue = qvalue,
        pvalue = pvalue,
        cutoff_analysis = TRUE,
        ...)
      
      # output files
      summitsFile <- paste0(coverage.file, "_summits.bed")
      narrowPeaksFile <- paste0(coverage.file, "_peaks.narrowPeak")
      xlsFile <- paste0(coverage.file, "_peaks.xls")
      
      # read Summits!
      peaks <- data.table::fread(summitsFile, select = c(1,2,3,5))
      peaks <- GRanges(peaks$V1, IRanges(peaks$V2 + 1, peaks$V3), score = peaks$V5)
      peaks <- sort(sortSeqlevels(peaks))
      
      # remove Files
      r2 <- suppressWarnings(file.remove(summitsFile, narrowPeaksFile, xlsFile))
      
      return(peaks)
    },
    BPPARAM = BPPARAM
  ))
  
  group.names <- unique(mae[[pseudobulk.experiment.name]]$group)
  names(group.names) <- group.names
  cluster.peak.sets <- bptry(bplapply(group.names, function(x) {
    return(.identifyReproduciblePeaks(groupPeaks[mae[[pseudobulk.experiment.name]]$ID[mae[[pseudobulk.experiment.name]]$group == x]], by = "score", reproducibility = reproducibility, extendSummits=extendSummits))
  }, BPPARAM = BPPARAM))
  
  res.peak.set <- .identifyReproduciblePeaks(cluster.peak.sets, by = "groupScoreQuantile")
  
  ## create peak matrix ##
  
  # check if the hdf5 files already exist
  fragment.files <- mae$fragment_file[grep('tsv.gz',mae$fragment_file)]
  output.file.names <- paste0(output.dir,'/PeakMatrix_',names(fragment.files),'.h5')
  names(output.file.names) <- names(fragment.files)
  if(any(file.exists(output.file.names))) {
    stop(paste0(output.file.names[which(file.exists(output.file.names))[1]],' already exists. We do not want to overwrite the file in case it is being used. Either remove the file if you think it is safe to do so or specify a different output.dir.'))
  }
  
  res.peak.set <- sortSeqlevels(res.peak.set)
  
  # For each fragment.file create a matrix hdf5 file
  # Parallelizing per sample
  peak.res.list <- bptry(bplapply(
    names(fragment.files),
    function(sample.name,fragment.files,output.file.names,regions) {
      matrix.res <- saveRegionMatrix(
        fragment.file = as.character(fragment.files[sample.name]),
        output.file = as.character(output.file.names[sample.name]),
        output.name = 'peak_matrix',
        regions = res.peak.set,
        barcodes = as.character(sub('.*#','',rownames(colData(mae[[sc.experiment.name]])))[mae[[sc.experiment.name]]$Sample == sample.name])
      )
      return(matrix.res)
    }, 
    fragment.files=fragment.files, 
    output.file.names=output.file.names, 
    regions=res.peak.set,
    BPPARAM = BPPARAM
  ))
  names(peak.res.list) <- names(fragment.files)
  
  peak.sce <- .getSCEFromH5List(peak.res.list, res.peak.set)
  
  exp.list <- experiments(mae)
  exp.list[['PeakMatrix']] <- peak.sce
  
  el <- ExperimentList(exp.list)
  maplist <- lapply(exp.list, function(se) {
    if('Sample' %in% colnames(colData(se))) {
      data.frame(primary = se$Sample, colname = colnames(se), stringsAsFactors = FALSE)
    } else {
      data.frame(primary = se$ID, colname = colnames(se), stringsAsFactors = FALSE)
    }
  })
  sampMap <- listToMap(maplist)
  
  # Create and annotate the MultiAssayExperiment
  new.mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
  metadata(new.mae) <- metadata(mae)
  
  return(new.mae)
}

#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges<-
#' @importFrom methods as
.getSCEFromH5List <- function(h5.res.list, grs) {
  
  # Combined the per sample results into one matrix
  if(length(h5.res.list) == 1) {
    mat <- h5.res.list[[1]]
  } else {
    mat <- AmalgamatedArray(h5.res.list, along=2)
  }
    
  # Map the cells back to samples and update colnames to include sample names
  cell.to.sample <- unlist(sapply(names(h5.res.list),function(x) {rep(x,ncol(h5.res.list[[x]]))}))
  
  # Create a SingleCellExperiment 
  mat.list <- list(counts=mat)
  se <-  SummarizedExperiment(mat.list)
  rowRanges(se) <- grs
  colnames(se) <- paste0(cell.to.sample,'#',colnames(se))
  colData(se)$Sample <- as.character(cell.to.sample)
  sce <- as(se, 'SingleCellExperiment')
  
  return(sce)
}

#.fastAnnoPeaks <- function(
#    peaks = NULL, 
#    BSgenome = NULL, 
#    geneAnnotation = NULL, 
#    promoterRegion = c(2000, 100)
#){
#  
#  #Validate
#  peaks <- .validGRanges(peaks)
#  peakSummits <- GenomicRanges::resize(peaks,1,"center")
#  geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
#  geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
#  geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
#  BSgenome <- validBSgenome(BSgenome)
#  
#  #First Lets Get Distance to Nearest Gene Start
#  distPeaks <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
#  mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
#  mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
#  promoters <- extendGR(GenomicRanges::resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
#  op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
#  og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
#  oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
#  type <- rep("Distal", length(peaks))
#  type[which(og & oe)] <- "Exonic"
#  type[which(og & !oe)] <- "Intronic"
#  type[which(op)] <- "Promoter"
#  mcols(peaks)$peakType <- type
#  
#  #First Lets Get Distance to Nearest TSS's
#  .logMessage("Annotating Peaks : TSS", logFile = logFile)
#  distTSS <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
#  mcols(peaks)$distToTSS <- mcols(distTSS)$distance
#  if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
#    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distTSS)]
#  }else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
#    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distTSS)]
#  }
#  
#  #Get NucleoTide Content
#  nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
#  mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
#  mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
#  peaks
#  
#}

#' @importFrom S4Vectors queryHits subjectHits 
#' @importFrom IRanges overlapsAny
.identifyReproduciblePeaks <- function(
    peak.sets,
    reproducibility = 0.51,
    extendSummits = 250,
    by = 'score'
){
  
  for(i in names(peak.sets)) {
    peak.sets[[i]] <- GenomicRanges::resize(peak.sets[[i]], extendSummits * 2 + 1, "center")
    peak.sets[[i]] <- sortSeqlevels(peak.sets[[i]])
    peak.sets[[i]]$replicateScoreQuantile <- round(trunc(BiocGenerics::rank(peak.sets[[i]][,by]))/length(peak.sets[[i]]),3)
    peak.sets[[i]]$GroupReplicate <- sub('.*/','',i)
  }
  
  combined.all.peaks <- Reduce("c", peak.sets)
  
  combined.all.peaks <- unique(combined.all.peaks)
  peak.overlaps <- findOverlaps(combined.all.peaks)
  peak.overlaps <- peak.overlaps[which(queryHits(peak.overlaps) < subjectHits(peak.overlaps))]
  overlapping.peaks.idx <- unique(c(queryHits(peak.overlaps),subjectHits(peak.overlaps)))
  non.overlapping.peaks <- combined.all.peaks[-overlapping.peaks.idx]
  overlapping.peaks <- combined.all.peaks[overlapping.peaks.idx]
  overlapping.peaks <- overlapping.peaks[order(overlapping.peaks[,by],decreasing = TRUE)]
  
  overlapping.peaks <- unique(overlapping.peaks)
  
  hits.res <- queryHits(findOverlaps(overlapping.peaks))
  selected.peaks <- overlapping.peaks[hits.res[match(unique(hits.res),hits.res)]]
  
  combined.peaks <- c(non.overlapping.peaks,selected.peaks)

  overlapMat <- Reduce('cbind',lapply(split(Reduce("c", peak.sets), Reduce("c", peak.sets)$GroupReplicate), function(x){
    overlapsAny(combined.peaks, x)}))
  
  if(length(peak.sets) > 1){
    combined.peaks$Reproducibility <- rowSums(overlapMat)
    combined.peaks$ReproducibilityPercent <- round(rowSums(overlapMat) / ncol(overlapMat) , 3)
    idxPass <- which(combined.peaks$Reproducibility >= reproducibility)
    ncombined.peaks <- combined.peaks[idxPass]
  }else{
    combined.peaks$Reproducibility <- rep(NA, length(combined.peaks))
  }
  
  combined.peaks$groupScoreQuantile <- round(trunc(rank(combined.peaks$replicateScoreQuantile))/length(combined.peaks),3)
  
  combined.peaks
}