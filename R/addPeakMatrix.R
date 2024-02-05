#' Create a peak matrix
#'
#' Call peaks using MACSr and then create a peak matrix experiment in the MultiAssayExperiment.
#'
#' @param mae \linkS4class{MultiAssayExperiment}
#' @param genome.size Number specifying the effective genome size or the size of hte genome that is mappable to use in MACSr peak calling. See [MACSr::callpeak] for more details.
#' @param pseudobulk.experiment.name Experiment name for the pseudobulk \linkS4class{SummarizedExperiment} containing colData pointing to coverage files. This experiment can be created from arbalist::addGroupCoverages.
#' @param sc.experiment.name Experiment name for \linkS4class{SingleCellExperiment} with cell barcodes as the rownames
#' @param reproducibility A string that indicates how peak reproducibility should be handled. This string is dynamic and can be a function of n where n is the number of samples being assessed. For example, reproducibility = "2" means at least 2 samples must have a peak call at this locus and reproducibility = "(n+1)/2" means that the majority of samples must have a peak call at this locus.
#' @param shift Number of baseparis to shift each Tn5 insertion. When combined with extsize this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACSr::callpeak.
#' @param extsize Number of basepairs to extend the fragment after shift has been applied. When combined with extsize this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACSr::callpeak.
#' @param method Method to use for significance testing in MACS2. Options are "p" for p-value and "q" for q-value. When combined with cutOff this gives the method and significance threshold for peak calling (see MACS2 documentation)
#' @param cutOff Numeric significance cutOff for the testing method indicated by method (see MACSr::callpeak documentation).
#' @param nomodel Whether or not to build the shifting model during MACSr::callpeak
#' @param nolambda If True, MACS will use fixed background lambda as local lambda for every peak region.
#' @param extendSummits Number of basepairs to extend peak summits (in both directions) to obtain final fixed-width peaks. For example, extendSummits = 250 will create 501-bp fixed-width peaks from the 1-bp summits.
#' @param output.dir Directory where hdf5 files should be output while creating the peak matrix \linkS4class{SingleCellExperiment}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @return \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocParallel bptry bpparam bplapply
#' @importFrom alabaster.matrix AmalgamatedArray
#' @importFrom MACSr callpeak
#' @importFrom data.table fread
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
    BPPARAM = bpparam()
){
  
  # Find the single cell experiment result
  sce <- findSCE(mae,sc.experiment.name)$sce

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
        cutoff_analysis = TRUE)
      
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
  
  groupPeaks <- groupPeaks[which(unlist(lapply(groupPeaks,length)) > 0)]
  
  group.names <- unique(mae[[pseudobulk.experiment.name]]$group)
  group.names <- names(which(sapply(group.names,function(x) {any(grep(x,names(groupPeaks)))})))
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
        barcodes = as.character(sub('.*#','',rownames(colData(sce)))[sce$Sample == sample.name])
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
  
  new.mae <- c(mae, 'PeakMatrix'=peak.sce)
  
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
#' @importFrom GenomicRanges resize
.identifyReproduciblePeaks <- function(
    peak.sets,
    reproducibility = 0.51,
    extendSummits = 250,
    by = 'score'
){
  
  n <- length(peak.sets)
  
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
  
  if(length(overlapping.peaks.idx) == 0) {
    combined.peaks <- Reduce("c", peak.sets)
  } else {
    non.overlapping.peaks <- combined.all.peaks[-overlapping.peaks.idx]
    overlapping.peaks <- combined.all.peaks[overlapping.peaks.idx]
    overlapping.peaks <- overlapping.peaks[order(overlapping.peaks[,by],decreasing = TRUE)]
    
    overlapping.peaks <- unique(overlapping.peaks)
    
    hits.res <- queryHits(findOverlaps(overlapping.peaks))
    selected.peaks <- overlapping.peaks[hits.res[match(unique(hits.res),hits.res)]]
    
    combined.peaks <- c(non.overlapping.peaks,selected.peaks)
  }

  overlapMat <- Reduce('cbind',lapply(split(Reduce("c", peak.sets), Reduce("c", peak.sets)$GroupReplicate), function(x){
    overlapsAny(combined.peaks, x)}))
  
  if(length(peak.sets) > 1){
    combined.peaks$Reproducibility <- rowSums(overlapMat)
    combined.peaks$ReproducibilityPercent <- round(rowSums(overlapMat) / ncol(overlapMat) , 3)
    min.rep <- eval(parse(text=reproducibility))
    idxPass <- which(combined.peaks$Reproducibility >= min.rep)
    ncombined.peaks <- combined.peaks[idxPass]
  }else{
    combined.peaks$Reproducibility <- rep(NA, length(combined.peaks))
  }
  
  combined.peaks$groupScoreQuantile <- round(trunc(rank(combined.peaks$replicateScoreQuantile))/length(combined.peaks),3)
  
  combined.peaks
}