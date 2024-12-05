#' Create a peak matrix
#'
#' Call peaks using MACSr and then create a peak matrix experiment in the MultiAssayExperiment.
#'
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param genome.size Integer scalar specifying the effective genome size or the size of the genome that is mappable to use in MACSr peak calling. See [MACSr::callpeak] for more details.
#' @param pseudobulk.experiment.name String containing Experiment name for the pseudobulk \linkS4class{SummarizedExperiment} containing colData pointing to coverage files. This experiment can be created from arbalist::addGroupCoverages.
#' @param sc.experiment.name String containing experiment name for \linkS4class{SingleCellExperiment} with cell barcodes as the rownames.
#' @param gene.grs A \linkS4class{GRanges} specifying gene coordinates to use for peak to gene annotation. If the The start coordinate will be used as the transcription start site.
#' @param exon.grs A \linkS4class{GRanges} specifying exon coordinates to use for peak to gene annotation.
#' @param reproducibility String that indicates how peak reproducibility should be handled. This string is dynamic and can be a function of n where n is the number of samples being assessed. For example, reproducibility = "2" means at least 2 samples must have a peak call at this locus and reproducibility = "(n+1)/2" means that the majority of samples must have a peak call at this locus.
#' @param shift Integer scalar specifying how many base pairs to shift each Tn5 insertion. When combined with extsize this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACSr::callpeak.
#' @param extsize Integer scalar specifying how many base pairs to extend the fragment after shift has been applied. When combined with extsize this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACSr::callpeak.
#' @param method String containing the method to use for significance testing in MACS2. Options are "p" for p-value and "q" for q-value. When combined with cut.off this gives the method and significance threshold for peak calling (see MACS2 documentation).
#' @param cut.off Numeric scalar specifying significance cut off for the testing method indicated by method (see MACSr::callpeak documentation).
#' @param nomodel Logical specifying whether or not to build the shifting model during MACSr::callpeak.
#' @param nolambda If True, MACS will use fixed background lambda as local lambda for every peak region.
#' @param extendSummits Integer scalar specifying how many base pairs to extend peak summits (in both directions) to obtain final fixed-width peaks. For example, extendSummits = 250 will create 501-bp fixed-width peaks from the 1-bp summits.
#' @param output.dir String containing the directory where hdf5 files should be output while creating the peak matrix \linkS4class{SingleCellExperiment}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @return \linkS4class{MultiAssayExperiment} with a new peak matrix \linkS4class{SingleCellExperiment}
#' 
#' @author Natalie Fox
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom BiocParallel bptry bpparam bplapply
#' @importFrom MACSr callpeak
#' @importFrom data.table fread
#' @export
addPeakMatrix <- function(
    mae,
    genome.size,
    pseudobulk.experiment.name = 'TileMatrix500_pseudobulk',
    sc.experiment.name = 'TileMatrix500',
    gene.grs = NULL,
    exon.grs = NULL,
    reproducibility = "2",
    shift = -75,
    extsize = 150,
    method = "q",
    cut.off = 0.1,
    nomodel = TRUE,
    nolambda = TRUE,
    extendSummits = 250,
    output.dir = tempdir(),
    BPPARAM = bpparam()
){
  
  if(!is.numeric(genome.size)) {
    stop('genome.size must be numeric')
  }
  
  # Find the single cell experiment result
  sce <- findSCE(mae, experiment.name = sc.experiment.name)$sce

  ## Calculate peak set ##
  
  # set Significance Threshold Parameters
  qvalue <- NULL
  pvalue <- NULL
  if(tolower(method) == "p"){
    pvalue <- cut.off
  }else{
    qvalue <- cut.off
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
  
  groupPeaks <- groupPeaks[which(vapply(groupPeaks,length,FUN.VALUE=1L) > 0)]
  
  group.names <- unique(mae[[pseudobulk.experiment.name]]$group)
  group.names <- names(which(vapply(group.names,function(x) {any(grep(x,names(groupPeaks)))},FUN.VALUE=TRUE)))
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
  mainExpName(peak.sce) <- 'PeakMatrix'
  
  new.mae <- c(mae, 'PeakMatrix'=peak.sce)
  
  # Add gene info relative to the peaks
  new.mae <- .peakGeneAnno(new.mae, gene.grs = gene.grs, exon.grs = exon.grs)
  
  return(new.mae)
}

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

#' @importFrom S4Vectors mcols subjectHits mcols<-
#' @importFrom GenomicRanges GRanges distanceToNearest resize
#' @importFrom IRanges overlapsAny
#' @importFrom BiocGenerics start<- end<-
.peakGeneAnno <- function(    
    mae,
    peak.experiment.name = 'PeakMatrix',
    experiment.name.for.gene.grs = 'GeneExpressionMatrix',
    gene.grs = NULL,
    exon.grs = NULL,
    promoter.region = c(2000, 100)
    ) {
  
  peaks.grs <- rowRanges(mae[[peak.experiment.name]])
  
  if(is.null(gene.grs)) {
    gene.sce.list <- findSCE(mae, experiment.name.for.gene.grs)
    if(is.null(gene.sce.list)) {
      warning(paste0('Cannot find a ',experiment.name.for.gene.grs,' to get gene coordinates from so please specify gene.grs if you want peak to gene annotation'))
      return(mae)
    }
    gene.grs$name <- rowData(gene.sce.list$sce)$name
    gene.grs <- GRanges(rowData(gene.sce.list$sce)$interval[rowData(gene.sce.list$sce)$interval != 'NA'])
    gene.grs$name <- rowData(gene.sce.list$sce)$name[rowData(gene.sce.list$sce)$interval != 'NA']
    if(is.null(exon.grs)) {
      exon.grs <- rowRanges(gene.sce.list$sce)
    }
    if(length(gene.grs) == 0) {
      warning(paste0('There are no interval specified in the rowData of ',experiment.name.for.gene.grs,' so please specify gene.grs if you want peak to gene annotation'))
      return(mae)
    }
  }
  if(!'name' %in% colnames(mcols(gene.grs))) {
    if('symbol' %in% colnames(mcols(gene.grs))) {
      gene.grs$name <- gene.grs$symbol
    } else if(grep('name',colnames(mcols(gene.grs)))) {
      gene.grs$name <- mcols(gene.grs)[grep('name',colnames(mcols(gene.grs)))[1]]
    } else if(grep('id',colnames(mcols(gene.grs)))) {
      gene.grs$name <- mcols(gene.grs)[grep('id',colnames(mcols(gene.grs)))[1]]
    } else {
      warning('there is no names for the genes so skipping gene annotation')
      return(mae)
    }
  }
  
  dist.to.peak.start <- GenomicRanges::distanceToNearest(peaks.grs, GenomicRanges::resize(gene.grs, 1, "start"), ignore.strand = TRUE)
  mcols(peaks.grs)$distToGeneStart <- mcols(dist.to.peak.start)$distance
  mcols(peaks.grs)$nearestGene <- gene.grs$name[subjectHits(dist.to.peak.start)]
  promoters <- GenomicRanges::resize(gene.grs, 1, "start")
  idx.grs.reverse <- BiocGenerics::which(strand(promoters) == "-")
  idx.grs.forward <- BiocGenerics::which(strand(promoters) != "-")
  start(promoters)[idx.grs.forward] <- start(promoters)[idx.grs.forward] - promoter.region[1]
  end(promoters)[idx.grs.forward] <- end(promoters)[idx.grs.forward] + promoter.region[2]
  end(promoters)[idx.grs.reverse] <- end(promoters)[idx.grs.reverse] + promoter.region[1]
  start(promoters)[idx.grs.reverse] <- start(promoters)[idx.grs.reverse] - promoter.region[2]
  logical.prompoter.overlap <- overlapsAny(peaks.grs, promoters, ignore.strand = TRUE)
  logical.gene.overlap <- overlapsAny(peaks.grs, gene.grs, ignore.strand = TRUE)
  logical.exon.overlap <- overlapsAny(peaks.grs, exon.grs, ignore.strand = TRUE)
  peak.type <- rep("Distal", length(peaks.grs))
  peak.type[which(logical.gene.overlap & logical.exon.overlap)] <- "Exonic"
  peak.type[which(logical.gene.overlap & !logical.exon.overlap)] <- "Intronic"
  peak.type[which(logical.prompoter.overlap)] <- "Promoter"
  mcols(peaks.grs)$peakType <- peak.type

  rowRanges(mae[[peak.experiment.name]]) <- peaks.grs
  
  return(mae)
  
}
