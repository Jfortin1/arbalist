#' @importFrom SummarizedExperiment colData
addGroupCoverages <- function(
    mae,
    experiment.name = "TileMatrix500",
    group.by.colname = "Clusters",
    useLabels = TRUE,
    sampleLabels = "Sample",
    min.cells = 40,
    max.cells = 500,
    max.frags = 25*10^6,
    min.reps = 2,
    max.reps = 5,
    sample.ratio = 0.8,
    kmerLength = 6,
    return.groups = FALSE,
    force = FALSE
){
  
  # Find the experiment result
  sce <- NULL;
  sce.idx <- NULL;
  alt.exp.name <- NULL;
  if(experiment.name %in% names(mae)) {
    sce.idx <- which(names(mae) == experiment.name)
    sce <- mae[[sce.idx]]
  } else {
    for(i in seq_len(length(mae))) {
      if(experiment.name %in% altExpNames(mae[[i]])) {
        sce <- altExp(mae[[i]],experiment.name)
        # record where the iterative LSI was saved so we can put the UMAP in the same place
        sce.idx <- i
        alt.exp.name <- experiment.name
        break
      }
    }
  }
  if(is.null(sce.idx)) {
    stop(paste0(experiment.name,' is not found in mae'))
  }
  
  # Samples Passing Min Filter
  samplesPassFilter <- sum(num.cells.per.sample >= min.cells)
  samplesThatCouldBeMergedToPass <- floor(sum(num.cells.per.sample[num.cells.per.sample < min.cells]) / min.cells)
  
  num.samples <- length(unique(colData(sce)$Sample))
  num.cells <- ncol(sce)
  num.cells.per.sample <- table(colData(sce)$Sample)
  num.cells.per.group <- table(colData(sce)[,group.by.colname])
  num.cells.per.group.and.sample <- table(paste0(colData(sce)[,group.by.colname],colData(sce)$Sample))
  cells.per.group <- split(colnames(sce),colData(sce)[,group.by.colname])
  cells.per.group.and.sample <- sapply(cells.per.group,function(x) {split(x,sub('#.*','',x))},simplify = FALSE)
  
  num.cells.to.sample.from.group.and.sample <- ceiling(min.cells/num.samples)
  
  rep.num <- max( min.reps, min( max.reps, floor(min(num.cells.per.group)/min.cells)))
  
  selected.cells <- list()
  for(group.name in names(cells.per.group.and.sample)) {
    num.cells.to.sample.from.group.and.sample <- ceiling(min.cells/num.samples)
    
    x <- cells.per.group.and.sample[[group.name]]
    
    selected.cells[[group.name]] <- list()
    
    # There aren't enough cells from individual samples within the group so merge samples for sampling
    if(any(sapply(x,length) < num.cells.to.sample.from.group.and.sample)) {
      x.merged <- x[which(sapply(x,length) >= num.cells.to.sample.from.group.and.sample)]
      x.merged[[length(x.merged)+1]] <- unlist(x[which(sapply(x,length) < num.cells.to.sample.from.group.and.sample)])
      x <- x.merged
      if(any(sapply(x,length) < num.cells.to.sample.from.group.and.sample)) {
        x <- list(unlist(x))
      }
      num.cells.to.sample.from.group.and.sample <- ceiling(min.cells/length(x))
    }
    
    for(x2 in x) {
      # create replicates while selecting different cells for each replicate (as much as possible)
      selected.cells.per.group.and.sample <- sample(x2, num.cells.to.sample.from.group.and.sample, replace = length(x2) < num.cells.to.sample.from.group.and.sample)
      if(length(selected.cells[[group.name]]) == 0) {
        selected.cells[[group.name]][[1]] <- selected.cells.per.group.and.sample
      } else {
        selected.cells[[group.name]][[1]] <- c(selected.cells[[group.name]][[1]],selected.cells.per.group.and.sample)
      }
      cells.left <- setdiff(x2,selected.cells[[group.name]][[1]])
      for(j in 2:rep.num) {
        if(length(cells.left) < num.cells.to.sample.from.group.and.sample) {
          cells.left <- x2
        }
        selected.cells.per.group.and.sample <- sample(cells.left, num.cells.to.sample.from.group.and.sample, replace = length(x2) < num.cells.to.sample.from.group.and.sample)
        if(j > length(selected.cells[[group.name]])) {
          selected.cells[[group.name]][[j]] <- selected.cells.per.group.and.sample
        } else {
          selected.cells[[group.name]][[j]] <- c(selected.cells[[group.name]][[j]],selected.cells.per.group.and.sample)
        }
        cells.left <- setdiff(cells.left,selected.cells[[group.name]][[j]])
      }
    }
    for(j in 1:rep.num) {
      selected.cells[[group.name]][[j]] <- as.character(selected.cells[[group.name]][[j]])
    }
  }

  x <- assay(sce)
  # error that type 'S4' is non subsettable when trying to fetch data using tatami from a DelayedArray
  # so converting to sparseMatrix as a temporary work around until issue with seed handling in tatami_r fixed
  if(!is(x,'sparseMatrix')) {
    x <- as(x, 'sparseMatrix')
  }
  ptr <- beachmat::initializeCpp(x)
  replicates.matrix <- matrix(NA, nrow = nrow(x), ncol = rep.num * length(selected.cells))
  replicates.matrix.coldata <- data.frame(
    group = rep(names(selected.cells),each=rep.num),
    rep = rep(1:rep.num, length(selected.cells)),
    coverage.file = rep(NA,ncol(replicates.matrix))
  )
  replicates.matrix.coldata$coverage.file <- paste0(replicates.matrix.coldata$group,'_pseudobulk_',replicates.matrix.coldata$rep,'.txt')
  
  for(i in 1:nrow(replicates.matrix.coldata)) {
      output.file <- replicates.matrix.coldata$coverage.file[i]
      pseudobulk.cells <- unique(selected.cells[[replicates.matrix.coldata$group[i]]][[replicates.matrix.coldata$rep[i]]])
      #create_pseuobulk_file(mae$fragment_file,output.file,sub('.*#','',pseudobulk.cells))
      ptr.subset <- apply_subset(ptr, which(colnames(sce) %in% pseudobulk.cells), FALSE)
      replicates.matrix[,i] <- aggregate_counts(ptr.subset, rep(1,length(pseudobulk.cells))-1L, nthreads = 1)
  }
  replicates.matrix <- as(replicates.matrix, 'sparseMatrix')

  # Add replicates matrix as SummarizedExperiment to MAE
  
  # save the new pseudobulk replicate selection to the MAE
  for(i in names(selected.cells)) {
    for(j in 1:length(selected.cells[[i]])) {
      rep.colname <- paste0('Psuedobulk_Rep',j)
      if(is.null(alt.exp.name)) {
        if(! rep.colname %in% colnames(colData(mae[[sce.idx]]))) {
          colData(mae[[sce.idx]])[,rep.colname] <- FALSE
        }
        colData(mae[[sce.idx]])[selected.cells[[i]][[j]],rep.colname] <- TRUE
      } else {
        if(! rep.colname %in% colnames(colData(altExp(mae[[sce.idx]], alt.exp.name)))) {
          colData(altExp(mae[[sce.idx]], alt.exp.name))[,rep.colname] <- FALSE
        }
        colData(altExp(mae[[sce.idx]], alt.exp.name))[selected.cells[[i]][[j]],rep.colname] <- TRUE
      }
    }
  }
  
  mae
}

#####################################################################################################
# Write Coverage To Bed File for MACS2
#####################################################################################################

.writeCoverageToBed <- function(coverageFile = NULL, out = NULL, excludeChr = NULL, logFile = NULL){
  rmf <- .suppressAll(file.remove(out))
  allChr <- .availableSeqnames(coverageFile, "Coverage")
  if(!is.null(excludeChr)){
    allChr <- allChr[allChr %ni% excludeChr]
  }
  if(length(allChr)==0){
    stop("No Chromosomes in Coverage after Excluding Chr!")
  }
  ##Note that there was a bug with data.table vs data.frame with stacking
  for(x in seq_along(allChr)){
      iS <- .getCoverageInsertionSites(coverageFile = coverageFile, chr = allChr[x])
      iS <- data.table(seqnames = allChr[x], start = iS - 1L, end = iS)
      if(!any(is.na(iS$start))) {
        data.table::fwrite(iS, out, sep = "\t", col.names = FALSE, append = TRUE)
      } else {
        message(paste0("Warning - No insertions found on seqnames ", allChr[x], " for coverageFile ", coverageFile,"."))
      }
  }
  out
}
  