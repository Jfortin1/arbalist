#' Create Pseudobulk replicates
#'
#' Select cells for pseudobulk replicates, create coverage files and add a pseudobulk experiment to the MultiAssayExperiment.
#'
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param experiment.name String containing experiment name for creating the pseudobulk counts.
#' @param group.by.colname String containing experiment colData colname for the cell groups for pseudobulks.
#' @param sampleLabels String containing experiment colData colname for sample labels.
#' @param min.cells Integer scalar specifying minimum number of cells to select per group for pseudobulk creation.
#' @param max.cells Integer scalar specifying maximum number of cells to select per group for pseudobulk creation.
#' @param max.frags Integer scalar specifying maximum number of fragments per cell group to add to the pseudobulk coverage file.
#' @param min.reps Integer scalar specifying minimum number of pseudobulk replicates to generate.
#' @param max.reps Integer scalar specifying maximum number of pseudobulk replicates to generate.
#' @param coverage.file.path String containing the output directory for generated coverage files.
#' @param skip.se.creation Logical specifying whether to skip generating the pseudobulk summarized experiment.
#' @param skip.coverage.file.creation Logical specifying whether to skip generating coverage files.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @return \linkS4class{MultiAssayExperiment} with a new pseudobulk \linkS4class{SummarizedExperiment} and
#' columns are added to the \linkS4class{SingleCellExperiment} colData specifying which pseudobulk cell composition.
#' Lastly, pseudobulk coverage files containing the fragments coordinates are written to file and file names are specified in the pseudobulk \linkS4class{SummarizedExperiment} colData.
#' 
#' @author Natalie Fox
#' @importFrom BiocParallel bptry bplapply bpparam
#' @importFrom SummarizedExperiment colData rowRanges
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList colData colData<- listToMap experiments
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom beachmat initializeCpp flushMemoryCache
#' @export
addGroupCoverages <- function(
    mae,
    experiment.name = "TileMatrix500",
    group.by.colname = "Clusters",
    sampleLabels = "Sample",
    min.cells = 40,
    max.cells = 500,
    max.frags = 25*10^6,
    min.reps = 2,
    max.reps = 5,
    coverage.file.path = getwd(),
    skip.se.creation = FALSE,
    skip.coverage.file.creation = FALSE,
    BPPARAM = bpparam()
){

  # Find the experiment result
  sce.list <- findSCE(mae,experiment.name)
  sce <- sce.list$sce
  sce.idx <- sce.list$sce.idx
  alt.exp.name <- sce.list$alt.exp.name

  # count the number of cells per sample and per group
  num.samples <- length(unique(colData(sce)[,sampleLabels]))
  num.cells.per.group <- table(colData(sce)[,group.by.colname])
  
  # find the group names with at least min.cells
  num.cells.per.group <- table(colData(sce)[,group.by.colname])[which(num.cells.per.group >= min.cells)]
  
  # decide the number of cells to sample per group
  num.cells.to.sample <- max(min(min(num.cells.per.group),max.cells),min.cells)
  
  # find the number of cells to sample per group and sample
  num.cells.to.sample.from.group.and.sample <- ceiling(num.cells.to.sample/num.samples)
  num.cells.per.group.and.sample <- table(paste0(colData(sce)[,group.by.colname],colData(sce)[,sampleLabels])[colData(sce)[,group.by.colname] %in% names(num.cells.per.group)])
  overall.num.cells.to.sample.from.group.and.sample <- max(min(min(num.cells.per.group.and.sample),num.cells.to.sample.from.group.and.sample),ceiling(min.cells/num.samples))
  
  # decide the number of pseudobulk reps to create
  rep.num <- max( min.reps, min( max.reps, floor(min(num.cells.per.group)/num.cells.to.sample)))

  # list the samples per group and per group and sample
  cells.per.group <- split(
    colnames(sce)[colData(sce)[,group.by.colname] %in% names(num.cells.per.group)],
    colData(sce)[colData(sce)[,group.by.colname] %in% names(num.cells.per.group),group.by.colname])
  cells.per.group.and.sample <- sapply(cells.per.group,function(x) {split(x,sub('#.*','',x))},simplify = FALSE)
  
  selected.cells <- list()
  for(group.name in names(cells.per.group.and.sample)) {

    # reset in case the number of cells to sample was change for the last group
    num.cells.to.sample.from.group.and.sample <- overall.num.cells.to.sample.from.group.and.sample
    
    # get the cells to sample
    x <- cells.per.group.and.sample[[group.name]]
    
    selected.cells[[group.name]] <- list()
    
    # If there aren't enough cells from individual samples within the group, merge samples for sampling
    x.per.sample.length <- vapply(x,length,FUN.VALUE = 1L)
    if(any(x.per.sample.length < num.cells.to.sample.from.group.and.sample)) {
      x.merged <- x[which(x.per.sample.length >= num.cells.to.sample.from.group.and.sample)]
      x.merged[[length(x.merged)+1]] <- unlist(x[which(x.per.sample.length < num.cells.to.sample.from.group.and.sample)])
      x <- x.merged
      if(any(vapply(x,length,FUN.VALUE = 1L) < num.cells.to.sample.from.group.and.sample)) {
        x <- list(unlist(x))
      }
      num.cells.to.sample.from.group.and.sample <- ceiling(min.cells/length(x))
    }
    
    # for each sample, select cells for the pseudobulk replicate
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
    # make sure the results are character vectors
    for(j in 1:rep.num) {
      selected.cells[[group.name]][[j]] <- as.character(selected.cells[[group.name]][[j]])
    }
  }
  
  # save the new pseudobulk replicate selection to the MAE
  for(i in names(selected.cells)) {
    for(j in 1:length(selected.cells[[i]])) {
      rep.colname <- paste0('Pseudobulk_Rep',j)
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

  x <- assay(sce)

  # set up metadata for the pseudobulk replicates
  replicates.matrix <- matrix(NA, nrow = nrow(x), ncol = rep.num * length(selected.cells))
  replicates.matrix.coldata <- data.frame(
    group = rep(names(selected.cells),each=rep.num),
    rep = rep(1:rep.num, length(selected.cells)),
    num.cells = rep(NA,ncol(replicates.matrix)),
    coverage.file = rep(NA,ncol(replicates.matrix))
  )
  replicates.matrix.coldata$ID <- paste0(replicates.matrix.coldata$group,'_pseudobulk_rep',replicates.matrix.coldata$rep)
  if(is.na(coverage.file.path) || coverage.file.path == '') {
    replicates.matrix.coldata$coverage.file <- paste0(replicates.matrix.coldata$ID,'.txt')
  } else {
    replicates.matrix.coldata$coverage.file <- paste0(coverage.file.path,'/',replicates.matrix.coldata$ID,'.txt')
  }

  # create coverage files
  .create_pseudobulk_file_parallel_func <- function(j, replicates.matrix.coldata, mae) {
    output.file <- replicates.matrix.coldata$coverage.file[j]
    pseudobulk.cells <- unique(selected.cells[[replicates.matrix.coldata$group[j]]][[replicates.matrix.coldata$rep[j]]])
    create_pseudobulk_file(mae$fragment_file, output.file, sub('.*#','',pseudobulk.cells))
  }
  if(!skip.coverage.file.creation) {
    res.list <- bptry(bplapply(seq(nrow(replicates.matrix.coldata)), .create_pseudobulk_file_parallel_func, replicates.matrix.coldata=replicates.matrix.coldata, mae=mae, BPPARAM = BPPARAM))
  }

  if(!skip.se.creation) {
    beachmat::flushMemoryCache()
    ptr <- beachmat::initializeCpp(x, memorize=FALSE)
    # Create matrix for pseudobulk experiment
    for(i in 1:nrow(replicates.matrix.coldata)) {
      output.file <- replicates.matrix.coldata$coverage.file[i]
      pseudobulk.cells <- unique(selected.cells[[replicates.matrix.coldata$group[i]]][[replicates.matrix.coldata$rep[i]]])
      ptr.subset <- apply_subset(ptr, which(colnames(sce) %in% pseudobulk.cells), FALSE)
      replicates.matrix[,i] <- aggregate_counts(ptr.subset, rep(1,length(pseudobulk.cells))-1L, nthreads = 1, binarize = FALSE)
      replicates.matrix.coldata$num.cells <- length(pseudobulk.cells)
    }
    replicates.matrix <- as(replicates.matrix, 'sparseMatrix')
    
    # Add replicates matrix as SummarizedExperiment to MAE
    pseudobulk.se <-  SummarizedExperiment(list(counts=replicates.matrix))
    rowRanges(pseudobulk.se) <- rowRanges(sce)
    colData(pseudobulk.se) <- DataFrame(replicates.matrix.coldata)
    colnames(pseudobulk.se) <- replicates.matrix.coldata$ID
    
    exp.list <- experiments(mae)
    exp.list[[paste0(experiment.name,'_pseudobulk')]] <- pseudobulk.se
    
    el <- ExperimentList(exp.list)
    maplist <- lapply(exp.list, function(se) {
      if(sampleLabels %in% colnames(colData(se))) {
        data.frame(primary = colData(se)[,sampleLabels], colname = colnames(se), stringsAsFactors = FALSE)
      } else {
        data.frame(primary = se$ID, colname = colnames(se), stringsAsFactors = FALSE)
      }
    })
    sampMap <- listToMap(maplist)
    
    # Create and annotate the MultiAssayExperiment
    new.mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
    pseudobulk.colData.filler <- matrix(NA,ncol=ncol(colData(mae)),nrow=ncol(pseudobulk.se))
    colnames(pseudobulk.colData.filler) <- colnames(colData(mae))
    rownames(pseudobulk.colData.filler) <- colnames(pseudobulk.se)
    
    colData(new.mae) <- rbind(colData(mae),pseudobulk.colData.filler)
    metadata(new.mae) <- metadata(mae)
    
    beachmat::flushMemoryCache()
    
    return(new.mae)
  } else {
    return(mae)
  }
}