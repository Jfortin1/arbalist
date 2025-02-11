#' Create Pseudobulk replicates
#'
#' Select cells for pseudobulk replicates, create coverage files and add a pseudobulk experiment to the MultiAssayExperiment.
#'
#' @param mae \linkS4class{MultiAssayExperiment}.
#' @param experiment.name String containing experiment name for creating the pseudobulk counts.
#' @param group.by.colname String containing experiment colData colname for the cell groups for pseudobulks.
#' @param sample.colname String containing experiment colData colname for sample labels.
#' @param min.cells Integer scalar specifying minimum number of cells to select per group for pseudobulk creation.
#' @param max.cells Integer scalar specifying maximum number of cells to select per group for pseudobulk creation.
#' @param max.frags Integer scalar specifying maximum number of fragments per cell group to add to the pseudobulk coverage file.
#' @param min.reps Integer scalar specifying minimum number of pseudobulk replicates to generate.
#' @param max.reps Integer scalar specifying maximum number of pseudobulk replicates to generate.
#' @param sampling.ratio Numeric scalar specifying the fraction of the total cells that can be sampled to generate any given pseudobulk replicate.
#' @param coverage.file.path String containing the output directory for generated coverage files.
#' @param skip.se.creation Logical specifying whether to skip generating the pseudobulk summarized experiment.
#' @param skip.coverage.file.creation Logical specifying whether to skip generating coverage files.
#' @param num.threads Integer scalar specifying the number of threads.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @return \linkS4class{MultiAssayExperiment} with a new pseudobulk \linkS4class{SummarizedExperiment} and
#' columns are added to the \linkS4class{SingleCellExperiment} colData specifying which pseudobulk cell composition.
#' Lastly, pseudobulk coverage files containing the fragments coordinates are written to file and file names are specified in the pseudobulk \linkS4class{SummarizedExperiment} colData.
#' 
#' @author Natalie Fox
#' @importFrom BiocParallel bptry bplapply bpparam
#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<- rowRanges rowRanges<-
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList listToMap experiments
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom SingleCellExperiment altExp altExp<-
#' @importFrom beachmat initializeCpp flushMemoryCache tatami.subset tatami.row.sums
#' @export
addGroupCoverages <- function(
    mae,
    experiment.name = "TileMatrix500",
    group.by.colname = "Clusters",
    sample.colname = "Sample",
    min.cells = 40,
    max.cells = 500,
    max.frags = 25*10^6,
    min.reps = 2,
    max.reps = 5,
    sampling.ratio = 0.8,
    coverage.file.path = getwd(),
    skip.se.creation = FALSE,
    skip.coverage.file.creation = FALSE,
    num.threads = 1,
    BPPARAM = bpparam()
){
  # Find the experiment result
  sce.list <- findSCE(mae,experiment.name)
  sce <- sce.list$sce
  sce.idx <- sce.list$sce.idx
  alt.exp.name <- sce.list$alt.exp.name
  
  selected.cells <- selectCellsForPseudobulks(
    sce,
    group.by.colname = group.by.colname,
    sample.colname = sample.colname,
    min.cells = min.cells,
    max.cells = max.cells,
    min.reps = min.reps,
    max.reps = max.reps,
    sampling.ratio = sampling.ratio
  )

  # save the new pseudobulk replicate selection to the MAE
  max.rep.num <- 1
  for(i in names(selected.cells)) {
    if(length(selected.cells[[i]]) > max.rep.num) {
      max.rep.num <- length(selected.cells[[i]])
    }
    for(j in seq_len(length(selected.cells[[i]]))) {
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
  replicates.matrix <- matrix(NA, nrow = nrow(x), ncol = max.rep.num * length(selected.cells))
  replicates.matrix.coldata <- data.frame(
    group = rep(names(selected.cells),each=max.rep.num),
    rep = rep(seq_len(max.rep.num), length(selected.cells)),
    num.cells = rep(NA,ncol(replicates.matrix)),
    coverage.file = rep(NA,ncol(replicates.matrix)),
    check = rep(FALSE, ncol(replicates.matrix))
  )
  replicates.matrix.coldata$ID <- paste0(replicates.matrix.coldata$group,'_pseudobulk_rep',replicates.matrix.coldata$rep)
  if(is.na(coverage.file.path) || coverage.file.path == '') {
    replicates.matrix.coldata$coverage.file <- paste0(replicates.matrix.coldata$ID,'.txt')
  } else {
    replicates.matrix.coldata$coverage.file <- paste0(coverage.file.path,'/',replicates.matrix.coldata$ID,'.txt')
  }
  for(i in seq_len(length(selected.cells))) {
    replicates.matrix.coldata$check[replicates.matrix.coldata$group == names(selected.cells)[i] & replicates.matrix.coldata$rep %in% seq_along(selected.cells[[i]])] <- TRUE
  }
  replicates.matrix <- replicates.matrix[,replicates.matrix.coldata$check]
  replicates.matrix.coldata <- replicates.matrix.coldata[replicates.matrix.coldata$check,colnames(replicates.matrix.coldata) != 'check']

  # create coverage files
  .create.pseudobulk.file.parallel.func <- function(j, replicates.matrix.coldata, mae) {
    output.file <- replicates.matrix.coldata$coverage.file[j]
    pseudobulk.cells <- unique(selected.cells[[replicates.matrix.coldata$group[j]]][[replicates.matrix.coldata$rep[j]]])
    create_pseudobulk_file(mae$fragment_file[unique(colData(sce)[pseudobulk.cells,sample.colname])], output.file, sub('.*#','',pseudobulk.cells))
  }
  if(!skip.coverage.file.creation) {
    res.list <- bptry(bplapply(seq_len(nrow(replicates.matrix.coldata)), .create.pseudobulk.file.parallel.func, replicates.matrix.coldata=replicates.matrix.coldata, mae=mae, BPPARAM = BPPARAM))
  }

  if(!skip.se.creation) {
    beachmat::flushMemoryCache()
    if(is(x,"DelayedArray")) {
      ptr <- beachmat.hdf5::initializeCpp(x, memorize=TRUE)
    } else {
      ptr <- beachmat::initializeCpp(x, memorize=TRUE)
    }
    # Create matrix for pseudobulk experiment
    for(i in seq_len(nrow(replicates.matrix.coldata))) {
      output.file <- replicates.matrix.coldata$coverage.file[i]
      pseudobulk.cells <- unique(selected.cells[[replicates.matrix.coldata$group[i]]][[replicates.matrix.coldata$rep[i]]])
      ptr.subset <- tatami.subset(ptr, subset = which(colnames(sce) %in% pseudobulk.cells), by.row = FALSE)
      replicates.matrix[,i] <- tatami.row.sums(ptr.subset, num.threads = 1)
      replicates.matrix.coldata$num.cells <- length(pseudobulk.cells)
    }
    replicates.matrix <- as(replicates.matrix, 'sparseMatrix')
    
    # Add replicates matrix as SummarizedExperiment to MAE
    pseudobulk.se <-  SummarizedExperiment(list(counts = replicates.matrix))
    rowRanges(pseudobulk.se) <- rowRanges(sce)
    colData(pseudobulk.se) <- DataFrame(replicates.matrix.coldata)
    colnames(pseudobulk.se) <- replicates.matrix.coldata$ID
    
    exp.list <- experiments(mae)
    exp.list[[paste0(experiment.name,'_pseudobulk')]] <- pseudobulk.se
    
    el <- ExperimentList(exp.list)
    maplist <- lapply(exp.list, function(se) {
      if(sample.colname %in% colnames(colData(se))) {
        data.frame(primary = colData(se)[,sample.colname], colname = colnames(se), stringsAsFactors = FALSE)
      } else {
        data.frame(primary = se$ID, colname = colnames(se), stringsAsFactors = FALSE)
      }
    })
    sampMap <- listToMap(maplist)
    
    # Create and annotate the MultiAssayExperiment
    new.mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
    pseudobulk.colData.filler <- matrix(NA, ncol = ncol(colData(mae)), nrow = ncol(pseudobulk.se))
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

selectCellsForPseudobulks <- function(
    sce,
    group.by.colname = "Clusters",
    sample.colname = "Sample",
    min.cells = 40,
    max.cells = 500,
    min.reps = 2,
    max.reps = 5,
    sampling.ratio = 0.8
){
  
  # count the number of cells per sample and per group
  num.samples <- length(unique(colData(sce)[,sample.colname]))
  num.cells.per.group <- table(colData(sce)[,group.by.colname])
  
  if(is.null(max.reps)) {
    max.reps <- num.samples
  }
  
  sampling.info <- data.frame(cell = colnames(sce), sample = colData(sce)[,sample.colname], group = colData(sce)[,group.by.colname])
  sample.group.counts <- table(sampling.info[,-1])
  
  # find the group names with at least min.cells
  num.cells.per.group <-  num.cells.per.group[which(num.cells.per.group >= min.cells)]
  
  selected.cells <- list()
  unused.cells <- list()
  
  # select the groups with enough cells per sample to meet the min reps with min cells sampled without replacement
  for(group.name in colnames(sample.group.counts)[apply(sample.group.counts > min.cells, 2, sum) >= min.reps]) {
    group.cells <- sampling.info[sampling.info$group == group.name,]
    selected.cells[[group.name]] <- list()
    for(sample.name in names(sort(table(group.cells$sample),decreasing = TRUE))[seq_len(min(max.reps,sum(table(group.cells$sample) >= min.cells)))]) {
      if(max.cells < sum(group.cells$sample == sample.name)) {
        selected.cells[[group.name]][[sample.name]] <- sample(group.cells$cell[group.cells$sample == sample.name], max.cells)
      } else {
        selected.cells[[group.name]][[sample.name]] <- group.cells$cell[group.cells$sample == sample.name]
      }
    }
  }
  # select the groups with enough cells to meet the min reps with min cells sampled without replacement if we combine samples
  sample.group.counts <- sample.group.counts[,setdiff(colnames(sample.group.counts),names(selected.cells)), drop=FALSE]
  for(group.name in colnames(sample.group.counts)[colSums((sample.group.counts <= 40)*sample.group.counts + min.cells*(sample.group.counts > min.cells)) > min.cells*min.reps]) {
    group.cells <- sampling.info[sampling.info$group == group.name,]
    selected.cells[[group.name]] <- list()
    unused.cells[[group.name]] <- c()
    for(sample.name in names(sort(table(group.cells$sample),decreasing = TRUE))) {
      if(length(selected.cells[[group.name]]) >= max.reps) {
        # if we have enough replicates, then move on to the next group
        break
      }
      if(sum(group.cells$sample == sample.name) > min.cells) {
        # this sample has enough cells to be a replicate on its own
        if(min.cells < sum(group.cells$sample == sample.name)) {
          selected.cells[[group.name]][[sample.name]] <- sample(group.cells$cell[group.cells$sample == sample.name], min.cells)
          unused.cells[[group.name]] <- c(unused.cells[[group.name]],setdiff(group.cells$cell[group.cells$sample == sample.name],selected.cells[[group.name]][[sample.name]]))
        } else {
          selected.cells[[group.name]][[sample.name]] <- group.cells$cell[group.cells$sample == sample.name]
        }
      } else {
        # need to combine samples to create the replicate
        last.rep.num <- length(selected.cells[[group.name]])
        if(last.rep.num == 0 || length(selected.cells[[group.name]][[last.rep.num]]) >= min.cells) {
          # start a new replicate
          selected.cells[[group.name]][[sample.name]] <- group.cells$cell[group.cells$sample == sample.name]
        } else {
          # check if we need to sample or if all cells are required in the replicate
          num.cell.in.rep.already <- length(selected.cells[[group.name]][[last.rep.num]])
          if(min.cells < (sum(group.cells$sample == sample.name) + num.cell.in.rep.already)) {
            selected.cells[[group.name]][[last.rep.num]] <- c(selected.cells[[group.name]][[last.rep.num]],sample(group.cells$cell[group.cells$sample == sample.name], min.cells - num.cell.in.rep.already))
            unused.cells[[group.name]] <- c(unused.cells[[group.name]],setdiff(group.cells$cell[group.cells$sample == sample.name],selected.cells[[group.name]][[sample.name]]))
          } else {
            selected.cells[[group.name]][[last.rep.num]] <- c(selected.cells[[group.name]][[last.rep.num]],group.cells$cell[group.cells$sample == sample.name])
          }
          names(selected.cells[[group.name]])[last.rep.num] <- paste0(names(selected.cells[[group.name]])[last.rep.num],'_',sample.name)
        }
      }
    }
    if(length(selected.cells[[group.name]][[length(selected.cells[[group.name]])]]) < min.cells) {
      last.rep.num <- length(selected.cells[[group.name]])
      if(length(unused.cells[[group.name]]) > (min.cells - length(selected.cells[[group.name]][[last.rep.num]]))) {
        selected.cells[[group.name]][[last.rep.num]] <- c(selected.cells[[group.name]][[last.rep.num]], sample(unused.cells[[group.name]], min.cells - length(selected.cells[[group.name]][[last.rep.num]])))
#      } else if(length(unused.cells[[group.name]]) > (min.cells - length(selected.cells[[group.name]][[last.rep.num]]))) {
#        selected.cells[[group.name]][[last.rep.num]] <- c(selected.cells[[group.name]][[last.rep.num]], unused.cells[[group.name]])
      } else {
        selected.cells[[group.name]] <- selected.cells[[group.name]][-last.rep.num]
      }
    }
  }
  sample.group.counts <- sample.group.counts[,setdiff(colnames(sample.group.counts),names(selected.cells)), drop=FALSE]
  # use different cells but from the same sample to complete reps (sampling without replacement)
  for(group.name in colnames(sample.group.counts)[colSums(sample.group.counts) >= min.cells * min.reps]) {
    group.cells <- sampling.info[sampling.info$group == group.name,]
    group.cells$sample <- factor(group.cells$sample, levels=names(sort(table(group.cells$sample), decreasing = TRUE)))
    group.cells <- group.cells[order(group.cells$sample),]
    selected.cells[[group.name]] <- list()
    num.cells.to.sample <- min(max.cells,max(min.cells,floor(nrow(group.cells)/min.reps)))
    for(rep.num in seq_len(min(max.reps,max(sum(table(group.cells$sample) >= min.cells), min.reps)))) {
      if(num.cells.to.sample < nrow(group.cells)) {
        selected.cells[[group.name]][[rep.num]] <- sample(group.cells$cell, num.cells.to.sample)
        group.cells <- group.cells[!group.cells$cell %in% selected.cells[[group.name]][[rep.num]],]
      } else {
        selected.cells[[group.name]][[rep.num]] <- group.cells$cell
      }
    }
  }
  sample.group.counts <- sample.group.counts[,setdiff(colnames(sample.group.counts),names(selected.cells)), drop=FALSE]
  # use replacement to create min.reps with min.cells as long as we have enough cells adjusting with the sampling ratio
  for(group.name in colnames(sample.group.counts)[colSums(sample.group.counts) >= (min.cells * min.reps * sampling.ratio + min.cells - min.cells * sampling.ratio)]) {
    group.cells <- sampling.info[sampling.info$group == group.name,]
    group.cells$sample <- factor(group.cells$sample, levels=names(sort(table(group.cells$sample), decreasing = TRUE)))
    group.cells <- group.cells[order(group.cells$sample),]
    selected.cells[[group.name]] <- list()
    # select the sampling ratio amount of cells without replacement
    for(rep.num in seq_len(min.reps)) {
      selected.cells[[group.name]][[rep.num]] <- sample(group.cells$cell, floor(min.cells*sampling.ratio))
      group.cells <- group.cells[!group.cells$cell %in% selected.cells[[group.name]][[rep.num]],]
    }
    # fill in the rest of the replicates using the same cells for all replicates so there will be overlap between replicates
    for(rep.num in seq_len(min.reps)) {
      selected.cells[[group.name]][[rep.num]] <- c(selected.cells[[group.name]][[rep.num]],sample(group.cells$cell, min.cells - length(selected.cells[[group.name]][[rep.num]])))
    }
  }
  sample.group.counts <- sample.group.counts[,setdiff(colnames(sample.group.counts),names(selected.cells)), drop=FALSE]
  for(group.name in colnames(sample.group.counts)[colSums(sample.group.counts) >= min.cells * min.reps * sampling.ratio]) {
    group.cells <- sampling.info[sampling.info$group == group.name,]
    group.cells$sample <- factor(group.cells$sample, levels=names(sort(table(group.cells$sample), decreasing = TRUE)))
    group.cells <- group.cells[order(group.cells$sample),]
    selected.cells[[group.name]] <- list()
    # select the sampling ratio amount of cells without replacement
    for(rep.num in seq_len(min.reps)) {
      selected.cells[[group.name]][[rep.num]] <- sample(group.cells$cell, min.cells)
    }
  }
  sample.group.counts <- sample.group.counts[,setdiff(colnames(sample.group.counts),names(selected.cells)), drop=FALSE]
  
  if(ncol(sample.group.counts) > 0) {
    warning("Because there wasn't enough cells, no pseudobulk was created for ",paste(colnames(sample.group.counts), collapse=', '))
  }
  
  return(selected.cells)
}