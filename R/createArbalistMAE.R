#' Import a single cell ATAC-seq MultiAssayExperiment
#' 
#' Import scATAC-seq or multiome results into a MultiAssayExperiment.
#' 
#' @param sample.names String vector containing the sample names. Must be the same length as fragment.files & filtered.feature.matrix.files.
#' @param fragment.files String vector containing the fragment file names and paths. ex. 10x cellranger result output atac_fragments.tsv.gz or fragments.tsv.gz
#' @param filtered.feature.matrix.files String vector containing the filtered feature matrix file names and paths. ex. 10x cellranger result output filtered_feature_bc_matrix.h5 or filtered_tf_bc_matrix.h5
#' @param barcode.annotation.files String vector containing barcode annotation to put in the SingleCellExperiment colData. ex. 10x cellranger result output per_barcode_metrics.csv or singlecell.csv
#' @param sample.annotation.files String vector containing sample annotation to put in the MultiAssayExperiment colData. ex. 10x cellranger result output summary.csv
#' @param output.dir String containing the directory where files should be output while creating the \linkS4class{MultiAssayExperiment}.
#' @param multiome Logical whether to use createMultiomeRNASCE on the filtered.feature.matrix.files to extrct the RNA features and create a SingleCellExperiment. If NULL, then will become TRUE if filtered.feature.matrix.files contain "filtered_feature_bc_matrix.h5" otherwise FALSE.
#' @param min.frags Number specifying the minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses. Cells containing greater than or equal to min.frags total fragments will be retained.
#' @param max.frags Number specifying the maximum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses. Cells containing less than or equal to max.frags total fragments will be retained.
#' @param tile.size Number specifying size of the tiles across the genome in base pairs.
#' @param seq.lengths Named integer vector containing the lengths of the reference sequences used for alignment.
#' @param gene.grs Genomic Ranges specifying gene coordinates for creating the gene score matrix. If NULL, then the geneset will be selected based on the genome version.
#' @param use.alt.exp Logical for selecting the MultiAssayExperiment structure. TRUE means that there will only be one experiment in the MultiAssayExperiment and all other experiments will be in alternative experiments. This option is only available if the columns are the same for all Matrices. FALSE means that each Matrix will be a separate experiment in the MAE.
#' @param main.exp.name String containing the name of the experiment that will be the main experiment when use.alt.exp is TRUE.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#'
#' @return A \linkS4class{MultiAssayExperiment}
#' 
#' @author Natalie Fox
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList colData colData<- listToMap
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment altExp<-
#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<-
#' @importFrom BiocParallel bpparam
#' @export
createArbalistMAE <- function(
  sample.names,
  fragment.files,
  filtered.feature.matrix.files,
  barcode.annotation.files = NULL,
  sample.annotation.files = NULL,
  output.dir = tempdir(),
  multiome = NULL,
  min.frags = 1000,
  max.frags = Inf,
  tile.size = 500,
  seq.lengths = NULL,
  gene.grs = NULL,
  use.alt.exp = FALSE,
  main.exp.name = 'TileMatrix500',
  BPPARAM = bpparam()
) {
  
  if(length(fragment.files) != length(sample.names)) {
    stop('fragment.files and sample names need to be the same length.')
  } else if(length(filtered.feature.matrix.files) != length(sample.names)) {
    stop('filtered.feature.matrix.files and sample names need to be the same length.')
  } else if(!is.null(barcode.annotation.files) && length(barcode.annotation.files) != length(sample.names)) {
    stop('barcode.annotation.files and sample names need to be the same length.')
  }
  
  names(fragment.files) <- sample.names
  
  # Get filtered barcodes and annotation for the barcodes
  barcodes.list <- list()
  barcode.anno.list <- list()
  barcode.anno <- NULL
  for(i in seq_along(filtered.feature.matrix.files)) {
    filtered.file <- filtered.feature.matrix.files[i]
    barcodes.list[[sample.names[i]]] <- h5read(filtered.feature.matrix.files[i],'/matrix/barcodes')
    barcode.anno.file <- barcode.annotation.files[i]
    if(!is.null(barcode.anno.file) && any(file.exists(barcode.anno.file))) {
      per.sample.barcodes <- read.csv(barcode.anno.file[file.exists(barcode.anno.file)])
      rownames(per.sample.barcodes) <- paste0(sample.names[i] ,'#', per.sample.barcodes$barcode)
      if(is.null(barcode.anno)) {
        barcode.anno <- per.sample.barcodes
      } else {
        barcode.anno <- rbind(barcode.anno, per.sample.barcodes)
      }
    }
  }
  
  if(is.null(names(fragment.files))) {
    stop('Please add sample names as the names on fragment.files')
  }
  
  # Get the ATAC-seq Experiment results from the fragment files
  all.exp <- .getExpListFromFragments(
    fragment.files,
    output.dir,
    tile.size,
    seq.lengths,
    gene.grs,
    barcodes.list = barcodes.list,
    BPPARAM = BPPARAM
  )

  # Add a SingleCellExperiment for RNA results if this is a multiome result
  if(is.null(multiome) && any(grep("filtered_feature_bc_matrix.h5",filtered.feature.matrix.files))) {
    multiome <- TRUE
  } else {
    multiome <- FALSE
  }
  if(multiome) {
    all.exp[['GeneExpressionMatrix']] <- createMultiomeRNASCE(
      h5.files = filtered.feature.matrix.files,
      sample.names = sample.names
    )
  }

  # If asked for, use alternative experiments for the MAE structure
  if(use.alt.exp) {
    if(!main.exp.name %in% names(all.exp) & length(exp.list) > 1) {
      cell.union <- colnames(all.exp[[1]])
      cell.intersect <- colnames(all.exp[[1]])
      for(i in 2:length(all.exp)) {
        cell.union <- union(cell.union, colnames(all.exp[[i]]))
        cell.intersect <- intersect(cell.intersect, colnames(all.exp[[i]]))
      }
      if(length(cell.union) == length(cell.intersect)) {
        # Convert the list of SCEs to main/alternative experiments in one SCE
        main.sce <- exp.list[[main.exp.name]]
        for(matrix.name in setdiff(names(exp.list), main.exp.name)) {
          altExp(main.sce, matrix.name) <- exp.list[[matrix.name]]
        }
        mainExpName(main.sce) <- main.exp.name
        if('GeneExpressionMatrix' %in% names(exp.list)) {
          exp.list <- list(multiome = main.sce)
        } else {
          exp.list <- list(atac = main.sce)
        }
      }
    }
  }
  
  # Create the MultiAssayExperiment from the list of Experiments
  mae <- .getMAEFromExpList(all.exp)
  
  # Add barcode information to the SE colData
  if(!is.null(barcode.anno)) {
    for(i in intersect(names(mae),c('TileMatrix500','GeneScoreMatrix','GeneAccessibilityMatrix'))) {
      if(any(grep('atac_',colnames(barcode.anno)))) {
        barcode.anno.atac.col <- barcode.anno[rownames(colData(mae[[i]])), setdiff(colnames(barcode.anno)[grep('atac_',colnames(barcode.anno))],c('atac_barcode','atac_peak_region_fragments','atac_peak_region_cutsites'))]
        colnames(barcode.anno.atac.col) <- sub('atac_','',colnames(barcode.anno.atac.col))
        colData(mae[[i]]) <- cbind(colData(mae[[i]]), barcode.anno.atac.col)
      } else {
        if(any(which(colnames(barcode.anno) == 'passed_filters')) && !any(which(colnames(barcode.anno) == 'fragments')))  {
          colnames(barcode.anno)[which(colnames(barcode.anno) == 'passed_filters')] <- 'fragments'
        }
        colData(mae[[i]]) <- cbind(colData(mae[[i]]), barcode.anno[rownames(colData(mae[[i]])),])
      }
    }
    for(i in intersect(names(mae),c('GeneExpressionMatrix'))) {
      barcode.anno.gex.col <- barcode.anno[rownames(colData(mae[[i]])), grep('gex_',colnames(barcode.anno))]
      colnames(barcode.anno.gex.col) <- sub('gex_','',colnames(barcode.anno.gex.col))
      colData(mae[[i]]) <- cbind(colData(mae[[i]]), barcode.anno.gex.col)
    }
  }
  
  # Add sample summary information to the MAE colData
  summary.info.all <- NULL
  for(i in seq_along(sample.annotation.files)) {
    summary.file.name <- sample.annotation.files[i]
    if(file.exists(summary.file.name)) {
      summary.info <- read.csv(summary.file.name)
      rownames(summary.info) <- summary.info$Sample.ID
      if(!is.null(summary.info.all) && ncol(summary.info.all) != ncol(summary.info)) {
        summary.info.all <- NULL
        break
      }
      summary.info.all <- rbind(summary.info.all, summary.info)
    }
  }
  if(!is.null(summary.info.all)) {
    colData(mae) <- cbind(colData(mae),summary.info.all[,!colnames(summary.info.all) %in% 'Sample.ID'])
  }
  # add fragment files to MAE colData
  colData(mae)$fragment_file <- fragment.files[rownames(colData(mae))]
  
  return(mae)
}

#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList colData colData<- listToMap
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment altExps
.getMAEFromExpList <- function(exp.list) {
  # Create the sample map for the MultiAssayExperiment
  el <- ExperimentList(exp.list)
  maplist <- lapply(exp.list, function(se) {
    data.frame(primary = se$Sample, colname = colnames(se), stringsAsFactors = FALSE)
  })
  sampMap <- listToMap(maplist)
  
  # Create and annotate the MultiAssayExperiment
  mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
  
  return(mae)
}