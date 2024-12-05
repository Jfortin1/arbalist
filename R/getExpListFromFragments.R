#' Create a list of single cell experiments
#' 
#' Create a list of single cell experiments from fragment files
#' 
#' @param fragment.files Vector of strings specifying fragment files. Vector names need to be sample names.
#' @param output.dir String containing the directory where files should be output while creating the \linkS4class{MultiAssayExperiment}.
#' @param gene.grs Genomic Ranges specifying gene coordinates for creating the gene score matrix. If NA, the gene accessibility matrix will not be created.
#' @param barcodes.list A List with samples as names and the values a vector of barcodes for that sample. If NULL, all barcodes from the fragment file will be used.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.
#' @inheritParams saveTileMatrix
#' 
#' @return A list of experiments
#' 
#' @author Natalie Fox
#' @importFrom SingleCellExperiment mainExpName<-
#' @importFrom BiocParallel bptry bplapply bpparam
getExpListFromFragments <- function(
  fragment.files,
  output.dir = tempdir(),
  tile.size = 500,
  seq.lengths = NULL,
  gene.grs = NULL,
  barcodes.list = NULL,
  BPPARAM = bpparam()
) {
  
  if(is.null(names(fragment.files))) {
    stop('Please add sample names as the names on fragment.files')
  }
  
  all.exp <- list()
  
  # check if the hdf5 files already exist
  matrix.name <- paste0('TileMatrix',tile.size)
  output.file.names <- paste0(output.dir,'/',matrix.name,'_',names(fragment.files),'.h5')
  if(any(table(output.file.names) > 1)) {
    # if two file names are the same then add a random component to the file name.
    output.file.names <- tempfile(pattern=paste0(matrix.name,'_',names(fragment.files),'_'),tmpdir = output.dir, fileext='.h5')
  }
  names(output.file.names) <- names(fragment.files)
  if(any(file.exists(output.file.names))) {
    stop(paste0(output.file.names[which(file.exists(output.file.names))[1]],' already exists. We do not want to overwrite the file in case it is being used. Either remove the file if you think it is safe to do so or specify a different output.dir.'))
  }
  
  # check that fragment headers contained required information
  if(is.null(seq.lengths)) {
    info <- .processFragmentHeader(fragment.files[1])
    if(! 'reference_path' %in% names(info)) {
      stop('the fragment file header does not have reference_path information so please specify the seq.lengths argument')
    }
  }
  
  # For each fragment.file create a tile matrix hdf5 file
  # Parallelizing per sample
  res.list <- bptry(bplapply(seq_along(fragment.files), .saveTileMatrixCall, fragment.files=fragment.files, output.file.names=output.file.names, tile.size=tile.size, seq.lengths=seq.lengths, barcodes.list=barcodes.list, BPPARAM = BPPARAM))
  tile.res.list <- sapply(res.list,function(x){x$counts})
  names(tile.res.list) <- names(fragment.files)
  tile.grs <- res.list[[1]]$tiles
  # check that the tiles are the same for all samples
  for(i in setdiff(seq_along(res.list),1)) {
    if(length(tile.grs) != length(res.list[[i]]$tiles) || !all(tile.grs == res.list[[i]]$tiles)) {
      stop('Tile Matrix GRanges do not match')
    }
  }
  all.exp[[matrix.name]] <- .getSCEFromH5List(tile.res.list, tile.grs)

  if(!is.null(gene.grs)) {
    matrix.name <- 'GeneAccessibilityMatrix'
    
    # check if the hdf5 files already exist
    output.file.names <- paste0(output.dir,'/',matrix.name,'_',names(fragment.files),'.h5')
    names(output.file.names) <- names(fragment.files)
    if(any(table(output.file.names) > 1)) {
      # if two file names are the same then add a random component to the file name.
      output.file.names <- tempfile(pattern=paste0(matrix.name,'_',names(fragment.files),'_'),tmpdir = output.dir, fileext='.h5')
    }
    names(output.file.names) <- names(fragment.files)
    if(any(file.exists(output.file.names))) {
      stop(paste0(output.file.names[which(file.exists(output.file.names))[1]],' already exists. We do not want to overwrite the file in case it is being used. Either remove the file if you think it is safe to do so or specify a different output.dir.'))
    }

    # For each fragment.file create a gene score matrix hdf5 file
    # Parallelizing per sample
    gs.res.list <- bptry(bplapply(seq_along(fragment.files), .saveRegionMatrixCall, fragment.files=fragment.files, output.file.names=output.file.names, regions=gene.grs, barcodes.list=barcodes.list, BPPARAM = BPPARAM))
    names(gs.res.list) <- names(fragment.files)
    
    # Create a SingleCellExperiment for the Gene regions (ATAC-seq results)
    all.exp[[matrix.name]] <- .getSCEFromH5List(gs.res.list, gene.grs)
  }
  
  # Fill in the experiment name on each SingleCellExperiment
  for(exp.name in names(all.exp)) {
    mainExpName(all.exp[[exp.name]]) <- exp.name
  }
  
  return(all.exp)
}

#' @importFrom alabaster.matrix AmalgamatedArray
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges<- colData<-
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

.saveTileMatrixCall <- function(
  sample.name,
  fragment.files,
  output.file.names,
  tile.size = 500,
  seq.lengths = NULL,
  barcodes.list = NULL
) {
  tile.res <- saveTileMatrix(
    as.character(fragment.files[sample.name]),
    output.file = as.character(output.file.names[sample.name]),
    output.name = 'tile_matrix',
    tile.size = tile.size,
    seq.lengths = seq.lengths,
    barcodes = as.character(barcodes.list[[sample.name]])
  )
  return(tile.res)
}

.saveRegionMatrixCall <- function(
  sample.name,
  fragment.files,
  output.file.names,
  regions,
  barcodes.list = NULL
) {
  matrix.res <- saveRegionMatrix(
    as.character(fragment.files[sample.name]),
    output.file = as.character(output.file.names[sample.name]),
    output.name = 'gene_matrix',
    regions = regions,
    barcodes = as.character(barcodes.list[[sample.name]])
  )
  return(matrix.res)
}
