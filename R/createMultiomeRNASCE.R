#' Import gene expression results from multiome Cell Ranger results
#' 
#' Import RNA results from multiome Cell Ranger results into a SingleCellExperiment.
#' 
#' @param h5.files Vector of strings specifying filtered_feature_bc_matrix.h5 path. ex. could just be filtered_feature_bc_matrix.h5. Vector must be the same length as sample.names.
#' @param sample.names Vector of strings specifying sample names. Vector must be the same length as h5.files.
#' @param feature.type String specifying the feature type to select from filtered_feature_bc_matrix.h5.
#' @param filter.features.without.intervals Logical whether to remove features from the h5.files that do not have interval specified. Often these are mitochondria genes.
#'
#' @return A \linkS4class{SingleCellExperiment}
#' 
#' @author Natalie Fox
#' @importFrom BiocGenerics which
#' @importFrom GenomicRanges GRanges
#' @importFrom Matrix sparseMatrix
#' @importFrom rhdf5 h5read
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment mainExpName<-
#' @importFrom SummarizedExperiment SummarizedExperiment cbind rowData rowRanges<- 
#' @importFrom utils object.size read.csv
#' @export
createMultiomeRNASCE <- function(
  h5.files,
  sample.names = NULL,
  feature.type = 'Gene Expression',
  filter.features.without.intervals = FALSE
){

  # collect potential RNA files for each
  names(h5.files) <- sample.names

  se.all.samples <- NULL
  for(y in 1:length(h5.files)) {
    feature.matrix <- h5.files[y]
    sample.name <- sample.names[y]
    
    # load the Cell Ranger results (per sample matrix)
    barcodes <- h5read(feature.matrix, "/matrix/barcodes")
    data <- h5read(feature.matrix, "/matrix/data")
    indices <- h5read(feature.matrix, "/matrix/indices")
    indptr <- h5read(feature.matrix, "/matrix/indptr")
    shape <- h5read(feature.matrix, "/matrix/shape")
    features <- h5read(feature.matrix, "/matrix/features")
    
    sparse.feature.matrix <- sparseMatrix(
      i = indices, 
      p = indptr, 
      x = data, 
      dims = shape,
      index1 = FALSE
    )
    
    # combine sample name and barcode for cell identifiers/column names
    colnames(sparse.feature.matrix) <- paste0(sample.name, "#", barcodes)
    
    features <- as.data.frame(features[sapply(features,length) == nrow(sparse.feature.matrix)])
    
    # create a per sample SummarizedExperiment
    se <- SummarizedExperiment(assays = SimpleList(counts = sparse.feature.matrix), rowData = features)
    rownames(se) <- features$id
    
    # select the RNA features from the Cell Ranger results
    if('feature_type' %in% colnames(rowData(se))){
      se <- se[which(rowData(se)$feature_type == feature.type)]
    }

    gc()
    
    # combine sample results together
    if(is.null(se.all.samples)) {
      se.all.samples <- se
    } else {
      se.all.samples <- SummarizedExperiment::cbind(se.all.samples,se)
    }
  }

  colData(se.all.samples)$Sample <- sub('#.*$', '', colnames(se.all.samples))
  sce <- as(se.all.samples, 'SingleCellExperiment')
  
  # if intervals are specified in the feature information, then set the rowRanges for the experiment
  if(is.null(rowRanges(se)) && "interval" %in% colnames(rowData(se))) {
    na.features <- which(rowData(se)$interval == 'NA')
    if(!any(na.features)) {
      rowRanges(se) <- GRanges(rowData(se)$interval)
    } else if(filter.features.without.intervals) {
      warning(paste0(
        'Removing the following rows so that we can specify the rowRanges for all features: ',
        paste0(rowData(se)[na.features,"name"],collapse=', '), 
        ". Set the filter.features.without.intervals argument to FALSE to skip adding rowRanges."
        ))
      rowRanges(se) <- GRanges(rowData(se)$interval[rowData(se)$interval != 'NA'])
    }
  }
  
  mainExpName(sce) <- 'GeneExpressionMatrix'
  
  return(sce)
}
