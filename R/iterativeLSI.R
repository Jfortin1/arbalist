#' Calculate iterative LSI dimensionality reduction
#'
#' @param x input matrix for creating iterative LSI embeddings, assumes ATAC-seq so chooses top accessible features.
#' @param col.subset Integer vector specifying the subset of the columns to include in the iterative LSI calculations.
#' @param rank Integer scalar specifying the rank for irlba_realized.
#' @param iterations Integer scalar specifying number of LSI iterations to perform.
#' @param num.features Integer scalar specifying the number of accessible features to select when selecting the most accessible or most variable features.
#' @param lsi.method Number or string indicating the order of operations in the TF-IDF normalization. Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
#' @param cluster.method String containing cluster method. Current options: "seurat", "scran".
#' @param correlation.cutoff 	Numeric scalar specifying the cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the correlation.cutoff, it will be excluded from analysis.
#' @param scale.to Numeric scalar specifying the center for TF-IDF normalization.
#' @param num.threads Integer scalar specifying the number of threads to be used for parallel computing.
#' @param seed Numeric scalar to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param total.features Integer scalar specifying how many features to consider for use in LSI after ranking the features by the total number of insertions. These features are the only ones used throughout the variance identification and LSI.
#' @param filter.quantile Numeric scalar between 0 and 1 inclusive that indicates the quantile above which features should be removed based on insertion counts prior.
#' @param outlier.quantiles Numeric vector specifying the lower and upper quantiles for filtering cells based on the number of accessible regions.
#' @param binarize Logical specifying whether to binarize the matrix while creating iterative LSI reduced dimensions.
#' @return A list is returned containing:
#' \itemize{
#' \item \code{embedding}, a matrix containing the iterativeLSI embedding
#' \item \code{projection}, 
#' \item \code{subset}, a vector of indices specifying the selected subset of features
#' }
#' 
#' @author Natalie Fox
#' @export
#' @importFrom matrixStats rowVars
#' @importFrom beachmat initializeCpp flushMemoryCache
#' @importFrom stats kmeans cor
#' @importFrom utils head
iterativeLSI <- function(
    x,
    col.subset = NULL,
    rank = 30,
    iterations = 2,
    num.features = 25000,
    lsi.method = 1,
    cluster.method = "Seurat",
    correlation.cutoff = 0.75,
    scale.to = 10000,
    num.threads = 1,
    seed = 5,
    total.features = 500000,
    filter.quantile = 0.995,
    outlier.quantiles = c(0.02, 0.98),
    binarize = TRUE
) {
  
  # select the top features based on accessibility of the binarized features
  beachmat::flushMemoryCache()
  ptr <- beachmat::initializeCpp(x, memorize = TRUE)
  if(is.null(col.subset)) {
    col.subset <- seq(1, ncol(x))
  } else {
    col.subset <- sort(unique(col.subset),decreasing = TRUE)
    if(length(col.subset) <= ncol(x)) {
     ptr <- apply_subset(ptr, col.subset, row = FALSE)
    } else {
      stop('There are more entries in the column subset to select than the number of columns.')
    }
  }
  stats <- lsi_matrix_stats(ptr, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row
  row.accessibility <- stats$frequency
  rm.top <- floor((1-filter.quantile) * total.features)
  row.subset <- head(order(row.accessibility, decreasing = TRUE), num.features + rm.top )[-seq_len(rm.top)]
  
  # TF-IDF normalization (log(TF-IDF) method) and compute LSI
  lsi.res <- .computeLSI(x, ptr, col.subset = col.subset, row.subset = row.subset, lsi.method = lsi.method, scale.to = scale.to, num.dimensions = rank, outlier.quantiles = outlier.quantiles, seed = seed, binarize = binarize, num.threads = num.threads)
  embedding <- lsi.res$matSVD
  for (i in seq_len(iterations-1)) {
    # drop embeddings that are correlated with the library size
    embedding <- embedding[,which(cor(embedding, stats$sums[colnames(x)[col.subset] %in% rownames(embedding)]) <= correlation.cutoff),drop=FALSE]
    
    # find cell clusters
    cluster.output <- cluster.matrix(embedding, method = cluster.method)
    if(length(table(cluster.output)) == 1) {
      warning('Data is not splitting into clusters so we cannot calculate iterativeLSI')
      return(NULL)
    }
    # aggregate counts for each cluster and find the most variable features
    aggregated <- aggregate_counts(ptr, as.integer(factor(cluster.output)) - 1L, nthreads = num.threads, binarize = binarize)
    sums <- colSums(aggregated)
    center_sf <- function(y) {
      center <- if (is.null(scale.to)) mean(y) else scale.to
      out <- y / center
      out[out == 0] <- 1 # avoid risk of divide-by-zero for libraries with column sums of zero
      out
    }
    sf2 <- center_sf(sums)
    normalized <- log1p(t(t(aggregated) / sf2))
    row.vars <- rowVars(normalized)
    keep <- head(order(row.vars, decreasing=TRUE), num.features)
    row.subset <- sort(keep)

    # TF-IDF normalization (log(TF-IDF) method) and compute LSI
    lsi.res <- .computeLSI(x, ptr, col.subset = col.subset, row.subset = row.subset, lsi.method = lsi.method, scale.to = scale.to, num.dimensions = rank, outlier.quantiles = outlier.quantiles, seed = seed, binarize = binarize, num.threads = num.threads)
    embedding <- lsi.res$matSVD
  }
  beachmat::flushMemoryCache()
  
  list(embedding = embedding, projection = lsi.res$svd$v, subset = row.subset, details = lsi.res)
}

#' @importFrom S4Vectors SimpleList
#' @importFrom stats quantile
.computeLSI <- function(
    x,
    ptr,
    col.subset = seq(1, ncol(x)),
    row.subset = NULL,
    lsi.method = 1,
    scale.to = 10^4,
    num.dimensions = 50,
    outlier.quantiles = c(0.02, 0.98),
    seed = 1,
    verbose = FALSE,
    binarize = TRUE,
    num.threads = 1
){
  
  set.seed(seed)
  if(!is.null(row.subset)) {
    ptr <- apply_subset(ptr, row.subset, row = TRUE)
  }
  # compute column sums and remove columns with only zero values
  if(binarize) {
    ptr.t <- apply_transpose(ptr)
    stats <- lsi_matrix_stats(ptr.t, nthreads = num.threads) # because of transpose: sums = rowSums, frequency = # non-zero per columns
    col.sums <- stats$frequency
  } else {
    stats <- lsi_matrix_stats(ptr, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row
    col.sums <- stats$sums
  }
  if(any(col.sums == 0)){
    warning(paste0('removing ',length(which(col.sums == 0)),' columns that do not have non-zero values for the features being used.'))
    idx <- which(col.sums != 0)
    col.subset <- col.subset[idx]
    ptr <- apply_subset(ptr, idx, row = FALSE)
    col.sums <- col.sums[idx]
  }

  # filter outliers
  cn <- colnames(x)[col.subset]
  filter.outliers <- FALSE
  if(!is.null(outlier.quantiles)){
    qCS <- quantile(col.sums, probs = c(min(outlier.quantiles), max(outlier.quantiles)))
    idx.outlier <- which(col.sums <= qCS[1] | col.sums >= qCS[2])
    idx.keep <- setdiff(seq_along(col.sums), idx.outlier)
    if(length(idx.outlier) > 0){
      col.idx.outlier <- col.subset[idx.outlier]
      ptr <- apply_subset(ptr, idx.keep, row = FALSE)
      col.idx2 <- col.subset[idx.keep][head(seq_len(length(idx.keep)),100)]
      col.subset <- col.subset[idx.keep]
      col.sums <- col.sums[idx.keep]
      filter.outliers <- TRUE
    }
  }
  # clean up zero rows
  if(binarize) {
    stats <- lsi_matrix_stats(ptr, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row
    row.sums <- stats$frequency
  } else {
    ptr.t <- apply_transpose(ptr)
    stats <- lsi_matrix_stats(ptr.t, nthreads = num.threads) # because of transpose: sums = rowSums, frequency = # non-zero per columns
    row.sums <- stats$sums
  }
  idx <- which(row.sums > 0)
  row.subset <- row.subset[idx]
  row.sums <- row.sums[idx]
  # apply normalization
  if(is.null(row.subset)) {
    mat <- .apply.tf.idf.normalization(x[, col.subset], ncol(x), row.sums, lsi.method = lsi.method, scale.to = scale.to, binarize = binarize)
  } else {
    mat <- .apply.tf.idf.normalization(x[row.subset, col.subset], ncol(x), row.sums, lsi.method = lsi.method, scale.to = scale.to, binarize = binarize)
  }
  # calculate SVD then LSI
  svd <- irlba::irlba(mat, num.dimensions, num.dimensions)
  svdDiag <- matrix(0, nrow = num.dimensions, ncol = num.dimensions)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
  # Return Object
  out <- SimpleList(
    matSVD = matSVD, 
    row.sums = row.sums, 
    ncol = length(col.sums),
    svd = svd,
    scale.to = scale.to,
    num.dimensions = num.dimensions,
    lsi.method = lsi.method,
    outliers = NA,
    date = Sys.Date(),
    seed = seed, 
    row.subset = row.subset
  )
  
  if(filter.outliers){
    
    # Quick Check LSI-Projection Works
    pCheck <- .projectLSI(mat = x[row.subset, col.idx2, drop = FALSE], lsi.res = out, verbose = verbose)
    pCheck2 <- out$matSVD[rownames(pCheck), ]
    pCheck3 <- unlist(lapply(seq_len(ncol(pCheck)), function(x){
      cor(pCheck[,x], pCheck2[,x])
    }))
    if(min(pCheck3) < 0.95){
      warning("Warning with LSI-projection! Cor less than 0.95 of re-projection.")
    }
    # Project LSI Outliers
    out$outliers <- colnames(x)[col.idx.outlier]
    outlierLSI <- .projectLSI(mat = x[row.subset, col.idx.outlier, drop = FALSE], lsi.res = out, verbose = verbose)
    allLSI <- rbind(out$matSVD, outlierLSI)
  #  allLSI <- allLSI[cn, , drop = FALSE] #Re-Order to original
    out$matSVD <- allLSI
  }
  gc()
  
  out
}

#' @importFrom Matrix diag
.projectLSI <- function(mat, lsi.res, return.model = FALSE, verbose = FALSE, binarize = TRUE){   
  
  require(Matrix)
  set.seed(lsi.res$seed)
  
  # sparse matrix in memory is returned from .apply.tf.idf.normalization
  mat <- .apply.tf.idf.normalization(mat, lsi.res$ncol, lsi.res$row.sums, scale.to = lsi.res$scale.to, lsi.method = lsi.res$lsi.method, binarize = binarize) 
  
  # Clean Up Matrix
  idxNA <- Matrix::which(is.na(mat), arr.ind = TRUE)
  if(length(idxNA) > 0){
    mat[idxNA] <- 0
  }
  
  # Calc V
  V <- Matrix::t(mat) %*% lsi.res$svd$u %*% Matrix::diag(1/lsi.res$svd$d)
  
  # LSI Diagonal
  svdDiag <- matrix(0, nrow = lsi.res$num.dimensions, ncol = lsi.res$num.dimensions)
  diag(svdDiag) <- lsi.res$svd$d
  matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
  matSVD <- as.matrix(matSVD)
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
  
  if(return.model){
    X <- lsi.res$svd$u %*% diag(lsi.res$svd$d) %*% t(V)
    return(list(matSVD = matSVD, V = V, X = X))
  }else{
    return(matSVD)
  }
  
}

#' @importFrom Matrix Diagonal
# num.col and row.sums might not match mat. For example, when normalizing for the projection, num.col and run.sums will match the LSI result instead of mat.
.apply.tf.idf.normalization <- function(mat, num.col, row.sums, scale.to = 10^4, lsi.method = 1, binarize = TRUE) {
  
  if(!is(mat,'sparseMatrix')) {
    mat <- as(mat, 'sparseMatrix')
  }
  if(binarize) {
    mat@x[mat@x > 0] <- 1
  }
  # compute Term Frequency
  col.sums <- Matrix::colSums(mat)
  if(any(col.sums == 0)){
    exclude.zero.columns <- which(col.sums == 0)
    mat <- mat[,-exclude.zero.columns]
    col.sums <- col.sums[-exclude.zero.columns]
  }
  # divide each value by its column sum
  mat@x <- mat@x / rep.int(col.sums, Matrix::diff(mat@p)) 
  
  if(lsi.method == 1 | tolower(lsi.method) == "tf-logidf"){
    # Adapted from Casanovich et al.
    # compute Inverse Document Frequency
    # LogIDF
    idf <- as(log(1 + num.col / row.sums), "sparseVector")
  }else if(lsi.method == 2 | tolower(lsi.method) == "log(tf-idf)"){
    # Adapted from Stuart et al.
    # compute Inverse Document Frequency
    # IDF
    idf   <- as(num.col / row.sums, "sparseVector")
  }else if(lsi.method == 3 | tolower(lsi.method) == "logtf-logidf"){
    # LogTF
    mat@x <- log(mat@x + 1)
    # compute Inverse Document Frequency
    # LogIDF
    idf   <- as(log(1 + num.col / row.sums), "sparseVector")
  }else{
    stop("lsi.method unrecognized please select valid method!")
  }
  
  # compute TF-IDF Matrix
  # TF-IDF
  mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% mat
  
  if(lsi.method == 2 | tolower(lsi.method) == "log(tf-idf)"){
    # Log transform TF-IDF
    mat@x <- log(mat@x * scale.to + 1)  
  }
  
  gc()
  
  return(mat)
}