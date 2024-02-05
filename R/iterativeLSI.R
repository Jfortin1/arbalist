#' Create iterative LSI embeddings
#'
#' @param x input matrix for creating iterative LSI embeddings, assumes ATAC-seq so chooses top accessible features
#' @param rank Number specifying the rank for irlba_realized
#' @param iterations Number of LSI iterations to perform.
#' @param num.features Number of accessible features to select when selecting the most accessible or most variable features.
#' @param lsi.method Number or string indicating the order of operations in the TF-IDF normalization. Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf"
#' @param cluster.method String specifying cluster method. Currently "kmeans" is the only supported option.
#' @param cluster.k Number of clusters to use
#' @param correlation.cutoff 	Numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded from analysis.
#' @param scale.to Number specifying the center for TF-IDF normalization.
#' @param num.threads Number of threads to be used for parallel computing.
#' @param seed Number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param total.features Number of features to consider for use in LSI after ranking the features by the total number of insertions. These features are the only ones used throughout the variance identification and LSI.
#' @param filter.quantile Number 0,1 that indicates the quantile above which features should be removed based on insertion counts prior
#'
#' @return A list is returned containing:
#' \itemize{
#' \item \code{embedding}, a matrix containing the iterativeLSI embedding
#' \item \code{projection}, 
#' \item \code{subset}, a vectory of indices specifying the selected subset of features
#' }
#' 

#' @author Natalie Fox
#' @examples
#' \dontrun{
#' mae <- maw.scatac::importScAtac('FRS15024')
#' res <- iterativeLSI(MultiAssayExperiment::assay(mae[[1]]))
#'}
#' @export
#' @importFrom matrixStats rowVars
#' @importFrom beachmat initializeCpp flushMemoryCache
#' @importFrom stats kmeans cor
#' @importFrom utils head
iterativeLSI <- function(
    x, 
    rank = 30,
    iterations = 2,
    num.features = 25000,
    lsi.method = 1,
    cluster.method = "kmeans",
    cluster.k = 20,
    correlation.cutoff = 0.75,
    scale.to = 10000,
    num.threads = 1,
    seed = 5,
    total.features = 500000,
    filter.quantile = 0.995,
    outlier.quantiles = c(0.02, 0.98)
) {

  if (cluster.method == "kmeans") {
    cluster.method <- function(x) kmeans(x, centers=cluster.k)$cluster
  }
  
  # error that type 'S4' is non subsettable when trying to fetch data using tatami from a DelayedArray
  # so converting to sparseMatrix as a temporary work around until issue with seed handling in tatami_r fixed
#  if(!is(x,'sparseMatrix')) {
#    x <- as(x, 'sparseMatrix')
#  }
  
  # select the top features based on accessibility of the binarized features
  beachmat::flushMemoryCache()
  ptr <- beachmat::initializeCpp(x, memorize=TRUE)
  stats <- lsi_matrix_stats(ptr, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row
  row.accessibility <- stats$frequency
  rm.top <- floor((1-filter.quantile) * total.features)
  top.idx <- head(order(row.accessibility, decreasing=TRUE), num.features + rm.top )[-seq_len(rm.top)]
  
  # subset to the most accessible features
  mat.subset <- x[top.idx,]
  
  # TF-IDF normalization (log(TF-IDF) method) and compute LSI
  lsi.res <- .computeLSI(mat.subset, lsi.method = lsi.method, scale.to=scale.to, outlier.quantiles=outlier.quantiles, seed=seed)
  embedding <- lsi.res$matSVD
  
  for (i in seq_len(iterations-1)) {
    # drop embeddings that are correlated with the library size
    embedding <- embedding[,which(cor(embedding, stats$sums[colnames(mat.subset) %in% rownames(embedding)]) <= correlation.cutoff),drop=FALSE]
    
    # find cell clusters
    cluster.output <- cluster.method(embedding)
    
    # aggregate counts for each cluster and find the most variable features
    aggregated <- aggregate_counts(ptr, as.integer(factor(cluster.output)) - 1L, nthreads = num.threads)
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
    indices <- sort(keep) 
    
    # subset to the most varaible features
    mat.subset <- x[indices,]
    
    # TF-IDF normalization (log(TF-IDF) method) and compute LSI
    lsi.res <- .computeLSI(mat.subset, lsi.method=lsi.method, scale.to=scale.to, outlier.quantiles=outlier.quantiles, seed=seed)
    embedding <- lsi.res$matSVD
  }
  beachmat::flushMemoryCache()
  
  list(embedding = embedding, projection = lsi.res$v, subset = indices)
}

#' @importFrom S4Vectors SimpleList
.computeLSI <- function(
    mat, 
    lsi.method = 1,
    scale.to = 10^4,
    num.dimensions = 50,
    outlier.quantiles = c(0.02, 0.98),
    seed = 1, 
    verbose = FALSE
){
  
  set.seed(seed)
  
  if(!is(mat,'sparseMatrix')) {
    mat <- as(mat, 'sparseMatrix')
  }
  
  # make sure the matrix is binarized
  mat@x[mat@x > 0] <- 1

  # compute column sums and remove columns with only zero values
  col.sums <- Matrix::colSums(mat)
  if(any(col.sums == 0)){
    exclude.zero.columns <- which(col.sums==0)
    warning(paste0('removing ',length(exclude.zero.columns),' columns that do not have non-zero values for the features being used.'))
    mat <- mat[,-exclude.zero.columns, drop = FALSE]
    col.sums <- col.sums[-exclude.zero.columns]
  }
  
  # filter outliers
  cn <- colnames(mat)
  filter.outliers <- 0
  if(!is.null(outlier.quantiles)){
    qCS <- quantile(col.sums, probs = c(min(outlier.quantiles), max(outlier.quantiles)))
    idx.outlier <- which(col.sums <= qCS[1] | col.sums >= qCS[2])
    if(length(idx.outlier) > 0){
      matO <- mat[, idx.outlier, drop = FALSE]
      mat <- mat[, -idx.outlier, drop = FALSE]
      mat2 <- mat[, head(seq_len(ncol(mat)), 100), drop = FALSE] # A 2nd Matrix to Check Projection is Working
      col.sums <- col.sums[-idx.outlier]
      filter.outliers <- 1       
    }
  }
  
  # clean up zero rows
  row.sums <- Matrix::rowSums(mat)
  idx <- which(row.sums > 0)
  mat <- mat[idx, ]
  row.sums <- row.sums[idx]
  
  # apply normalization
  mat <- .apply.tf.idf.normalization(mat, ncol(mat), row.sums, lsi.method = lsi.method, scale.to = scale.to) 
  
  # calculate SVD then LSI
  svd <- irlba::irlba(as.matrix(mat), num.dimensions, num.dimensions)
  svdDiag <- matrix(0, nrow=num.dimensions, ncol=num.dimensions)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
  
  # Return Object
  out <- SimpleList(
    matSVD = matSVD, 
    row.sums = row.sums, 
    ncol = length(col.sums),
    idx = idx, 
    svd = svd,
    scale.to = scale.to,
    num.dimensions = num.dimensions,
    lsi.method = lsi.method,
    outliers = NA,
    date = Sys.Date(),
    seed = seed
  )
  
  if(filter.outliers == 1){
    # Quick Check LSI-Projection Works
    pCheck <- .projectLSI(mat = mat2, lsi.res = out, verbose = verbose)
    pCheck2 <- out[[1]][rownames(pCheck), ]
    pCheck3 <- unlist(lapply(seq_len(ncol(pCheck)), function(x){
      cor(pCheck[,x], pCheck2[,x])
    }))
    if(min(pCheck3) < 0.95){
      warning("Warning with LSI-projection! Cor less than 0.95 of re-projection.")
    }
    # Project LSI Outliers
    out$outliers <- colnames(matO)
    outlierLSI <- .projectLSI(mat = matO, lsi.res = out, verbose = verbose)
    allLSI <- rbind(out[[1]], outlierLSI)
    allLSI <- allLSI[cn, , drop = FALSE] #Re-Order Correctly to original
    out[[1]] <- allLSI
  }
  
  rm(mat)
  gc()
  
  out
}

.projectLSI <- function(mat, lsi.res, return.model = FALSE, verbose = FALSE){   
  
  require(Matrix)
  set.seed(lsi.res$seed)
  
  # Get Same Features - ie subsetting by Non-Zero features in inital Matrix
  mat <- mat[lsi.res$idx,]
  
  if(!is(mat,'sparseMatrix')) {
    mat <- as(mat, 'sparseMatrix')
  }
  
  # Binarize Matrix
  mat@x[mat@x > 0] <- 1     

  mat <- .apply.tf.idf.normalization(mat, lsi.res$ncol, lsi.res$row.sums, scale.to = lsi.res$scale.to, lsi.method = lsi.res$lsi.method) 
  
  # Clean Up Matrix
  idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
  if(length(idxNA) > 0){
    mat[idxNA] <- 0
  }
  
  # Calc V
  V <- Matrix::t(mat) %*% lsi.res$svd$u %*% Matrix::diag(1/lsi.res$svd$d)
  
  # LSI Diagonal
  svdDiag <- matrix(0, nrow=lsi.res$num.dimensions, ncol=lsi.res$num.dimensions)
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

# num.col and row.sums might not match mat. For example, when normalizing for the projection, mun.col na drun.sums will match the LSI result instead of mat.
.apply.tf.idf.normalization <- function(mat, num.col, row.sums, scale.to = 10^4, lsi.method = 1) {
  if(!is(mat,'sparseMatrix')) {
    mat <- as(mat, 'sparseMatrix')
  }
  
  # compute Term Frequency
  col.sums <- Matrix::colSums(mat)
  if(any(col.sums == 0)){
    exclude.zero.columns <- which(col.sums==0)
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
  mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
  
  if(lsi.method == 2 | tolower(lsi.method) == "log(tf-idf)"){
    # Log transform TF-IDF
    mat@x <- log(mat@x * scale.to + 1)  
  }
  
  gc()
  
  return(mat)
}

