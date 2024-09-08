#' Calculate iterative LSI dimensionality reduction
#'
#' @param x input matrix for creating iterative LSI embeddings, assumes ATAC-seq so chooses top accessible features.
#' @param col.subset Integer vector specifying the subset of the columns to include in the iterative LSI calculations.
#' @param rank Integer scalar specifying the rank for irlba_realized.
#' @param iterations Integer scalar specifying number of LSI iterations to perform.
#' @param first.selection String specifying either "Top" for the most accessible features or "Var" for the most variable features.
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
#' @importFrom stats kmeans cor
#' @importFrom utils head
#' @importFrom SparseArray rowSums rowVars colSums
iterativeLSI.sparse.in.mem <- function(
    x,
    rank = 30,
    iterations = 2,
    first.selection = "Var", # or "Top"
    num.features = 25000,
    lsi.method = 1,
    cluster.method = "Seurat",
    correlation.cutoff = 0.75,
    scale.to = 10000,
    num.threads = 4,
    seed = 5,
    total.features = 500000,
    filter.quantile = 0.995,
    outlier.quantiles = c(0.02, 0.98),
    binarize = FALSE
) {
  col.names <- colnames(x)
  
  # select the top features based on accessibility of the binarized features (ex for ATAC data) or variance (ex for RNA data)
  if(binarize) {
    x@x[x@x > 0] <- 1
  }
  if(first.selection == "Top") {
    x2 <- x
    x2@x[x2@x > 0] <- 1
    row.order.stat <- SparseArray::rowSums(x2)
  } else if (first.selection == "Var") {
    #row.order.stat <- apply(x,1,var)
    row.order.stat <- SparseArray::rowVars(x)
  } else {
    stop('first.selection not "Top" or "Var".')
  }
  
  col.sums <- SparseArray::colSums(x)
  
  rm.top <- floor((1-filter.quantile) * total.features)
  row.subset <- head(order(row.order.stat, decreasing = TRUE), num.features + rm.top )[-seq_len(rm.top)]

  # TF-IDF normalization (log(TF-IDF) method) and compute LSI
  lsi.res <- .computeLSI.in.mem(x[row.subset,], lsi.method = lsi.method, scale.to = scale.to, num.dimensions = rank, outlier.quantiles = outlier.quantiles, seed = seed, num.threads = num.threads)
  embedding <- lsi.res$matSVD
    
  for (i in seq_len(iterations-1)) {
    # drop embeddings that are correlated with the library size
    embedding <- embedding[,which(cor(embedding, col.sums[rownames(embedding)]) <= correlation.cutoff),drop=FALSE]
    
    # find cell clusters
    cluster.output <- cluster.matrix(embedding, method = cluster.method)
    if(length(table(cluster.output)) == 1) {
      warning('Data is not splitting into clusters so we cannot calculate iterativeLSI')
      return(NULL)
    }
    # aggregate counts for each cluster and find the most variable features
    #aggregated <- bplapply(sort(unique(cluster.output)), function(x, mat, col.groups){
    #  return(rowSums(x[,col.groups = x]))
    #}, x, cluster.output, BPPARAM = bpparam())
    aggregated <- sapply(sort(unique(cluster.output)), function(clust.name, mat, col.groups){
      return(SparseArray::rowSums(mat[,col.groups == clust.name]))
    }, x, cluster.output)
    
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
    lsi.res <- .computeLSI.in.mem(x[row.subset,], lsi.method = lsi.method, scale.to = scale.to, num.dimensions = rank, outlier.quantiles = outlier.quantiles, seed = seed, num.threads = num.threads)
    embedding <- lsi.res$matSVD
  }
  beachmat::flushMemoryCache()
  
  list(embedding = embedding, projection = lsi.res$svd$v, subset = row.subset, details = lsi.res)
}

#' @importFrom S4Vectors SimpleList
#' @importFrom stats quantile
#' @importFrom beachmat tatami.column.sums tatami.subset tatami.row.sums tatami.realize tatami.dim
.computeLSI.in.mem <- function(
    x,
    lsi.method = 1,
    scale.to = 10^4,
    num.dimensions = 50,
    outlier.quantiles = c(0.02, 0.98),
    seed = 1,
    verbose = FALSE,
    num.threads = 4
){
  set.seed(seed)
  # compute column sums and remove columns with only zero values
  col.sums <- SparseArray::colSums(x)
  if(any(col.sums == 0)){
    warning(paste0('removing ',length(which(col.sums == 0)),' columns that do not have non-zero values for the features being used.'))
    x <- x[,which(col.sums != 0)]
    col.sums <- col.sums[which(col.sums != 0)]
  }

  # filter outliers
  filter.outliers <- FALSE
  if(!is.null(outlier.quantiles)){
    qCS <- quantile(col.sums, probs = c(min(outlier.quantiles), max(outlier.quantiles)))
    idx.outlier <- which(col.sums <= qCS[1] | col.sums >= qCS[2])
    idx.keep <- setdiff(seq_along(col.sums), idx.outlier)
    if(length(idx.outlier) > 0){
      x.outliers <- x[,idx.outlier]
      x <- x[, idx.keep]
      idx.keep.subset <- head(seq_len(length(idx.keep)),100)
      x.keep.subset <- x[,idx.keep.subset]
      col.sums <- col.sums[idx.keep]
      filter.outliers <- TRUE
    }
  }
  # clean up zero rows
  row.sums <- SparseArray::rowSums(x)
  row.subset <- which(row.sums > 0)
  if(any(row.subset)) {
    x<- x[row.subset,]
    row.sums <- row.sums[row.subset]
  }
  
  # apply normalization
  mat <- .apply.tf.idf.normalization.in.mem(x, length(col.sums), row.sums, lsi.method = lsi.method, scale.to = scale.to)
  
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
    x.keep.subset <- x.keep.subset[row.subset,]
    pCheck <- .projectLSI.in.mem(x.keep.subset, lsi.res = out, verbose = verbose)
    pCheck2 <- out$matSVD[idx.keep.subset, ]
    pCheck3 <- unlist(lapply(seq_len(ncol(pCheck)), function(x){
      cor(pCheck[,x], pCheck2[,x])
    }))
    if(min(pCheck3) < 0.95){
      warning("Warning with LSI-projection! Cor less than 0.95 of re-projection.")
    }
    # Project LSI Outliers
    out$outliers <- idx.outlier
    x.outliers <- x.outliers[row.subset,]
    outlierLSI <- .projectLSI.in.mem(x.outliers, lsi.res = out, verbose = verbose)
    allLSI <- rbind(out$matSVD, outlierLSI)
    out$matSVD <- allLSI
  }
  gc()
  
  out
}

#' @importFrom Matrix diag
.projectLSI.in.mem <- function(x, lsi.res, return.model = FALSE, verbose = FALSE, num.threads = 4){   
  
  require(Matrix)
  set.seed(lsi.res$seed)
  
  # sparse matrix in memory is returned from .apply.tf.idf.normalization
  mat <- .apply.tf.idf.normalization.in.mem(x, lsi.res$ncol, lsi.res$row.sums, scale.to = lsi.res$scale.to, lsi.method = lsi.res$lsi.method) 
  
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
.apply.tf.idf.normalization.in.mem <- function(mat, num.col, row.sums, scale.to = 10^4, lsi.method = 1) {
  
  if(!is(mat,'sparseMatrix')) {
    mat <- as(mat, 'sparseMatrix')
  }
#  if(binarize) {
#    mat@x[mat@x > 0] <- 1
#  }
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