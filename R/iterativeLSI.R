#' Calculate iterative LSI dimensionality reduction
#'
#' @param x input matrix for creating iterative LSI embeddings, assumes ATAC-seq so chooses top accessible features.
#' @param cell.names String vector specifying cell names. Length should match the number of columns in x.
#' @param sample.names String vector specifying sample names. Length should match the number of columns in x.
#' @param cell.depths String vector specifying cell depths. Length should match the number of columns in x.
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
#' @param num.cells.to.sample Scalar integer specifying the number of cells to sample iterations prior to the last in order to perform sub-sampled LSI and sub-sampled clustering.
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
#' @importFrom beachmat initializeCpp flushMemoryCache tatami.compare tatami.subset tatami.arith  tatami.row.sums tatami.column.sums tatami.dim
#' @importFrom stats kmeans cor
#' @importFrom utils head
iterativeLSI <- function(
    x,
    cell.names,
    sample.names,
    cell.depths,
    col.subset = NULL,
    rank = 30, 
    iterations = 2,
    first.selection = "Top", # or "Var"
    num.features = 25000,
    lsi.method = 2,
    cluster.method = "Seurat", 
    correlation.cutoff = 0.75,
    scale.to = 10000,
    num.threads = 4,
    seed = 1,
    total.features = 500000,
    filter.quantile = 0.995,
    outlier.quantiles = c(0.02, 0.98),
    binarize = TRUE,
    num.cells.to.sample = 10000
) {
  if(num.features < 1000) {
    stop('Please specify a number of features equal to or more than 1000 to num.features.')
  }
  
  col.names <- colnames(x)
  
  # select the top features based on accessibility of the binarized features (ex for ATAC data) or variance (ex for RNA data)
  beachmat::flushMemoryCache()
  if(is(x,"DelayedArray")) {
    ptr <- beachmat.hdf5::initializeCpp(x, memorize=TRUE)
  } else {
    ptr <- beachmat::initializeCpp(x, memorize = TRUE)
  }
  if(!is.null(col.subset)) {
    col.subset <- sort(unique(col.subset),decreasing = TRUE)
    if(length(col.subset) <= ncol(x)) {
     ptr <- tatami.subset(ptr, subset = col.subset, by.row = FALSE)
     col.names <- col.names[col.subset]
    } else {
      stop('There are more entries in the column subset to select than the number of columns.')
    }
  }
  if(tolower(first.selection) == "top") {
      row.order.stat <- tatami.row.sums(ptr, num.threads = num.threads)
  } else if (tolower(first.selection) == "var") {
    if(binarize) {
      stop("binarize must be FALSE if first.selection is Var")
    }
    row.means <- tatami.row.sums(ptr, num.threads = num.threads)/ncol(x)
    ptr2 <- tatami.arith(ptr, op = '-', val = row.means, by.row = TRUE, right = TRUE)
    ptr2 <- tatami.arith(ptr2, op = '^', val = 2, by.row = TRUE, right = TRUE)
    row.order.stat <- tatami.row.sums(ptr2, num.threads = num.threads)/(ncol(x)-1) # row variance
  } else {
    stop('first.selection not "Top" or "Var".')
  }
  
  ptr.orig.not.binarized <- ptr
  if(binarize) {
    ptr <- tatami.compare(ptr, op = '!=', val = 0, by.row = TRUE, right = TRUE)
  }
  
  col.sums <- tatami.column.sums(ptr, num.threads = num.threads)
  
  rm.top <- floor((1-filter.quantile) * total.features)
  row.subset <- sort(head(order(row.order.stat, decreasing = TRUE), num.features + rm.top )[-seq_len(rm.top)])
  ptr.subset <- tatami.subset(ptr, subset = row.subset, by.row = TRUE)
  
  # TF-IDF normalization (log(TF-IDF) method) and compute LSI
  lsi.res <- .computeLSI(
    ptr.subset, 
    cell.names = cell.names, 
    sample.names = sample.names, 
    cell.depths = cell.depths, 
    lsi.method = lsi.method, 
    scale.to = scale.to, 
    num.dimensions = rank, 
    outlier.quantiles = outlier.quantiles, 
    seed = seed, 
    num.threads = num.threads, 
    num.cells.to.sample = num.cells.to.sample, 
    project.all.cells = FALSE
    )
  embedding <- lsi.res$matSVD

  for (i in seq_len(iterations-1)) {
    # drop embeddings that are correlated with the library size
    embedding <- embedding[,which(cor(embedding, col.sums[cell.names %in% rownames(embedding)]) <= correlation.cutoff),drop=FALSE]
    
    # find cell clusters
    bias.vals <- log10(cell.depths + 1)
    names(bias.vals) <- cell.names
    cluster.output <- .clusterMatrix(embedding, method = cluster.method, bias.vals = bias.vals[rownames(embedding)])
    if(length(table(cluster.output)) == 1) {
      warning('Data is not splitting into clusters so we cannot calculate iterativeLSI')
      return(NULL)
    }

    # aggregate counts for each cluster and find the most variable features
    group.features.idx <- sort(head(order(row.order.stat, decreasing = TRUE), total.features))
    ptr.subset <- tatami.subset(ptr.orig.not.binarized, subset = group.features.idx, by.row = TRUE)
    ptr.subset <- tatami.subset(ptr.subset, subset = which(cell.names %in% rownames(embedding)), by.row = FALSE)
    aggregated <- .aggregateGroupCounts(ptr.subset, cluster.output[cell.names[which(cell.names %in% names(cluster.output))]], num.threads = num.threads)
    
    normalized <- log2(t(t(aggregated) / colSums(aggregated)) * scale.to + 1)
    row.vars <- rowVars(normalized)
    keep <- group.features.idx[head(order(row.vars, decreasing=TRUE), num.features)]
    row.subset <- sort(keep)
    ptr.subset <- tatami.subset(ptr, subset = row.subset, by.row = TRUE)

    # TF-IDF normalization (log(TF-IDF) method) and compute LSI
    lsi.res <- .computeLSI(
      ptr.subset, 
      cell.names = cell.names, 
      sample.names = sample.names, 
      cell.depths = cell.depths,
      lsi.method = lsi.method, 
      scale.to = scale.to, 
      num.dimensions = rank, 
      outlier.quantiles = outlier.quantiles, 
      seed = seed, 
      num.threads = num.threads, 
      num.cells.to.sample = if(i != (iterations-1)) num.cells.to.sample else NULL, 
      project.all.cells = TRUE
      )

    embedding <- lsi.res$matSVD
  }
  beachmat::flushMemoryCache()
  
  list(embedding = embedding, projection = lsi.res$svd$v, subset = row.subset, details = lsi.res)
}

#' @importFrom S4Vectors SimpleList
#' @importFrom stats quantile
#' @importFrom beachmat tatami.column.sums tatami.subset tatami.row.sums tatami.realize tatami.dim
.computeLSI <- function(
    ptr,
    cell.names,
    sample.names,
    cell.depths,
    lsi.method = 1,
    scale.to = 10^4,
    num.dimensions = 50,
    outlier.quantiles = c(0.02, 0.98),
    seed = 1,
    verbose = FALSE,
    num.threads = 4,
    num.cells.to.sample = 10000,
    project.all.cells = TRUE
){
  set.seed(seed)

  cell.depths <- log10(cell.depths + 1)
  if(!is.null(num.cells.to.sample)) {
    sampleN <- ceiling(num.cells.to.sample * table(sample.names) / length(sample.names))
    split.cells <- split(cell.names, sample.names)
    split.depth <- split(cell.depths, sample.names)
    
    sampled.cells <- cell.names[cell.names %in% unlist(sapply(1:length(split.cells), function(i) {
      cell.list <- split.cells[[i]]
      n <- sampleN[names(split.cells)[i]]
      
      if(!is.null(outlier.quantiles)){
        quant <- quantile(split.depth[[i]], probs = c(min(outlier.quantiles) / 2, 1 - ((1-max(outlier.quantiles)) / 2)))
        idx <- which(split.depth[[i]] >= quant[1] & split.depth[[i]] <= quant[2])
      }else{
        idx <- seq_along(cell.list)
      }
      if(length(idx) >= n){
        sample(x = cell.list[idx], size = n)
      } else if (length(cell.list) >= n){
        sample(x = cell.list, size = n)
      } else {
        cell.list
      }
    }))]
  } else {
    sampled.cells <- cell.names
  }
  
  # compute column sums and remove columns with only zero values
  col.sums <- tatami.column.sums(ptr, num.threads = num.threads)
  if(any(col.sums == 0)){
    warning(paste0('removing ',length(which(col.sums == 0)),' columns that do not have non-zero values for the features being used.'))
    idx <- which(col.sums != 0)
    ptr <- tatami.subset(ptr, subset = idx, by.row = FALSE)
    col.sums <- col.sums[idx]
  }

  # filter outliers
  idx.outlier <- NULL
  idx.keep <- seq_along(col.sums)
  idx.not.in.sampled <- NULL
  if(any(!cell.names %in% sampled.cells)) {
    idx.not.in.sampled <- which(!cell.names %in% sampled.cells)
    idx.keep <- setdiff(idx.keep, idx.not.in.sampled)
  }
  if(!is.null(outlier.quantiles)){
    qCS <- quantile(col.sums[idx.keep], probs = c(min(outlier.quantiles), max(outlier.quantiles)))
    idx.outlier <- idx.keep[which(col.sums[idx.keep] <= qCS[1] | col.sums[idx.keep] >= qCS[2])]
  }
  if(project.all.cells) {
    idx.outlier <- sort(union(idx.outlier, idx.not.in.sampled))
  }
  if(length(idx.outlier) > 0) {
    idx.keep <- setdiff(idx.keep, idx.outlier)
    ptr.outliers <- tatami.subset(ptr, subset = idx.outlier, by.row = FALSE)
    ptr <- tatami.subset(ptr, subset = idx.keep, by.row = FALSE)
    idx.keep.subset <- head(seq_len(length(idx.keep)),50)
    ptr.keep.subset <- tatami.subset(ptr, subset = idx.keep.subset, by.row = FALSE)
    col.sums <- col.sums[idx.keep]
  } else if(length(idx.keep) < length(col.sums)){
    ptr <- tatami.subset(ptr, subset = idx.keep, by.row = FALSE)
  }

  # clean up zero rows
  row.sums <- tatami.row.sums(ptr, num.threads = num.threads)
  row.subset <- which(row.sums > 0)
  if(any(row.subset)) {
    ptr <- tatami.subset(ptr, subset = row.subset, by.row = TRUE)
    row.sums <- row.sums[row.subset]
  }
  
  # apply normalization
  mat <- .applyTFIDFNormalization(ptr, length(col.sums), row.sums, lsi.method = lsi.method, scale.to = scale.to, num.threads = num.threads)
  colnames(mat) <- cell.names[idx.keep]
  
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
  
  if(length(idx.outlier) > 0){
    
    # Quick Check LSI-Projection Works
    ptr.keep.subset <- tatami.subset(ptr.keep.subset, subset = row.subset, by.row = TRUE)
    pCheck <- .projectLSI(ptr.keep.subset, lsi.res = out, verbose = verbose)
    #pCheck2 <- out$matSVD[rownames(pCheck), ]
    pCheck2 <- out$matSVD[idx.keep.subset, ]
    pCheck3 <- unlist(lapply(seq_len(ncol(pCheck)), function(x){
      cor(pCheck[,x], pCheck2[,x])
    }))
    if(min(pCheck3) < 0.95){
      warning("Warning with LSI-projection! Cor less than 0.95 of re-projection.")
    }
    # Project LSI Outliers
    out$outliers <- idx.outlier
    ptr.outliers <- tatami.subset(ptr.outliers, subset = row.subset, by.row = TRUE)
    outlierLSI <- .projectLSI(ptr.outliers, lsi.res = out, verbose = verbose)
    allLSI <- rbind(out$matSVD, outlierLSI)
    rownames(allLSI) <- cell.names[c(idx.keep,idx.outlier)]
    if(length(allLSI) == length(sampled.cells)) {
      allLSI <- allLSI[sampled.cells,]
    } else if(length(allLSI) == length(cell.names)) {
      allLSI <- allLSI[cell.names,]
    }
    out$matSVD <- allLSI
  }
  gc()
  
  out
}

#' @importFrom Matrix diag
.projectLSI <- function(ptr, lsi.res, return.model = FALSE, verbose = FALSE, num.threads = 4){   
  
  require(Matrix)
  set.seed(lsi.res$seed)
  
  # sparse matrix in memory is returned from .applyTFIDFNormalization
  mat <- .applyTFIDFNormalization(ptr, lsi.res$ncol, lsi.res$row.sums, scale.to = lsi.res$scale.to, lsi.method = lsi.res$lsi.method, num.threads = num.threads) 
  
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
#' @importFrom beachmat tatami.column.sums tatami.subset tatami.arith tatami.realize
# num.col and row.sums might not match mat. For example, when normalizing for the projection, num.col and run.sums will match the LSI result instead of mat.
.applyTFIDFNormalization <- function(ptr, num.col, row.sums, scale.to = 10^4, lsi.method = 2, num.threads = 4) {
  
  # compute Term Frequency
  col.sums <- tatami.column.sums(ptr, num.threads = num.threads)
  if(any(col.sums == 0)){
    non.zero.columns <- which(col.sums != 0)
    ptr <- tatami.subset(ptr, subset = non.zero.columns, by.row = FALSE)
    col.sums <- col.sums[non.zero.columns]
  }
  # divide each value by its column sum
  ptr <- tatami.arith(ptr, op = '/', col.sums, by.row = FALSE, right = TRUE)

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
  mat <- beachmat::tatami.realize(ptr, num.threads = num.threads)
  mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% mat
  
  if(lsi.method == 2 | tolower(lsi.method) == "log(tf-idf)"){
    # Log transform TF-IDF
    mat@x <- log(mat@x * scale.to + 1)  
  }
  
  gc()
  
  return(mat)
}

.aggregateGroupCounts <- function(ptr, col.groups, num.threads = 1) {
  group.labels <- sort(unique(col.groups))
  ptr.subset <- tatami.subset(ptr, subset = which(col.groups == group.labels[1]), by.row = FALSE)
  row.sums <- tatami.row.sums(ptr.subset, num.threads = num.threads)
  mat.res <- matrix(NA,ncol=length(group.labels), nrow=length(row.sums))
  colnames(group.labels)
  mat.res[,1] <- row.sums
  for(i in (seq_along(group.labels))[-1]) {
    ptr.subset <- tatami.subset(ptr, subset = which(col.groups == group.labels[i]), by.row = FALSE)
    mat.res[,i] <- tatami.row.sums(ptr.subset, num.threads = num.threads)
  }
  return(mat.res)
}