#' @export
#' @importFrom matrixStats rowVars
#' @importFrom beachmat initializeCpp
#' @importFrom stats kmeans
#' @importFrom BiocParallel bpparam
iterativeLSI <- function(
    x, 
    rank = 30,
    iterations = 2,
    num.features = 25000,
    cluster.method = "kmeans",
    cluster.k = 20,
    correlation.cutoff = 0.75,
    subsample = 100000,
    scale.to = 10000,
    num.threads = 1,
    seed = 5,
    total.features = 500000,
    filter.quantile = 0.995,
    BPPARAM = bpparam()
) {

#  # filter x based on row sums for top features or most variable features
#  idx.list <- list()
#  num.line.per.parallelization <- 10000
#  for(i in 1:ceiling(nrow(x)/num.line.per.parallelization)) {
#    idx.list[[i]] <- seq((i-1)*num.line.per.parallelization+1,i*num.line.per.parallelization)
#  }
#  idx.list[[length(idx.list)]] <- seq(min(idx.list[[length(idx.list)]]),nrow(x))
#  rowSumParallelized <- function(idxs, x, num.threads) {
#    x <- x[idxs,,drop=FALSE]
#    if(!is(x,'sparseMatrix')) {
#      x <- as(x, 'sparseMatrix')
#    }
#    ptr <- beachmat::initializeCpp(x)
#    stats <- lsi_matrix_stats(ptr, nthreads = num.threads)
#    return(stats$frequency)
#  }
#  res.list <- bptry(bplapply(idx.list, rowSumParallelized, x=x, num.threads=num.threads))
#  row.accessibility <- unlist(res.list)

  
  if (cluster.method == "kmeans") {
    cluster.method <- function(x) kmeans(x, centers=cluster.k)$cluster
  }
  
  # error that type 'S4' is non subsettable when trying to fetch data using tatami from a DelayedArray
  # so converting to sparseMatrix as a temporary work around until issue with seed handling in tatami_r fixed
  if(!is(x,'sparseMatrix')) {
    x <- as(x, 'sparseMatrix')
  }
  
  # select the top features based on accessibility of the binarized features 
  ptr <- beachmat::initializeCpp(x)
  stats <- lsi_matrix_stats(ptr, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row
  row.accessibility <- stats$frequency
  rm.top <- floor((1-filter.quantile) * total.features)
  top.idx <- head(order(row.accessibility, decreasing=TRUE), num.features + rm.top )[-seq_len(rm.top)]
  
  # subset to the most accessible features
  ptr.top <- apply_subset(ptr, sort(top.idx), row=TRUE)
  
  # TF-IDF normalization (log(TF-IDF) method)
  
  stats <- lsi_matrix_stats(ptr.top, nthreads = num.threads) # sums = colSums, frequency = # non-zero per row
  
  idf <- log((ncol(x) +1) / (stats$frequency+1)) # idf stands for inverse document frequency

  center_sf <- function(y) {
    center <- if (is.null(scale.to)) mean(y) else scale.to
    out <- y / center
    out[out == 0] <- 1 # avoid risk of divide-by-zero for libraries with column sums of zero.
    out
  }
  sums <- stats$sums
  sf <- center_sf(sums)

  # divide matrix by TF (term frequency) 
  transformed <- apply_division(ptr.top, sf, right=TRUE, row=FALSE)
  # then multiple by IDF
  transformed <- apply_multiplication(ptr.top, idf, row=TRUE)
  # Finally log transform
  transformed <- apply_log1p(transformed)
  
  # Single value decomposition to get embedding
  current <- irlba_realized(transformed, rank = rank, nthreads = num.threads, seed = seed)
  embedding <- t(t(current$u) * current$d)
  
  # filter to the total number of features specified by total.features
  ptr2 <- apply_subset(ptr, sort(head(order(row.accessibility, decreasing = TRUE), total.features)), row=TRUE)
  
  for (i in seq_len(iterations-1)) {
    # Drop embeddings that are correlated with the library size.
    # Hey, I don't make the rules, I just do what I'm told.
    embedding <- embedding[,which(cor(embedding, sums) <= correlation.cutoff),drop=FALSE]
    
    # Finding clusters.
    cluster.output <- cluster.method(embedding)
    
    # Aggregating counts and finding the most variable features.
    aggregated <- aggregate_counts(ptr2, as.integer(factor(cluster.output)) - 1L, nthreads = num.threads)
    sums <- colSums(aggregated)
    sf2 <- center_sf(sums)
    normalized <- log1p(t(t(aggregated) / sf2))
    row.vars <- rowVars(normalized)
    keep <- head(order(row.vars, decreasing=TRUE), num.features)
    
    indices <- sort(keep) 
    indicesm1 <- indices - 1L
    subset <- apply_subset(ptr2, indicesm1, row = TRUE)
    
    # Re-normalizing on the subset.
    sums <- lsi_matrix_stats(subset, nthreads = num.threads)$sums
    sf <- center_sf(sums)
    transformed <- apply_division(subset, sf, right = TRUE, row = FALSE)
    transformed <- apply_multiplication(transformed, idf[indices], row = TRUE)
    transformed <- apply_log1p(transformed)
    
    # Computing the next iteration.
    seed <- sample(.Machine$integer.max, 1)
    current <- irlba_realized(subset, rank = rank, nthreads = num.threads, seed = seed + i)
    embedding <- t(t(current$u) * current$d)
  }
  
  list(embedding = embedding, projection = current$v, subset = indices)
}