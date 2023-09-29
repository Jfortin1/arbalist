#' Create iterative LSI embeddings
#'
#' @param x input matrix for creating iterative LSI embeddings
#' @param rank Number specifying the rank for irlba_realized
#' @param iterations Number of LSI iterations to perform.
#' @param num.features Number of accessible features to select when selecting the most accessible or most variable features.
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
#' @importFrom beachmat initializeCpp
#' @importFrom stats kmeans cor
#' @importFrom utils head
iterativeLSI <- function(
    x, 
    rank = 30,
    iterations = 2,
    num.features = 25000,
    cluster.method = "kmeans",
    cluster.k = 20,
    correlation.cutoff = 0.75,
    scale.to = 10000,
    num.threads = 1,
    seed = 5,
    total.features = 500000,
    filter.quantile = 0.995
) {

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