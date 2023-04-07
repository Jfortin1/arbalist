#' @export
#' @importFrom matrixStats rowVars
iterativeLSI <- function(x, 
    rank = 30, 
    iterations = 2, 
    num.features = 25000, 
    cluster.method = "kmeans", 
    cluster.k = 20, 
    correlation.cutoff = 0.75,
    subsample = 100000,
    num.threads = 1) 
{
    if (cluster.method == "kmeans") {
        cluster.method <- function(x) kmeans(x, centers=cluster.k)$cluster
    }

    ptr <- initializeCpp(x)
    stats <- lsi_matrix_utils(ptr, nthreads = num.threads)
    idf <- log(stats$detected / ncol(X))

    center_sf <- function(x) {
        center <- if (is.null(scale.to)) mean(x) else scale.to
        out <- x / center
        out[out == 0] <- 1 # avoid risk of divide-by-zero for libraries with column sums of zero.
        out
    }
    sums <- stats$sums
    sf <- center_sf(sums)

    # Initial calculation on the full matrix. 
    transformed <- apply_division(ptr, sf, right = TRUE, row = FALSE)
    transformed <- apply_multiplication(transformed, idf, row = TRUE)
    transformed <- apply_log1p(transformed)
    current <- irlba_tatami(transformed, rank = rank, nthreads = num.threads, seed = seed)
    embedding <- t(t(current$u) * current$d)

    for (i in seq_len(iterations)) {
        # Drop embeddings that are correlated with the library size.
        # Hey, I don't make the rules, I just do what I'm told.
        embedding <- embedding[,which(cor(embedding, sums) <= correlation.cutoff),drop=FALSE]

        # Finding clusters.
        output <- cluster.method(embedding)

        # Aggregating counts and finding the most variable features.
        aggregated <- aggregate_counts(x, as.integer(factor(output)) - 1L, nthreads = num.threads)
        sums <- colSums(aggregated)
        sf <- center_sf(sums)
        normalized <- log1p(t(t(aggregated) / sf2))
        out <- rowVars(aggregated)
        keep <- head(order(aggregated, decreasing=TRUE), num.features)

        indices <- sort(keep) 
        indicesm1 <- indices - 1L
        subset <- apply_subset(x, indicesm1, row = TRUE)

        # Re-normalizing on the subset.
        sf <- center_sf(tatami_colsums(subset, nthreads = num.threads))
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
