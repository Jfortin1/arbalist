#' Create per sample TSS enrichment plot
#'
#' Create a per sample quality control plot looking at TSS Enrichment compared to number unique fragments.
#'
#' @param mae \linkS4class{MultiAssayExperiment} containing a "TileMatrix500" experiment with "TSSEnrichment" and "fragments" columns in the colData.
#' @param sample.name String specifying the sample name to plot.
#' @param tss.threshold Numeric scalar specifying the TSS enrichment score threshold or horizontal line to mark on the plot.
#' @param nfrag.threshold Integer scalar specifying the number of fragments threshold or vertical line to mark on the plot.
#' @param max.nfrag Integer scalar specifying the maximum number of fragments to plot. Cells with a greater number of fragments will be plotted at this value.
#'
#' @author Natalie Fox
#' @importFrom ggplot2 ggplot aes theme_classic xlab ylab ggtitle geom_hline geom_vline
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom viridis scale_color_viridis
#' @importFrom SummarizedExperiment colData
#' @export
plotTSSenrichmentVsNumFragments <- function(mae, sample.name, tss.threshold = NULL, nfrag.threshold = NULL, max.nfrag = NULL) {
  tile.matrix.coldata <- colData(findSCE(mae,'TileMatrix500')$sce)
  plot.data <- as.data.frame(tile.matrix.coldata[tile.matrix.coldata$Sample == sample.name,])
  plot.data$fragments <- log10(plot.data$fragments)

  if(!is.null(max.nfrag)) {
    plot.data$fragments <- pmin(plot.data$fragments, max.nfrag)
  }
  
  plot <- ggplot(plot.data,aes(y=TSSEnrichment,x=fragments)) + 
    geom_pointdensity() + theme_classic() + ylab('TSS Enrichment') +
    xlab(bquote(log[10] ~ "Number of Fragments")) + scale_color_viridis() + ggtitle(sample.name)
  if(!is.null(tss.threshold)) {
    plot <- plot + geom_hline(yintercept=tss.threshold, linetype="dashed")
  }
  if(!is.null(nfrag.threshold)) {
    plot <- plot + geom_vline(xintercept=log10(nfrag.threshold), linetype="dashed")
  }
  plot
}

#' Plot Fragment size distribution
#'
#' Plot the fragment size distribution for a samples.
#'
#' @param mae \linkS4class{MultiAssayExperiment} containing fragment_file col in the colData.
#' @param sample.name String specifying the sample name to plot.
#' @param xlim Numeric vector with a min and max for the x-axis limits.
#'
#' @author Natalie Fox
#' @importFrom ggplot2 ggplot aes theme_classic xlab ylab geom_line ggtitle xlim
#' @importFrom SummarizedExperiment colData
#' @export
plotFragmentSizeDistribution <- function(mae, sample.name, xlim=c(0,700)) {
  fragment.size.counts <- count_fragment_size_distributions(colData(mae)$fragment_file[which(rownames(colData(mae)) == sample.name)])
  nfrags <- sum(fragment.size.counts)
  plot.data <- data.frame(size=seq_along(fragment.size.counts),frag_percentage=fragment.size.counts/nfrags*100)[seq(1, max(xlim)),]
  ggplot(plot.data,aes(y=frag_percentage,x=size)) + 
    geom_line() + theme_classic() + xlim(xlim) + ylab('Fragments (%)') +
    xlab('Size of Fragments (bp)') + ggtitle(paste0(sample.name,"\nnFrags = ",round(nfrags/1000000, digits = 2), ' M'))
}

#' Plot Cluster Confusion Matrix
#'
#' Plot the association between slusters and samples using a confusion matrix.
#'
#' @param mae \linkS4class{MultiAssayExperiment} containing a "TileMatrix500" experiment with "Clusters" and "Sample" columns in the colData.
#' 
#' @author Natalie Fox
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment colData
#' @export
plotClusterConfusionMap <- function(mae) {
  tile.matrix.coldata <- colData(findSCE(mae,'TileMatrix500')$sce)
  cM <- table(as.data.frame(tile.matrix.coldata[,c('Clusters','Sample')]))
  cM <- cM / Matrix::rowSums(cM)
  p.confusion.map <- pheatmap(
    mat = as.matrix(cM),
    border_color = 'black'
  )
}
