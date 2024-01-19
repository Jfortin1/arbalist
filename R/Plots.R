#' Create per sample TSS enrichment plot
#'
#' Create a per sample quality control plot looking at TSS Enrichment compared to number unique fragments.
#'
#' @author Natalie Fox
#' @importFrom ggplot2 ggplot aes theme_classic xlab
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom viridis scale_color_viridis
#' @importFrom SummarizedExperiment colData
#' @export
plotTSSenrichmentVsNumFragments <- function(mae, sample.name) {
  tile.matrix.coldata <- colData(findSCE(mae,'TileMatrix500')$sce)
  plot.data <- as.data.frame(tile.matrix.coldata[tile.matrix.coldata$Sample == sample.name,])
  ggplot(plot.data,aes(y=TSSEnrichment,x=log10(fragments))) + 
    geom_pointdensity() + theme_classic() + ylab('TSS Enrichment') +
    xlab(bquote(log[10] ~ "Number of Fragments")) + scale_color_viridis() + ggtitle(sample.name)
}


#' Plot Fragment size distribution
#'
#' Plot the fragment size distribution for a samples.
#'
#' @author Natalie Fox
#' @importFrom ggplot2 ggplot aes theme_classic xlab
#' @importFrom SummarizedExperiment colData
#' @export
plotFragmentSizeDistribution <- function(mae, sample.name) {

}

#' Plot Cluster Confusion Matrix
#'
#' Plot the association between slusters and samples using a confusion matrix.
#'
#' @author Natalie Fox
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment colData
#' @export
plotClusterConfusionMap <- function(mae) {
  tile.matrix.coldata <- colData(findSCE(mae,'TileMatrix500')$sce)
  cM <- table(tile.matrix.coldata[,c('Clusters','Sample')])
  cM <- cM / Matrix::rowSums(cM)
  p.confusion.map <- pheatmap(
    mat = as.matrix(cM),
    border_color = 'black'
  )
}
