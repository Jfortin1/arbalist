#' Filter duplicate features
#' 
#' Keeps one feature from each set of duplicated features based on specified criteria, such as the feature with the highest values. Other duplicate features are removed.
#'
#' @param se \linkS4class{SummarizedExperiment}
#' @param mcol.name String specifying the colname for the experiment rowData
#' @param summary.stat Function to summarize each feature (row) of the experiment
#' @param selection.metric Function to select the row to keep when there are duplicate rows with the same mcol.name
#'
#' @importFrom SummarizedExperiment mcols
#' @export
filterDuplicateFeatures <- function(se, mcol.name = 'name', summary.stat = sum, selection.metric = max) {

  duplicate.values <- names(which(table(mcols(se)[,mcol.name]) > 1))
  non.duplicate.rows <- which(! mcols(se)[,mcol.name] %in% duplicate.values)
  duplicate.rows <- which(mcols(se)[,mcol.name] %in% duplicate.values)
  selected.duplicate.rows <- sapply(duplicate.values,function(i) {
    duplicate.rows <- which(mcols(se)[,mcol.name] %in% i)
    row.summary.stats <- apply(assay(se)[duplicate.rows,],1,summary.stat)
    return(duplicate.rows[which(row.summary.stats == selection.metric(row.summary.stats))[1]])
  })
  
  se <- se[sort(c(non.duplicate.rows,selected.duplicate.rows)),]
  
  return(se)
  
}