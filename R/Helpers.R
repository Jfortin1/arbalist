#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment altExp
#' @export
findSCE <- function(mae,experiment.name) {
  # Find the experiment result
  sce <- NULL;
  sce.idx <- NULL;
  alt.exp.name <- NULL;
  if(experiment.name %in% names(mae)) {
    sce.idx <- which(names(mae) == experiment.name)
    sce <- mae[[sce.idx]]
  } else {
    for(i in seq_len(length(mae))) {
      if(experiment.name %in% altExpNames(mae[[i]])) {
        sce <- altExp(mae[[i]],experiment.name)
        # record where the iterative LSI was saved so we can put the UMAP in the same place
        sce.idx <- i
        alt.exp.name <- experiment.name
        if (ncol(colData(mae[[sce.idx]])) > 0) {
          if(ncol(colData(sce)) == 0) {
            colData(sce) <- colData(mae[[sce.idx]])
          } else {
            colData(sce) <- cbind(colData(mae[[sce.idx]])[,-which(table(c(colnames(colData(mae[[sce.idx]])),colnames(colData(sce)))) == 2)],colData(sce))
          }
        }
        break
      }
    }
  }
  if(is.null(sce.idx)) {
    return(NULL)
  }
  return(list(sce=sce, sce.idx=sce.idx, alt.exp.name=alt.exp.name))
}

#' @importFrom SingleCellExperiment altExp altExpNames reducedDim reducedDimNames
findReducedDimRes <- function(mae, name.reduced.dim) {
  if(length(name.reduced.dim) > 1) {
    stop('Please only specify one name.reduced.dim')
  }
  # Find the IterativeLSI result
  dim.matrix <- NULL;
  dim.exp.idx <- NULL;
  dim.alt.exp.name <- NULL;
  exp.to.check <- seq_len(length(mae))
  if(!is.null(names(name.reduced.dim))) {
    exp.to.check <- which(names(mae) == names(name.reduced.dim))
  }
  for(exp.idx in exp.to.check) {
    if(name.reduced.dim %in% reducedDimNames(mae[[exp.idx]])) {
      dim.matrix <- reducedDim(mae[[exp.idx]], name.reduced.dim)
      # record where the iterative LSI was saved so we can put the UMAP in the same place
      dim.exp.idx <- exp.idx
      break
    }
    if(is.null(dim.matrix)) {
      for(name.alt.exp in altExpNames(mae[[1]])) {
        if(name.reduced.dim %in% reducedDimNames(altExp(mae[[exp.idx]], name.alt.exp))) {
          dim.matrix <- reducedDim(altExp(mae[[exp.idx]], name.alt.exp), name.reduced.dim)
          # record where the iterative LSI was saved so we can put the UMAP in the same place
          dim.exp.idx <- exp.idx
          dim.alt.exp.name <- name.alt.exp
          break
        }
      }
    }
  }
  if(is.null(dim.matrix)) {
    stop(paste0(name.reduced.dim,' is not found in mae'))
  }
  return(list(reduced.dim.name=name.reduced.dim, matrix=dim.matrix, exp.idx=dim.exp.idx, alt.exp.name=dim.alt.exp.name))
}
