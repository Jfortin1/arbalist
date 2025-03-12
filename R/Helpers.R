#' Return SingleCellExperiment with the experiment name
#' 
#' Checks through the MultiAssayExperiment both the main and alternative experiments and returns the experiment with the specified name. 
#'
#' @param mae \linkS4class{MultiAssayExperiment}
#' @param experiment.name String specifying the the experiment name to return. Name can either be in names() of the MAE or the altExp name in one of the experiments.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment altExp
#' @importFrom methods is
#' @export
findSCE <- function(mae, experiment.name) {
  # Find the experiment result
  sce <- NULL;
  sce.idx <- NULL;
  alt.exp.name <- NULL;
  if(experiment.name %in% names(mae)) {
    sce.idx <- which(names(mae) == experiment.name)
    sce <- mae[[sce.idx]]
  } else {
    for(i in seq_len(length(mae))) {
      if(is(mae[[i]],'SingleCellExperiment') && experiment.name %in% altExpNames(mae[[i]])) {
        sce <- altExp(mae[[i]],experiment.name)
        # record where the iterative LSI was saved so we can put the UMAP in the same place
        sce.idx <- i
        alt.exp.name <- experiment.name
        if (ncol(colData(mae[[sce.idx]])) > 0) {
          colname.overlap.count <- table(c(colnames(colData(mae[[sce.idx]])),colnames(colData(sce))))
          if(ncol(colData(sce)) == 0) {
            colData(sce) <- colData(mae[[sce.idx]])
          } else if(!any(colname.overlap.count > 1)) {
            colData(sce) <- cbind(colData(mae[[sce.idx]]),colData(sce))
          } else {
            colData(sce) <- cbind(colData(mae[[sce.idx]])[,-which(colname.overlap.count == 2)],colData(sce))
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

#' Return reduced dimensions matrix with the reduced dimension name
#' 
#' Checks through the MultiAssayExperiment both the main and alternative experiments and returns the reduced dimension with the specified name. 
#'
#' @param mae \linkS4class{MultiAssayExperiment}
#' @param name.reduced.dim String specifying the the reduced dimensions name to return. String can be a named vector where the name in the experiment in case there are reduced dimensions within multiple experiments with the same names.
#'
#' @importFrom SingleCellExperiment altExp altExpNames reducedDim reducedDimNames
#' @export
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
    stop(name.reduced.dim,' is not found in mae')
  }
  return(list(reduced.dim.name=name.reduced.dim, matrix=dim.matrix, exp.idx=dim.exp.idx, alt.exp.name=dim.alt.exp.name))
}
