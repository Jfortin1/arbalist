#' Add gene integration matrix
#'
#' Calculate gene integration matrix and add the matrix to the MultiAssayExperiment.
#' 
#' @param mae \linkS4class{MultiAssayExperiment} containing scATAC-seq gene score experiment and optionally the scRNA-seq experiment.
#' @param rna.sce \linkS4class{SingleCellExperiment} or \linkS4class{RangedSummarizedExperiment} or \linkS4class{RangedSummarizedExperiment} specifying scRNA-seq sata to integrated with the scATAC-seq data. This can be NULL if the RNA experiment is already in the mae.
#' @param gene.score.experiment.name String containing the experiment name of the gene score matrix calculated from ATAC data.
#' @param gene.expression.experiment.name String containing the experiment name of gene expression matrix (RNA data).
#' @param name.reduced.dim String containing the reduced dimension name. 
#' @param group.gex String specifying the column name for RNA groups. If the column is not in the gene.expression.experiment.name matrix, the specifying the correct experiment name as the name in a named vector.
#' @param group.list A list of strings specifying the cells to be used for RNA-ATAC integration. This is used to constrain integration across biologically relevant groups. Format of the list must be subgroups with RNA and ATAC sets. For example: list(groupA = list(ATAC = cellsATAC_A, RNA = cellsRNA_A), groupB = list(ATAC = cellsATAC_B, RNA = cellsRNA_B))
#' @param sampleCellsATAC Scalar integer specifying the number of scATAC-seq cells to use for integration.
#' @param sampleCellsRNA Scalar integer specifying the number of scRNA-seq cells to use for integration.
#' @param reduction String specifying the seurat method to use for integrating modalities. See Seuat::FindTransferAnchors() for details.
#' @param genesUse String vector specifying gene names to use for integration.
#' @param num.gene.features Scalar integer specifying the number of variable genes determined by `Seurat::FindVariableGenes()`
#'
#' @return \linkS4class{MultiAssayExperiment} with DoubletScore and DoubletEnrichment columns added to the experiment colData.
#'
#' @author Natalie Fox
#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<- rowRanges<- assay
#' @importFrom SingleCellExperiment mainExpName<-
#' @importFrom Seurat DefaultAssay CreateSeuratObject
#' @export
addGeneIntegrationMatrix <- function(
    mae,
    rna.sce = NULL,
    gene.score.experiment.name = 'GeneScoreMatrix',
    gene.expression.experiment.name = 'GeneExpressionMatrix',
    name.reduced.dim = c("TileMatrix500" = "iterativeLSI"),
    group.gex = c("GeneExpressionMatrix" = "Clusters"),
    group.list = NULL,
    sampleCellsATAC = 10000,
    sampleCellsRNA = 10000,
    reduction = "cca",
    genesUse = NULL,
    num.gene.features = 2000
) {

  # find the experiments in the MAE that we are going to use in 
  required.exp.names <- unique(c(
    gene.score.experiment.name,
    gene.expression.experiment.name,
    names(group.gex)
  ))
  sce.list <- list()
  if(!is.null(rna.sce)) {
    sce.list <- list('GeneExpressionMatrix' = list(sce = rna.sce, sce.idx = NULL, alt.exp.name = NULL))
    names(sce.list) <- gene.expression.experiment.name
  }
  for(i in setdiff(required.exp.names, names(sce.list))) {
    sce.list[[i]] <- findSCE(mae, i)
    if(is.null(sce.list[[i]])) {
      stop(paste0(i,' is not found in mae'))
    }
  }
  
  # retrieve the reduced dim for weighting
  reduced.dim.list <- list()
  for(i in 1:length(name.reduced.dim)) {
    reduced.dim.list[[i]] <- findReducedDimRes(mae,name.reduced.dim[i])
  }
  # if more than one reduced dimension specified then combine them into one set
  reduced.dim.matrix <- reduced.dim.list[[1]]$matrix
  if(length(reduced.dim.list) > 1) {
    colnames(reduced.dim.matrix) <- paste(names(name.reduced.dim[1]),name.reduced.dim[1],colnames(reduced.dim.matrix), sep="_")
    for(i in 2:length(reduced.dim.list)) {
      colnames(reduced.dim.list[[i]]$matrix) <- paste(names(name.reduced.dim[i]),name.reduced.dim[i],colnames(reduced.dim.list[[i]]$matrix), sep="_")
      if(nrow(reduced.dim.matrix) != nrow(reduced.dim.list[[i]]$matrix)) {
        stop('The reduced dimension matrices need to have the same number of rows. Try filtering your experiments to have the same number of cells.')
      }
      reduced.dim.matrix <- cbind(reduced.dim.matrix,reduced.dim.list[[i]]$matrix)
    }
    # Seurat does not like the reduced dimension names when combining matrices so renaming to all be the same.
    colnames(reduced.dim.matrix) <- paste0('LSI_',seq(1:ncol(reduced.dim.matrix)))
  }
  
  # if the experiment for the colData wasn't specified as a named vector then we will assume the columns are in the colData attached to the GeneExpression experiment
  if(is.null(names(group.gex))) {
    names(group.gex) <- gene.expression.experiment.name
  }

  # if groups are not specified than set up to have all cells in one group
  if(is.null(group.list)) {
    group.list <- SimpleList()
    group.list[['onlygroup']] <- SimpleList(
      ATAC = colnames(sce.list[[gene.score.experiment.name]]$sce),
      RNA = colnames(sce.list[[gene.expression.experiment.name]]$sce)
    )
  }
  
  # create the seurat object and add RNA data
  
  seuratRNA <- CreateSeuratObject(counts = assay(sce.list[[gene.expression.experiment.name]]$sce))
  seuratRNA$Group <- paste0(colData(sce.list[[names(group.gex)]]$sce)[, group.gex, drop = TRUE])
  Seurat::DefaultAssay(seuratRNA) <- "RNA"
  dfRNA <- DataFrame(row.names = colnames(seuratRNA), Group = seuratRNA$Group)
  seuratRNA <- seuratRNA[, unique(group.list$RNA)]
  seuratRNA <- Seurat::NormalizeData(object = seuratRNA, verbose = FALSE)
  
  # create integration blocks
  
  block.list <- SimpleList()
  
  for(i in seq_along(group.list)){
    
    gLi <- group.list[[i]]
    
    # ATAC
    
    if(length(gLi$ATAC) > sampleCellsATAC){
      
      cellsATAC <- sample(gLi$ATAC, length(gLi$ATAC))
      
      cutoffs <- unlist(lapply(seq_len(1000), function(x) length(gLi$ATAC) / x))
      blockSize <- ceiling(min(cutoffs[order(abs(cutoffs - sampleCellsATAC))[1]] + 1, length(gLi$ATAC)))
      
      # density Based Blocking
      nBlocks <- ceiling(length(gLi$ATAC) / blockSize)
      blocks <- SimpleList(lapply(seq_len(nBlocks), function(x){
        cellsATAC[seq(x, length(cellsATAC), nBlocks)]
      }))
      
    }else{
      blocks <- list(gLi$ATAC)
    }
  
    # RNA
    
    probRNA <- rep(1, length(gLi$RNA))
    block.list.i <- SimpleList(lapply(seq_along(blocks), function(x){
      SimpleList(
        ATAC = blocks[[x]],
        RNA = sample(x = gLi$RNA, size = min(sampleCellsRNA, length(gLi$RNA)) , prob = probRNA)
      )
    }))
    
    block.list <- c(block.list, block.list.i)
    
  }
  rm(group.list)
  
  ## integration
  sce.existing <- NULL
  for(i in seq_along(block.list)){
    
    blocki <- block.list[[i]]
    
    # subset RNA
    subRNA <- seuratRNA[, blocki$RNA]
    subRNA <- subRNA[rownames(subRNA) %in% rownames(sce.list[[gene.score.experiment.name]]$sce), ]
    subRNA <- Seurat::FindVariableFeatures(object = subRNA, nfeatures = num.gene.features, verbose = FALSE)
    subRNA <- Seurat::ScaleData(object = subRNA, verbose = FALSE)
    
    if(is.null(genesUse)){
      genesUse <- Seurat::VariableFeatures(object = subRNA)
    }
    
    # log-normalize
    mat <- assay(sce.list[[gene.score.experiment.name]]$sce)[,blocki$ATAC]
    if(is.null(rownames(mat))) {
      if(!is.null(rowRanges(sce.list[[gene.score.experiment.name]]$sce)$name)) {
        rownames(mat) <- rowRanges(sce.list[[gene.score.experiment.name]]$sce)$name
      } else {
        gene.grs <- rowRanges(sce.list[[gene.score.experiment.name]]$sce)
        rownames(mat) <- paste0(seqnames(gene.grs),':',start(gene.grs),'-',end(gene.grs))
      }
    }
    if(!is(mat,'sparseMatrix')) {
      mat <- as(mat,'sparseMatrix')
    }
    mat <- log(mat + 1) #use natural log
    seuratATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    
    # set default assay
    Seurat::DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = FALSE)
    
    # transfer anchors
    gc()
    transferAnchors <- Seurat::FindTransferAnchors(
      reference = subRNA, 
      query = seuratATAC, 
      reduction = reduction, 
      features = genesUse,
      verbose = FALSE
    )
    
    # create gene integration matrix and col data
    w <- Seurat::CreateDimReducObject(
      embeddings = reduced.dim.matrix[colnames(seuratATAC),,drop=FALSE], 
      key = "LSI_",
      assay = Seurat::DefaultAssay(seuratATAC)
    )
    rna.prediction.scores <- Seurat::TransferData(
      anchorset = transferAnchors,
      weight.reduction = w,
      verbose = FALSE,
      dims = seq_len(ncol(reduced.dim.matrix)),
      refdata = subRNA$Group
    )
    rna.labels <- Seurat::TransferData(
      anchorset = transferAnchors,
      weight.reduction = w,
      verbose = FALSE,
      dims = seq_len(ncol(reduced.dim.matrix)),
      refdata = colnames(subRNA)
    )[,1]
    matched.rna.matrix <- Seurat::TransferData(
      anchorset = transferAnchors,
      weight.reduction = w,
      verbose = FALSE,
      dims = seq_len(ncol(reduced.dim.matrix)),
      refdata = Seurat::GetAssayData(subRNA, assay = "RNA", slot = "data")
    )@data
    matchDF <- DataFrame(
      Sample = colnames(seuratATAC), 
      prediction_score = rna.prediction.scores$prediction.score.max,
      predicted_group = rna.prediction.scores$predicted.id,
      rna_predicted_cell = rna.labels
    )
    
    # create GeneIntegrationMatrix experiment for this block
    mat.list <- list(counts=matched.rna.matrix)
    se <-  SummarizedExperiment(mat.list)
    rowRanges(se) <- rowRanges(sce.list[[gene.expression.experiment.name]]$sce)[rownames(matched.rna.matrix)]
    rowData(se) <- rowData(sce.list[[gene.expression.experiment.name]]$sce)[rownames(matched.rna.matrix),]
    colData(se) <- matchDF
    colnames(se) <- colnames(matched.rna.matrix)
    sce <- as(se, 'SingleCellExperiment')
    mainExpName(sce) <- 'GeneIntegrationMatrix'
    
    # combine experiment amoung blocks
    if(is.null(sce.existing)) {
      sce.existing <- sce
    } else {
      sce.existing <- cbind(sce.existing, sce)
    }
  }

  # add the GeneIntegrationMatrix Experiment to the MAE
  
  exp.list <- experiments(mae)
  exp.list[['GeneIntegrationMatrix']] <- sce.existing
  
  el <- ExperimentList(exp.list)
  sampleLabels <- 'Sample'
  maplist <- lapply(exp.list, function(se) {
    if(sampleLabels %in% colnames(colData(se))) {
      data.frame(primary = colData(se)[,sampleLabels], colname = colnames(se), stringsAsFactors = FALSE)
    } else {
      data.frame(primary = se$ID, colname = colnames(se), stringsAsFactors = FALSE)
    }
  })
  sampMap <- listToMap(maplist)
  
  new.mae <- MultiAssayExperiment(el, sampleMap = sampMap, colData = DataFrame(row.names=unique(sampMap$primary)))
  colData.filler <- matrix(NA,ncol=ncol(colData(mae)),nrow=ncol(sce.existing))
  colnames(colData.filler) <- colnames(colData(mae))
  rownames(colData.filler) <- colnames(sce.existing)
  
  colData(new.mae) <- rbind(colData(mae),colData.filler)
  metadata(new.mae) <- metadata(mae)

  return(new.mae)
}
