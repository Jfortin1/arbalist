% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addGeneScoreMatrix.R
\name{addGeneScoreMatrix}
\alias{addGeneScoreMatrix}
\title{Add gene score matrix}
\usage{
addGeneScoreMatrix(
  mae,
  gene.grs,
  experiment.name = "TileMatrix500",
  extend.outside.gene = c(1000, 1e+05),
  gene.model = "exp(-abs(x)/5000) + exp(-1)",
  gene.scale.factor = 5,
  upstream = 5000,
  downstream = 0,
  num.threads = 4
)
}
\arguments{
\item{mae}{\linkS4class{MultiAssayExperiment}.}

\item{gene.grs}{String specifying the experiment name to use rowRanges start as for TSS. If this experiment name is not found in mae, then gene.grs must be specified. If both experiment.name.for.gene.grs and gene.grs are specified, then gene.grs takes precedence.}

\item{experiment.name}{String containing the experiment name of the tile matrix.}

\item{extend.outside.gene}{A numeric vector with two values specifying the range to extend outside the gene prior to overlapping genes with tiles.}

\item{gene.model}{A string giving a "gene model function" used for weighting peaks for gene score calculation. This string
should be a function of `x`, where `x` is the stranded distance from the transcription start site of the gene.}

\item{gene.scale.factor}{A numeric scaling factor to weight genes based on the inverse of there length i.e. [(Scale Factor)/(Gene Length)]. This
is scaled from 1 to the scale factor. Small genes will be the scale factor while extremely large genes will be closer to 1. This scaling helps with
the relative gene score value.}

\item{upstream}{Integer describing the number of bp upstream the gene to extend the gene body. This effectively makes the gene body larger as there are proximal peaks that should be weighted equally to the gene body.}

\item{downstream}{Integer describing the number of bp downstream the gene to extend the gene body.This effectively makes the gene body larger as there are proximal peaks that should be weighted equally to the gene body.}

\item{num.threads}{Integer scalar specifying the number of threads to use for tatami multiply.}
}
\value{
\linkS4class{MultiAssayExperiment} with DoubletScore and DoubletEnrichment columns added to the experiment colData.
}
\description{
Calculate gene score matrix to estimate gene expression and add the matrix to the MultiAssayExperiment.
}
\author{
Natalie Fox
}
