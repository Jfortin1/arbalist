% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getExpListFromFragments.R
\name{getExpListFromFragments}
\alias{getExpListFromFragments}
\title{Create a list of single cell experiments}
\usage{
getExpListFromFragments(
  fragment.files,
  output.dir = tempdir(),
  tile.size = 500,
  seq.lengths = NULL,
  gene.grs = NULL,
  barcodes.list = NULL,
  BPPARAM = bpparam()
)
}
\arguments{
\item{fragment.files}{Vector of strings specifying fragment files. Vector names need to be sample names.}

\item{output.dir}{String containing the directory where files should be output while creating the \linkS4class{MultiAssayExperiment}.}

\item{tile.size}{Integer scalar specifying the size of the tiles in base pairs.}

\item{seq.lengths}{Named integer vector containing the lengths of the reference sequences used for alignment.
Vector names should correspond to the names of the sequences, in the same order of occurrence as in the fragment file.
If \code{NULL}, this is obtained from the reference genome used by Cellranger (itself located by scanning the header of the fragment file).}

\item{gene.grs}{Genomic Ranges specifying gene coordinates for creating the gene score matrix. If NA, the gene accessibility matrix will not be created.}

\item{barcodes.list}{A List with samples as names and the values a vector of barcodes for that sample. If NULL, all barcodes from the fragment file will be used.}

\item{BPPARAM}{A \linkS4class{BiocParallelParam} object indicating how matrix creation should be parallelized.}
}
\value{
A list of experiments
}
\description{
Create a list of single cell experiments from fragment files
}
\author{
Natalie Fox
}
