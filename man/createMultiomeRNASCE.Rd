% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createMultiomeRNASCE.R
\name{createMultiomeRNASCE}
\alias{createMultiomeRNASCE}
\title{Import gene expression results from multiome Cell Ranger results}
\usage{
createMultiomeRNASCE(
  h5.files,
  sample.names = NULL,
  feature.type = "Gene Expression",
  filter.features.without.intervals = TRUE
)
}
\arguments{
\item{h5.files}{Vector of strings specifying filtered_feature_bc_matrix.h5 path. ex. could just be filtered_feature_bc_matrix.h5. Vector must be the same length as sample.names.}

\item{sample.names}{Vector of strings specifying sample names. Vector must be the same length as h5.files.}

\item{feature.type}{String specifying the feature type to select from filtered_feature_bc_matrix.h5.}

\item{filter.features.without.intervals}{Logical whether to remove features from the h5.files that do not have interval specified. Often these are mitochondria genes.}
}
\value{
A \linkS4class{SingleCellExperiment}
}
\description{
Import RNA results from multiome Cell Ranger results into a SingleCellExperiment.
}
\author{
Natalie Fox
}
