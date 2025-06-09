\name{Initialize_Self_Similarity}
\alias{Initialize_Self_Similarity}
\title{Initialize Self-Similarity Matrices for SAHA Analysis}
\usage{
  Initialize_Self_Similarity(ann, slot, custom_ds = NULL)
}
\description{
  Initializes and stores a self-similarity matrix within the \code{ann} object based on either average expression or marker incidence. This is used for downstream similarity and clustering analysis in the SAHA framework.
}
\arguments{
  \item{ann}{An object containing query data and annotations, typically output from earlier steps in the SAHA pipeline.}
  \item{slot}{A character string specifying which type of similarity to compute. Must be either \code{"AvgExp"} for average expression or \code{"Markers"} for binary marker presence across clusters.}
  \item{custom_ds}{Optional. A character vector of gene names to subset the similarity matrix calculation (only applies to \code{slot = "AvgExp"}).}
}
\details{
  \itemize{
    \item If \code{slot = "AvgExp"}, a correlation matrix is computed across clusters based on their average gene expression.
    \item If \code{slot = "Markers"}, a binary incidence matrix is constructed to indicate gene presence across clusters, then converted to a wide format.
  }

  The resulting matrix is stored in \code{ann@ann1$AvgExp_selfsim_matrix} or \code{ann@ann1$Marker_selfsim_matrix}, depending on the chosen slot.

  The function also appends a record of the operation to \code{ann@params$markers$selfsim_analysis} for reproducibility.
}
\value{
  An updated \code{ann} object with the corresponding self-similarity matrix stored in the appropriate slot.
}
\examples{
  # Initialize self-similarity using average expression
  ann <- Initialize_Self_Similarity(ann, slot = "AvgExp")
  
  # Initialize self-similarity using markers
  ann <- Initialize_Self_Similarity(ann, slot = "Markers")
  
  # Using a custom gene set for AvgExp similarity
  genes_of_interest <- c("GeneA", "GeneB", "GeneC")
  ann <- Initialize_Self_Similarity(ann, slot = "AvgExp", custom_ds = genes_of_interest)
}
