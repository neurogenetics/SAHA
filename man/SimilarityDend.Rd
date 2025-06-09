\name{SimilarityDend}
\alias{SimilarityDend}
\title{Generate a Dendrogram Based on Self-Similarity}
\usage{
  SimilarityDend(ann, sim_type, annotation)
}
\description{
  This function generates a dendrogram based on a self-similarity matrix derived from either marker-based or average expression data. Optionally, it renames the dendrogram labels using a provided annotation data frame.
}
\arguments{
  \item{ann}{An object containing SAHA analysis results, including self-similarity matrices.}
  \item{sim_type}{The type of similarity matrix to use. Must be either \code{"Markers"} or \code{"AvgExp"}.}
  \item{annotation}{Optional. A data frame used to relabel the dendrogram tips. Must contain a \code{cluster} column and a \code{best_match} column.}
}
\details{
  This function extracts the appropriate similarity matrix based on the \code{sim_type} parameter, generates a dendrogram using hierarchical clustering, and optionally renames the labels for clarity using the provided \code{annotation} data frame.
}
\value{
  A dendrogram plot is generated. The function is primarily used for visualization and does not return a value.
}
\examples{
  # Generate a dendrogram using marker-based similarity
  SimilarityDend(ann, sim_type = "Markers")
  
  # Generate a dendrogram with annotated cluster labels
  SimilarityDend(ann, sim_type = "AvgExp", annotation = cluster_annotations)
}
