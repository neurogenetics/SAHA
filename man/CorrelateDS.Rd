\name{CorrelateDS}
\alias{CorrelateDS}
\title{Marker-free annotation: Correlation}
\description{
Performs correlation between normalized gene expression profiles of every cell type in the database to every unidentified cluster in the query dataset.
}
\usage{
CorrelateDS(ann, corr_method = "pearson")
}
\arguments{
  \item{ann}{A SAHA object which has been normalized using NormalizeDS() for marker-free annotation.}
  \item{corr_method}{Character string indicating which correlation coefficient is to be computed. Options are "pearson" (default), "kendall", or "spearman".}
}
\value{
Returns the modified \code{ann} object with the correlation data frame stored within \code{ann@results$marker_free$corr}.
}
\details{
This function calculates correlation coefficients between all query and database cell type expression profiles stored in the \code{ann@results$marker_free$norm_merge} data frame within the \code{ann} object. The results are stored in a data frame named \code{corr} within \code{ann@results$marker_free}.
}
\examples{
# Using default Pearson correlation
CorrelateDS(ann)

# Using Spearman correlation
CorrelateDS(ann, corr_method = "spearman")
}
