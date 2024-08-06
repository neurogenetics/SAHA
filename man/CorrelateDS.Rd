\name{CorrelateDS}
\alias{CorrelateDS}
\title{Marker-free annotation: Correlation}
\description{Performs Pearson correlation between normalized gene expression profiles of every cell type in the database to every unidentified cluster in the query dataset.}
\usage{
CorrelateDS(ann)
}
\arguments{
  \item{ann}{A SAHA object which has been normalized using NormalizeDS() for marker-free annotation.}
}