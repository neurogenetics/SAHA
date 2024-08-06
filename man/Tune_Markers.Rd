\name{Tune_Markers}
\alias{Tune_Markers}
\title{Marker-based annotation: Tuning markers}
\description{Thresholds markers for query or database for marker-based annotation. Should be run separately for database and query markers.}
\usage{Tune_Markers(ann, 
method, 
method_value, 
method_var,
set)
}
\arguments{
  \item{ann}{SAHA object which has been initialized using Initialize_Markers() for marker-based annotation.}
  \item{method}{Method for thresholding markers. Options are "absolute" (i.e., top X number of markers) or "relative" (i.e., above Xth percentile of markers).}
  \item{method_value}{Value to use for thresholding markers. If method = "absolute", represents the cutoff in ranked marker list (i.e., 100 means use top 100 markers). If method = "relative", should be a number between 0 and 1 (i.e., 0.75 means use only above 75th percentile of markers).}
  \item{method_var}{Method by which to rank markers. Can be "avg_log2FC" or "p_val_adj".}
  \item{set}{Dataset to threshold ("query" or "db").}
}

\examples{
Tune_Markers(ann, method = "absolute", method_value = 100, method_var = "avg_log2FC", set = "db")
Tune_Markers(ann, method = "relative", method_value = 0.75, method_var = "avg_log2FC", set = "query")
}