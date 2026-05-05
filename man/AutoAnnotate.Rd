\name{AutoAnnotate}
\alias{AutoAnnotate}
\title{Automatic Annotation of Query Clusters}
\usage{
AutoAnnotate(ann,
             data_type = NULL)
}
\description{
Automatically annotates query clusters in a SAHA object using either marker-based statistics, average expression correlation, or a combination of both.

This function helps quickly identify likely cell types for each cluster to assist in downstream annotation. The user must specify the type of annotation method using the \code{data_type} argument:

\itemize{
  \item \strong{"Markers"}: Uses marker enrichment data (\code{ann@ann2}) to assign the best-matching cell type with the lowest p-value in each cluster (p < 0.05). Clusters without a significant match are labeled as \code{"INCONCLUSIVE"}.
  
  \item \strong{"AvgExp"}: Uses correlation between average expression profiles from \code{ann@results$marker_free$corr} to assign the reference cell type with the highest correlation.
  
  \item \strong{"Both"}: Compares both methods. If both agree, the consensus is assigned. If they differ, the result is \code{"DISAGREEMENT"}. If either method is \code{"INCONCLUSIVE"}, that overrides the final result.
}

If \code{data_type} is not specified, a message will be printed discouraging blind auto-annotation.
}

\arguments{
  \item{ann}{A SAHA object containing marker-based results and/or correlation data.}
  \item{data_type}{One of \code{"Markers"}, \code{"AvgExp"}, or \code{"Both"}.}
}

\value{
A data frame containing annotation results per cluster.
}

\seealso{
\code{\link{SemiAutoAnnotation}}
}

\examples{
# Marker-based annotation
auto <- AutoAnnotate(ann, data_type = "Markers")

# Correlation-based annotation
auto <- AutoAnnotate(ann, data_type = "AvgExp")

# Consensus from both methods
auto <- AutoAnnotate(ann, data_type = "Both")
}
