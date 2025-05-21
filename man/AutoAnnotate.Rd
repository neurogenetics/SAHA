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
  \item{\strong{"Markers"}: }{Uses marker enrichment data (\code{ann@ann2}) to assign the best-matching cell type with the lowest p-value in each cluster (p < 0.05). Returns a data frame with cluster, best match, proportion, and p-value. Clusters without a significant match are labeled as \code{"INCONCLUSIVE"}.}
  \item{\strong{"AvgExp"}: }{Uses correlation between average expression profiles from \code{ann@results$marker_free$corr} to assign the reference cell type with the highest correlation to each cluster. Returns a data frame with cluster, best match, and correlation value.}
  \item{\strong{"Both"}: }{Compares both methods. If the best matches from the marker-based and correlation-based methods agree, the consensus is matched (case- specific) and the final assignment is the agreed-upon cell type. If they differ, the consensus is \code{"DISAGREEMENT"}. If either of the methods have an \code{"INCONCLUSIVE"} annotation, this result will override all other actions in the final dataframe.}
}

If \code{data_type} is not specified, a message will be printed discouraging blind auto-annotation and encouraging thoughtful review.
}

\arguments{
  \item{ann}{A required SAHA object containing marker-based results and/or correlation data.
  \item{data_type}{Specifies which method used for auto-annotation. Must be one of \code{"Markers"}, \code{"AvgExp"}, or \code{"Both"}.}
}

\value{
Returns a data frame with marker-based, marker-free, consensus, and final annotation per cluster. The structure depends on the selected \code{data_type}.
}

\seealso{
For manual adjustment and refinement of annotations, see \code{\link{SemiAutoAnnotation}}.
}

\examples{
# Marker-based annotation
auto <- AutoAnnotate(ann, data_type = "Markers")

# Correlation-based annotation
auto <- AutoAnnotate(ann, data_type = "AvgExp")

# Consensus from both methods
auto <- AutoAnnotate(ann, data_type = "Both")
}
