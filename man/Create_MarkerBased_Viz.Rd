\name{Create_MarkerBased_Viz}
\alias{Create_MarkerBased_Viz}
\title{Creates marker-based visualizations for cell type annotation}
\description{
Generates various visualizations based on marker-based cell type annotation results.
Creates dot plots and a bar plot, with optional faceting by cell type class.
}
\usage{
Create_MarkerBased_Viz(ann, meta = NULL)
}
\arguments{
  \item{ann}{An SAHA analysis object containing annotation results (likely from \code{Annotate()} or similar).}
  \item{meta}{Optional data frame containing cell type metadata. If provided, enables faceting of the plots by cell type class.
    The data frame must include the following columns:
    \itemize{
      \item{\code{subclass_spa}:}{ Cell type names matching those in the \code{ann} object.}
      \item{\code{class}:}{ Broader cell type classification.}
      \item{\code{neurotransmitter}:}{ Even broader cell type classification.}
      \item{\code{class_color}:}{ Hex code for the color of the \code{class} classification.}
      \item{\code{neurotransmitter_color}:}{ Hex code for the color of the \code{neurotransmitter} classification.}
    }
  }
}
\value{
Returns the \code{ann} object, with the following plots added to the \code{results$marker_based} slot:
  \itemize{
    \item{\code{markers_barplot}:}{ A bar plot showing the total number of markers per cell type.}
    \item{\code{dotplot_all}:}{ A dot plot showing all markers for each cell type and cluster.}
    \item{\code{dotplot_sig}:}{ A dot plot showing only significant markers for each cell type and cluster.}
    \item{\code{dotplot_best}:}{ A dot plot showing the most significant marker for each cell type and cluster.}
  }
}
\details{
This function generates visualizations to aid in the interpretation of marker-based cell type annotation results. If the \code{meta} data frame is provided, the dot plots will be faceted by cell type class, providing additional context for the marker expression patterns.  The color scheme for the facets is determined by the \code{class_color} column in the \code{meta} data frame.  The barplot shows the total number of markers available for each cell type.  Dot plots show marker enrichment and proportion, both for all markers, significantly enriched markers, and the top significantly enriched marker.
}
\examples{
# Basic usage without metadata
# ann <- Create_MarkerBased_Viz(ann)

# With metadata for faceting
# ann <- Create_MarkerBased_Viz(ann, meta)
}
