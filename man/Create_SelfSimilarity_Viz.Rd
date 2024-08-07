\name{Create_SelfSimilarity_Viz}
\alias{Create_SelfSimilarity_Viz}
\title{Create: Heatmaps to Visualize Similarity Between Clusters in Query}
\usage{
Create_SelfSimilarity_Viz(ann,
   slot,
   assay_db="RNA")
}
\description{
Create heatmaps for similarity between all clusters of the query and store them in ann@results. When slot = "Markers" is sepecified the heatmap will represent the overlapping number of markers plotted on a log scale. When slot = "AvgExp" is specified the heatmap will plot the R^2 between given clusters based on all genes in ann@db$AvgExp.

Intended to be written directly into the ann object loaded.
}
\arguments{
\item{ann}{A SAHA object.}
\item{slot}{Specifies whether the self-similarity matrix is created from:
  * "Markers": ann@query$Markers, a loaded df containing marker genes formatted as output from Seurat::FindAllMarkers()
  * "AvgExp": ann@query$AvgExp, a loaded df containing avgerage expression formatted as output from Seurat::AverageExpression()
}
\item{assay_db}{Assay name in Seurat object. This is used to clip the beginning of the default names off the clusters in the AvgExp dataframe (E.g. RNA.g.0 -> 0). Options may include:
  * "RNA": (recommended) created from the RNA counts matrix in a Seurat object
  * "SCT": created from SCTransform corrected counts in a Seurat object
  * "Xenium": created from quantifications of Xenium assay type in a Seurat object
}
}
\examples{
#Create marker-based self similarity
ann=Create_SelfSimilarity_Viz(ann,slot = "Markers")
}
