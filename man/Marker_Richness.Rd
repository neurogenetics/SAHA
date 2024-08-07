\name{Marker_Richness}
\alias{Marker_Richness}
\title{Marker Richness: Number of Markers Per Cluster}
\usage{
Marker_Richness(ann,
   varfeat = NULL)
}
\description{
Creates a barplot to visualize the number of markers in the ann@ann1 (initialized and possibly tuned markers from ann@query$Markers).

Not intended to be written into the SAHA object and must be re-run to create the plot.
}
\arguments{
\item{ann}{A SAHA object.}
\item{varfeat}{An optional input containing a string of genes named with the same convention as rownames of query. Recommended usage would be to write the variable features from a Seurat object so that the number of markers could be visualized with respect to the variable features used to create clusters.}
}
\examples{
#Visualize all markers per cluster
MarkerRichness(ann)

##optional plotting of varfeat
t1=MarkerRichness(ann,varfeat = varfeat)
ggplot(t1,aes(x=cluster, y = genes_in_varfeat))+
  geom_bar(stat="identity")+
  geom_hline(yintercept = length(varfeat),linetype = 3)+
  theme_bw()

}
