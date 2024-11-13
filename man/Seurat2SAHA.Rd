\name{Seurat2SAHA}
\alias{Seurat2SAHA}
\title{Create SAHA inputs from a Seurat Objet}
\usage{
SAHA_input = Seurat2SAHA(obj, output="Both")
}
\description{
Function that accepts a seurat object (obj) and returns a list of outputs with a list of Variable Features, a dataframe of average expression, and a dataframe of markers.
}
\arguments{
\item{obj}{Required Seurat object, after clustering was completed. In order to receive all relevant input for SAHA, clustering must already be completed. By default it will pull all information from the RNA slot. Please refer to the "creating inputs" vignette.}
\item{outputs}{The specified outputs to be created, options include "Markers", "AvgExp", or "Both" (AvgExp + Markers). The function will also return Variable Features of the Seurat object.}
}
\examples{
SAHA_input = Seurat2SAHA(obj, output="Both")
SAHA_input = Seurat2SAHA(obj, output="AvgExp")
SAHA_input = Seurat2SAHA(obj, output="Markers")

##accessing inputs to be used in other SAHA functions
Create_SAHA_Object(...,query_AvgExp=SAHAinput$avgexp)
Create_SAHA_Object(...,query_Markers=SAHAinpust$markers)
Marker_Richness(...,varfeat = SAHAinpust$varfeat)

}
