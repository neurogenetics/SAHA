\name{Investigate_Marker_Based}
\alias{Investigate_Marker_Based}
\title{Semi-Automated Hand Annotation: Investigate_Marker_Based}
\usage{
Initialize_Markers(ann, 
	query_cluster, 
	db_cell_type, 
	plot)
}
\description{
Investigate_Marker_Based is used to investigate specific clusters of a given SAHA object after Marker-based visualizations have been created. This function can compare one cluster to a all available annotations using plot = `stacked` or compared to a specific annotaiton using plot = `venn`. Additionally, this function will print out overlapping markers from the db and query dataset within the R console. To run Investigate_Marker_Based(), make sure you have successfully ran Run_Marker_Based(). 
}
\arguments{
\item{ann}{Required argument(SAHA object) that is used to create cluster-specific visualizations.}
\item{query_cluster}{Required arugment that specifies a given cluster.}
\item{db_cell_type}{Optional arugment (for plot = 'venn') to make comparisons between specific cell types within the marker database given clusters.}
\item{plot}{Required argument to pick which type of investigation to run. Two available options:
* "stacked": a stacked barplot showing the amount of markers within each cluster found within the database.
* "venn": a venn diagram comparing the amount of overlapped markers between a specific cluster and the database.}
}
\examples{
#See overlap between cluster 4 and all possible annotations
Investigate_Marker_Based(ann,query_cluster = 4,plot = "stacked")
#See overlap between cluster 4 and annotation "295 CBN Dmbxl Gaba"
Investigate_Marker_Based(ann,query_cluster = 4,plot = "venn",db_cell_type = "295 CBN Dmbx1 Gaba")
}