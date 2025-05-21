\name{SemiAutoAnnote}
\alias{SemiAutoAnnote}
\title{Annotate: Semi-Automated Annotation By User Input}
\usage{
SemiAutoAnnotate(ann,
	data_type = NULL,
	refine = NULL)
}
\description{
SemiAutoAnnotation requires a SAHA object (ann) and data_type as inputs to function. This tool is used to view cluster identification quickly to aid in annotation. Based on the data_type selected, either Marker-based, Marker-free, or both workflows must be completed prior to semi-auto annotation.

If refine = 'null', then the function will ask the user to label every cluster in the query dataset. If the user would like to rely on automated annotations, but have the choice to refine annotations that do not agree between Marker-based and Marker-free, the output of AutoAnnotate() must be placed into the refine argument. 
}
\arguments{
\item{ann}{Required argument that takes a SAHA object generated after running either Create_MarkerFree_Viz(), Create_MarkerBased_Viz() or both.}
\item{data_type}{Required argument that specifies which type of data the user is querying. The three available options:
* 'Markers': Iterate through Marker-based dotplot on a per-cluster basis. 
* 'AvgExp': Iterate through Marker-free heatmap on a per-cluster basis.
* 'Both': Iterate through both the Marker-based dotplot and Marker-free heatmap on a per-cluster basis.
}
\item{refine}{Optional argument that specifies how the user will re-name their clusters. If left NULL, SemiAutoAnnotate() will prompt each cluster in the query dataset. If a dataframe produced by AutoAnnotate() is provided into this argument, only clusters with best_match = "INCONCLUSIVE" will be shown.}
}
\examples{
#SemiAutoAnnotate all clusters based on the Marker-based 
semi_auto=SemiAutoAnnotate(ann,data_type = "Markers")

#SemiAutoAnnotate only INCONCLUSIVE clusters based on Marker-based and Marker-free
auto=AutoAnnotate(ann,data_type="Both")
semi_auto=SemiAutoAnnotate(ann,data_type = "Both", refine = auto)

}