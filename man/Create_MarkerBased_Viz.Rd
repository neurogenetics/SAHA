\name{Create_MarkerBased_Viz}
\alias{Create_MarkerBased_Viz}
\title{Create: Marker-based annotation plots}
\usage{
Create_Marker_Based(ann,
	meta,
	facet)
}
\description{
Create_MarkerBased_Viz requires a SAHA object to create 4 different visualizations that are used to classify the user's clusters using a marker-based database. The results are then stored within the `results` slot of the SAHA object. This function also allows the option to turn on/off faceting for the figures:
* Plot 1 Description: Ggplots all insignificant and significantly annotated clusters.
* Plot 2 Description: Barplots the overall distribution of markers in the inputted database.
* Plot 3 Description: Ggplots only the significantly annotated clusters. 
* Plot 4 Description: Ggplots only the most significant annotated clusters. 

This function is intended to write directly into a SAHA object (see usage below).
}
\arguments{
\item{ann}{Required argument that takes a SAHA object generated from Run_Marker_Based().}
\item{meta}{Optional argument to load in dataframe providing additional information about the cell types. This metadata is used to initiate the main data frame (master_df) with class, class color, neurotransmitter, and neurotransmitter color information for each cell type.}
\item{facet}{Optional argument that can specify to turn on or off the faceting feature when creating the visualiations.}
}
\examples{
ann=Create_MarkerBased_Viz(ann,meta = meta,facet = TRUE)
}