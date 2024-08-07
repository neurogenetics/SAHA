\name{Initialize_Markers}
\alias{Initialize_Markers}
\title{Initialize: Create Thresholds for Markers}
\usage{
Initialize_Markers(ann, 
	p_thresh = 0.05, 
	FC_thresh = 1.5, 
	sens_thresh = 0.25, 
	spec_thresh = 0.75)
}

\description{
Initialize_Markers is used for quality control of the user's SAHA object (ann) after the object has been created. While the SAHA object will contain unfiltered markers (i.e. ann@query$Markers), this function creates the marker dataframes to be used in marker-based annotation and stores in ann2 slot (i.e. ann@ann2$query). There are recommended defaults for each argument as well as suggestions throughout the code itself. To run Initialize_Markers, the user will have needed to first create their SAHA object using Create_SAHA_Object(). Please check that the colnames of both the query and db match the expected (standard output from Seurat::FindAllMarkers) if errors occur.

This function is meant to be written directly into an existing SAHA object.
}
\arguments{
\item{ann}{Required SAHA object, containing Marker dataframes in both the ann@query and ann@db slots. The dataframes should contain columns `avg_log2FC`, `cluster`, `pct.1`, and `pct.2`, and `p_val_adj` that are used for thresholding.}
\item{p_thresh}{Optional arugment that sets a p-value threshold (`p_val_adj` column) on both query and db Marker dataframes. The default is set to 0.05}
\item{FC_thresh}{Optional argument that sets a log base 2 Fold Change threshold (`avg_log2FC` column) on both query and db Marker dataframes. The default is set to 1.5.}
\item{sens_thresh}{Optional argument that sets a sensitivity threshold (`pct.1` column representing the percent of on-cluster cells expressing the marker) on both query and db Marker dataframes. The default is set to 0.25. This can be interpreted as "at least 25 percent of a cluster must contain at least 1 count of a gene for it to be considered a marker."}
\item{spec_thresh}{Optional argument that sets a specificity threshold (`pct.2` column representing the percent of off-cluster cells expressing the marker) on both query and db Marker dataframes. The default is set to 0.75. This can be interpreted as "no more than75 percent of a cells not in the cluster of interest can contain >= 1 count of a gene for it to be considered a marker."}
}
\examples{
#Intialize with default thresholds
ann=Initialize_Markers(ann)

#Initialize with new thresholds
ann=Initialize_Markers(ann, p_thresh = 0.01, FC_thresh = 1.75, sens_thresh = 0.30, spec_thresh = 0.80)
}