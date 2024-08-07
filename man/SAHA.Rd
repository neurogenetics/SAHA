\name{SAHA}
\alias{SAHA}
\title{Semi-Automated Hand Annotation: Quickstart}
\usage{
SAHA(query,
   db,
   meta,
   data_type)
}
\description{
Quickstart will automatically run Marker-based or Marker-free pipelines from input of either marker- or avgexp- dataframes to the summarized plots.
}
\arguments{
\item{query}{A dataframe containing marker genes or summary expression profiles to be tested. This argument can take either a loaded data.frame or the path to a dataframe, assuming the format is:
* Markers: rows = gene*cluster summary statistics; columns = `avg_log2FC`, `cluster`, `pct.1`, and `pct.2`, and `p_val_adj`
* AvgExp: rows = gene names ; columns = {Assay_Name}.{g}.{Cluster_Name}}
\item{db}{A dataframe containing marker genes or summary expression profiles from a given database to be queried against. This argument can take either a loaded data.frame or the path to a dataframe, assuming the format is:
* Markers: rows = gene*cluster summary statistics; columns = `avg_log2FC`, `cluster`, `pct.1`, and `pct.2`, and `p_val_adj`
* AvgExp: rows = gene names ; columns = {Assay_Name}.{g}.{Cluster_Name} }
\item{meta}{Optional argument to load in dataframe providing additional information about the cell types. This metadata should contain information based on (or similar to) the ABC metadata: class, class color, neurotransmitter, and neurotransmitter color information for each cell type.}
\item{data_type}{Required argument that can specify whether "Markers" or "AvgExp" data is being loaded. Only one may be loaded at a time.}
}
\examples{
#Quickstart to Marker-based dotplot
SAHA(query = My_markers, db = ABC_ISOCTX_Markers,meta = meta, data_type = "Markers")
#Quickstart to Marker-free heatmap
SAHA(query = My_AE, db = ABC_ISOCTX_AE,meta = meta, data_type = "AvgExp")
}
