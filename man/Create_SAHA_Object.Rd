\name{Create_SAHA_Object}
\alias{Create_SAHA_Object}
\title{Create Semi Automated Hand Annotated Object}
\usage{
Create_SAHA_Object(query,
   db,
   data_type,
   existing=NULL)
}
\description{
Creating a SAHA object allows the user to consolidate all marker and average expression data in one place. By creating the object, this function will also create empty slots for all possible SAHA analysis.
}
\arguments{
\item{query}{A dataframe containing marker genes or summary expression profiles to be tested. This argument can take either a loaded data.frame or the path to a dataframe, assuming the format is:
* Markers: rows = gene*cluster summary statistics; columns = `avg_log2FC`, `cluster`, `pct.1`, and `pct.2`, and `p_val_adj`
* AvgExp: rows = gene names ; columns = {Assay_Name}.{g}.{Cluster_Name}}
\item{db}{A dataframe containing marker genes or summary expression profiles from a given database to be queried against. This argument can take either a loaded data.frame or the path to a dataframe, assuming the format is:
* Markers: rows = gene*cluster summary statistics; columns = `avg_log2FC`, `cluster`, `pct.1`, and `pct.2`, and `p_val_adj`
* AvgExp: rows = gene names ; columns = {Assay_Name}.{g}.{Cluster_Name} }
\item{data_type}{Required argument that can specify whether "Markers" or "AvgExp" data is being loaded. Only one may be loaded at a time.}
\item{existing}{Argument that takes an existing SAHA object for loading whichever data_type was loaded first. Intended to also be written to this SAHA object to load both modalities. See examples for usage.}
}
\examples{
#Loading only Markers data (limits capability to Marker-based only)
ann=Create_SAHA_object(query = My_markers,db = ABC_ISOCTX_Markers,data_type = "Markers")

#Loading Markers data, followed by AvgExp data
ann=Create_SAHA_object(query = My_markers,db = ABC_ISOCTX_Markers,data_type = "Markers")
ann=Create_SAHA_object(query = My_AvgExp,db = ABC_ISOCTX_AE,data_type = "AvgExp",existing = ann)

#Loading AvgExp data, followed by Marker data
ann=Create_SAHA_object(query = My_AvgExp,db = ABC_ISOCTX_AE,data_type = "AvgExp")
ann=Create_SAHA_object(query = My_markers,db = ABC_ISOCTX_Markers,data_type = "Markers",existing = ann)


}
