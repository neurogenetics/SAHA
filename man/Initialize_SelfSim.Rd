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
\item{query}{A dataframe containing marker genes or summary expression profiles. This argument can take either a loaded data.frame or the path to a dataframe, assuming the columns X, Y, Z are present.}
\item{db}{A dataframe containing marker genes or summary expression profiles from a given database to be tested. By default, a PangloaDB instance from YYYY-MM-DD will be loaded from the package data. When selecting a custom databse, this argument can take either a loaded data.frame or the path to a dataframe, assuming the columns X, Y, Z are present. }
\item{meta}{Optional argument that can specify to only test a subset of the Pangloa database. IE, the argument "brain" will subset only celltypes reported in the brain.}
\item{data_type}{Optional argument that can specify to only test a subset of the Pangloa database. IE, the argument "brain" will subset only celltypes reported in the brain.}
}
\examples{
SAHA(query = query_markers, db = db_markers , meta = meta, data_type = "Markers")
SAHA(query = query_avgexp, db = db_avgexp , meta = meta, data_type = "AvgExp")
}
