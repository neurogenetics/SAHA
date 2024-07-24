\name{SAHA}
\alias{SAHA}
\title{Semi-Automated Hand Annotation: Analyze}
\usage{
SAHA(query,
   database = "Panglao",
   nameoforgan = NULL,
   mode = "db",
   speciesquery = NULL,
   speciesdatabase = NULL)
}
\description{
Runs one of three flavors of Semi-Automated Hand Annotation on a query dataset against a database of interest.
}
\arguments{
\item{query}{A dataframe containing marker genes or summary expression profiles. This argument can take either a loaded data.frame or the path to a dataframe, assuming the columns X, Y, Z are present.}
\item{database}{A dataframe containing marker genes or summary expression profiles from a given database to be tested. By default, a PangloaDB instance from YYYY-MM-DD will be loaded from the package data. When selecting a custom databse, this argument can take either a loaded data.frame or the path to a dataframe, assuming the columns X, Y, Z are present. }
\item{nameoforgan}{Optional argument that can specify to only test a subset of the Pangloa database. IE, the argument "brain" will subset only celltypes reported in the brain.}
\item{arg}{...}
}
\examples{
ann <- SAHA()
}
