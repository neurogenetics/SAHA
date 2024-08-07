\name{SAHA_lookup_cluster}
\alias{SAHA_lookup_cluster}
\title{Semi-Automated Hand Annotation: SAHA_lookup_cluster}
\usage{
SAHA_lookup_cluster(cluster,
	ann)
}
\description{
SAHA_lookup_cluster is a nested function that looks through each cluster of the user's SAHA object data to structure the master dataframe. This function will return a summary of how many times a cluster appears in the database based on associated genes. Then, it rowbinds information into the master dataframe used for marker-based analysis. 
No user interaction required, change at your own risk.
}
\arguments{
\item{cluster}{Argument used to loop through the queried dataset.}
\item{ann}{Required argument that takes an initialized SAHA object(ann).}
}
\examples{
}