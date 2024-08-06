\name{Investigate_MarkerFree}
\alias{Investigate_MarkerFree}
\title{Marker-free annotation: Investigating specific clusters}
\description{Generates a scatter plot of scaled gene expression for a specified query cluster and a database cell type.}
\usage{Investigate_MarkerFree(ann,
query_cluster,
db_cell_type)
}
\arguments{
  \item{ann}{SAHA object which has been passed through Create_MarkerFree_Viz() for marker-free identification.}
  \item{query_cluster}{Number of unknown query cluster to investigate.}
  \item{db_cell_type}{Name of known database cell type to investigate.}
}

\examples{
Investigate_MarkerFree(ann, query_cluster = 4, db_cell_type = "037.DG.Glut")
}