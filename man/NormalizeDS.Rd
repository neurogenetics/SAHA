\name{NormalizeDS}
\alias{NormalizeDS}
\title{Marker-free annotation: Normalization}
\description{Scales gene expression data across cell types in the query and database for marker-free annotation.}
\usage{NormalizeDS(ann, 
assay_db,
assay_query)}
\arguments{
  \item{ann}{A SAHA object which has been downsampled using Downsample() for marker-free annotation.}
  \item{assay_db}{Assay type of the database. For scRNA-seq, usually "RNA". For spatial, can be "Xenium", "NanoString", etc. If using default Allen Brain Cell Atlas database, defaults to "RNA".}
  \item{assay_query}{Assay type of the query dataset. For scRNA-seq, usually "RNA". For spatial, can be "Xenium", "NanoString", etc.}
}

\examples{
NormalizeDS(ann, assay_db = "RNA", assay_query = "RNA")
NormalizeDS(ann, assay_db = "RNA", assay_query = "Xenium")
}