\name{NormalizeDS}
\alias{NormalizeDS}
\title{Marker-free annotation: Normalization}
\description{Scales gene expression data across cell types in the query and database for marker-free annotation.}
\usage{
NormalizeDS(ann,
            assay_db = "RNA",
            assay_query = assay_db,
            norm_method = "across_clust")
}
\arguments{
  \item{ann}{A SAHA object which has been downsampled using Downsample() for marker-free annotation.}
  \item{assay_db}{Assay type of the database. For scRNA-seq, usually "RNA". For spatial, can be "Xenium", "NanoString", etc. Defaults to "RNA".}
  \item{assay_query}{Assay type of the query dataset. For scRNA-seq, usually "RNA". For spatial, can be "Xenium", "NanoString", etc. Defaults to the value of \code{assay_db}.}
  \item{norm_method}{Normalization method to use. Options are "across_clust" (default) or "within_clust". We recommend normalizing expression across clusters for highly heterogenous databases (different biological cell types) and within cluster for more homogeneous databases (i.e. sub celltypes or cell states).}
}

\examples{
NormalizeDS(ann, assay_db = "RNA", assay_query = "RNA")
NormalizeDS(ann, assay_db = "RNA", assay_query = "Xenium")
NormalizeDS(ann, assay_db = "RNA", assay_query = "RNA", norm_method = "within_clust")
}
