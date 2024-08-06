\name{Create_MarkerFree_Viz}
\alias{Create_MarkerFree_Viz}
\title{Marker-free annotation: Creating marker-free annotation heatmap}
\description{Generates a heatmap of correlation of gene expression profiles between query dataset and database.}
\usage{Create_MarkerFree_Viz(ann,
facet,
meta,
ABC,
chemistry)
}
\arguments{
  \item{ann}{A SAHA object which has been correlated using CorrelateDS() for marker-free annotation.}
  \item{facet}{Logical, specifies whether to facet in additional metadata from the database about cell types.}
  \item{meta}{A dataframe containing metadata about database cell types. Only necessary if facet = TRUE. Dataframe requires several specific columns:
  "subclass_per": these are the names of the cell types from the database which you are testing against. Words must be separated by periods rather than spaces or underscores.
  "neurotransmitter": this is the highest-level (broadest) classification for each cell type identified in "subclass_per".
  "class": this is a more specific level (between subclass and neurotransmitter) of classification for each cell type identified in "subclass_per".
  "neurotransmitter_color": hex code of color that should be used for the "neurotransmitter" classification.
  "class_color": hex code of color that should be used for the "class" classification.}
  \item{ABC}{Logical, specifies if database is the Allen Brain Cell Atlas. Default = NULL}
  \item{chemistry}{Chemistry for ABC database. Options are "10Xv2", when database is "isocortex", and "10Xv3", for everything else.}
}

\examples{
Create_MarkerFree_Viz(ann, facet = TRUE, meta = meta, ABC = TRUE, chemistry = "10Xv2")
}