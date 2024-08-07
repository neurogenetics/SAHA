\name{Marker_Diversity}
\alias{Marker_Diversity}
\title{Marker Diveristy: Shannon Diversity Index Per Cluster}
\usage{
SAHA(ann)
}
\description{
BETA FUNCTION: Creates a scatterplot to visualize the diversity of markers in the ann@ann1 (initialized and possibly tuned markers from ann@query$Markers) per cluster. Low diveristy clusters indicate clusters that may be hard toannotate, while high diversity clusters indicate clusters that may map to multiple annotations.

Not intended to be written into the SAHA object and must be re-run to create the plot.
}
\arguments{
\item{ann}{A SAHA object.}
}
\examples{
Marker_Diversity(ann)
}
