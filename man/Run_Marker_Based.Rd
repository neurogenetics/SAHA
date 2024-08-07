\name{Run_Marker_Based}
\alias{Run_Marker_Based}
\title{Run: Marker Based Annotation of Cell Types}
\usage{
Run_Marker_Based(ann)
}
\description{
Run_Marker_Based requires a SAHA object (ann) to create a `master dataframe` that will be stored in saha_object@ann2. To run this step, the user will have needed to created a saha object using Create_SAHA_object().  
}
\arguments{
\item{ann}{Required argument that takes an initialized SAHA object that uses cell markers for each of the user's clusters and runs a Fisher's Exact Test against a marker database. To use the Allen Brain Cell Atlas, the db slot must contain markers obtained from the ABC atlas. Please see Create_SAHA_Object() for more information. Additionally, the master dataframe will be used for all downstream marker based analysis. The current build is capable of loading custom databases however there may be some errors with Visualization steps in downstream analysis.}
}
\examples{
ann=Run_Marker_Based(ann)
}