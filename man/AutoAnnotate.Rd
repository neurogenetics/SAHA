\name{AutoAnnotate}
\alias{AutoAnnotate}
\title{Annotate: Automatic Annotation By Best Scores}
\usage{
AutoAnnotate(ann,
             data_type = NULL)
}
\description{
AutoAnnotate requires a SAHA object (ann) and data_type as inputs to function. This tool is used to view cluster identification quickly to aid in annotation.
The user is required to pick a datatype of markers, average experssion or both.
If the user chooses "Markers":
The annotation will first find the most significant annotation within each cluster and assign all of these into a new dataframe. The output dataframe's structure is: Cluster number, best cell type annotation, proportion, and sinificance.
If the user chooses "AvgExp":
The annotation will use the highest correlation in terms of average expression to detect the best match. The output dataframe's structure is: Cluster number, best cell type annotation, and the percent correlation.
If the user chooses "Both":
The best annotation is decided by create a dataframe that compares both the Markers and AvgExp dataframe's respective best_matches to one another, forming a consensus. If this tool is unable to run or find a best match for either Markers or AvgExp, it labels the consensus as "INCONCLUSIVE". If the two matches disagree with one another, or one of the matches are inconclusive, the consensus is labelled as "DISAGREEMENT".
To further refine best matches, use SemiAutoAnnotation() (check the description).
}
\arguments{
\item{ann}{Required argument that takes a SAHA object generated from Run_Marker_Based().}
\item{datatype}{Required argument that specifies which type of data the user is querying. The three available options are 'Markers' ,'AvgExp', and 'Both.'
}
\examples{
auto=AutoAnnotate(ann,data_type = "Markers")
auto=AutoAnnotate(ann,data_type = "AvgExp")
auto=AutoAnnotate(ann,data_type = "Both")
}
