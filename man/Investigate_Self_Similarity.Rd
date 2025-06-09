\name{Investigate_Self_Similarity}
\alias{Investigate_Self_Similarity}
\title{Investigate: Compare Two Clusters Based on Shared Characteristics}
\usage{
Investigate_Self_Similarity(ann,
   cluster1,
   cluster2,
   shared_df=NULL)
}
\description{
After the calculation of self-similarity based on markers for a query dataset, Investigate_Self_Similarity() allows for the creation of either a venn diagram between two clusters or a joint dataframe of the markers of two clusters.
}
\arguments{
\item{cluster1}{Name of one cluster to be compared, name must match the exact name of the clusters in query dataset.}
\item{db}{Name of the second cluster to be compared, name must match the exact name of the clusters in query dataset. }
\item{shared_df}{By default this function returns a venn diagram, if shared_df=TRUE the function returns a joint dataframe of the markers specific to and shared between cluster1 and cluster2.}
}
\examples{
#Comparison between markers from cluster "1" and cluster "10".
Investigate_Self_Similarity(ann,cluster1 = 1, cluster2 = 10)
#Dataframe of markers making up "cluster "1" and cluster "10".
shared_marker_df=Investigate_Self_Similarity(ann,cluster1 = 1, cluster2 = 10,shared_df = TRUE)
}
