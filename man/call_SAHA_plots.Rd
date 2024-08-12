\name{call_SAHA_plots}
\alias{call_SAHA_plots}
\title{Plotting SAHA Outputs}
\usage{
call_SAHA_plots(ann,
   plot_type,
   data_type=NULL)
}
\description{
The static outputs of SAHA (Self-Similarity, Marker-Based, and Marker-Free) are stored in ann@results. To avoid needing to call these plots manually, the call_SAHA_plots() function allows simple plotting of the stored data.
}
\arguments{
\item{ann}{A SAHA object.}
\item{plot_type}{Plots the main output of one of three SAHA flavors:
   * "self-similarity": Heatmap comparing the pairwise comparison of the loaded query clusters based either on Markers (heatmap containing counts of shared markers) or AvgExp (heatmap of R^2 value based on the expression of shared genes).
   * "Marker-based": Dotplot representing the proportion of db annotation markers present in each query cluster.
   * "Marker-free": Heatmap representing the R^2 (Pearson's) of the Min-max scaled counts between each db annotation and query cluster given the shared genes.
}
\item{data_type}{Optional argument that if left NULL will plot whatever analysis has been run. If more than one flavor of SAHA is loaded into the object, the user will need to specific which data_type the plot will be based on:
   * "Markers": Marker dataframe formated as the default output of Seurat::FindAllMarkers()
   * "AvgExp": AvgExp dataframe containing the average expression of each gene (rows) per cluster (column).
}}
\examples{
#Call self-similarity heatmap based on markers
call_SAHA_plots(ann, plot_type = "self-similarity",data_type = "Markers")

#Call Marker-based dotplot
call_SAHA_plots(ann, plot_type = "Marker-based",data_type = "Markers")

#Call Marker-free heatmap
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")

}
