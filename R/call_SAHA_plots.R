#' Retrieves plot data from a SAHA analysis object
#'
#' This function returns plot data based on the specified plot type and data type.
#' If `data_type` is not provided, it attempts to determine the data type from the `ann` object.
#'
#' Supported `plot_type` values include:
#' - "self-similarity": Retrieves self-similarity heatmaps based on the provided or inferred data type.
#' - "Marker-based": Returns the dot plot summarizing marker-based annotation results.
#' - "Marker-free": Returns the heatmap summarizing marker-free (average expression) annotation results.
#'
#' @param ann A SAHA analysis object containing marker-based and/or average expression-based results.
#' 
#' @param plot_type A string specifying the type of plot to retrieve. Options include "self-similarity", "Marker-based", or "Marker-free".
#' 
#' @param data_type A string specifying the data modality ("Markers", "AvgExp", or "Both") to retrieve plots for. If not provided, the function attempts to infer it from `ann@data_type`.
#' 
#' @return A ggplot2-compatible object (e.g., heatmap or dot plot) ready for rendering, or prints a message if data is missing or unspecified.
#'
#' @export

call_SAHA_plots <- function(ann, plot_type, data_type=NULL){
   if (is.null(data_type)) {
      paste0("SAHA detected ", ann@data_type, " information loaded. If the appropriate results do not display please ensure the pipeline has been run & results have been saved in ann@results.")
      if (ann@data_type == "Markers" & plot_type == "self-similarity") {
         ann@results$self_similarity$similiarity_heatmap_markers
      }else if (ann@data_type == "AvgExp" & plot_type == "self-similarity"){
         ann@results$self_similarity$similiarity_heatmap_avgexp
      }else if (ann@data_type == "Markers" & plot_type == "Marker-based"){
         ann@results$marker_based$dotplot_all
      }else if (ann@data_type == "Markers" & plot_type == "Marker-free"){
         ann@results$marker_free$heatmap
      }else { print("Multiple data types detected, please specify data_type you would like to visualize.")}
   }else{
      paste0("You have selected to visualize ", data_type,".")
      if (data_type == "Markers" & plot_type == "self-similarity") {
         ann@results$self_similarity$similiarity_heatmap_markers
      }else if (data_type == "AvgExp" & plot_type == "self-similarity"){
         ann@results$self_similarity$similiarity_heatmap_avgexp
      }else if (data_type == "Markers" & plot_type == "Marker-based"){
         ann@results$marker_based$dotplot_all
      }else if (data_type == "AvgExp" & plot_type == "Marker-free"){
         ann@results$marker_free$heatmap
      }
   }
}
