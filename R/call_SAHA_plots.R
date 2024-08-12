#' Retrieves plot data from a SAHA analysis object
#'
#' This function returns plot data based on the specified plot type and data type.
#' If `data_type` is not provided, it attempts to determine the data type from the `ann` object.
#'
#' @param ann A SAHA analysis object containing results.
#' @param plot_type A character string specifying the desired plot type (e.g., "self-similarity", "Marker-based").
#' @param data_type An optional character string specifying the data type to visualize (e.g., "Markers", "AvgExp").
#' @return A data frame or matrix containing the plot data.
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
