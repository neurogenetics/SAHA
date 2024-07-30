# SAHA
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#################################################
#Rework so this is the quickstart option########
#################################################

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
