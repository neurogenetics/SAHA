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

SimilarityDend <- function(ann,sim_type,annotation){
   if (sim_type == "Markers") {
      ht=ann@results$self_similarity$similiarity_heatmap_markers
   }else if (sim_type == "AvgExp"){
      ht=ann@results$self_similarity$similiarity_heatmap_avgexp
   }
   dend=suppressWarnings(row_dend(ht))

   # Assuming you have a dendrogram object named 'dend' and a dataframe named 'mapping' with columns 'old_names' and 'new_names'

   # Create a named vector for label mapping
   label_map <- setNames(annotation$best_match, annotation$cluster)

   # Rename labels in the dendrogram
   labels(dend) <- sapply(labels(dend), function(x) {
      ifelse(is.na(label_map[x]), x, label_map[x])
   })

   # Save previous graphical parameters
   old_par <- par(no.readonly = TRUE)
   # Set specific parameters for the dendrogram plot
   par(mar = c(8, 4, 1, 2) + 0.1)  # Adjust margins as needed
   return(plot(dend))
   par(old_par)

}
