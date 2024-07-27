# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
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

Create_SelfSimilarity_Viz <- function(ann, slot) {
   # Ensure data is numeric
   if(slot=="Markers"){
      df <- as.matrix(ann@ann1$Marker_selfsim_matrix)
      df[is.na(df)] <- 0  # Replace NA with 0

      # Calculate shared features
      overlap_matrix <- df %*% t(df)

      # Convert to matrix
      mat <- as.matrix(overlap_matrix)

      # Create row and column annotations (optional)
      row_ann <- rowAnnotation(cluster = factor(rownames(mat)))
      col_ann <- columnAnnotation(cluster = factor(colnames(mat)))

      mat2 <- log10(mat+1)

      # Create the heatmap
      p1 <- Heatmap(mat2,
                    name = paste0("Shared Markers \n Log10(Overlap + 1)"),
                    column_title = "Cluster",
                    row_title = "Cluster",
                    top_annotation = NULL,
                    left_annotation = NULL,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
      ann@results$self_similarity = list(similiarity_heatmap = p1)
      return(ann)
   }else if(slot == "AvgExp"){
      cor_mat=ann@ann1$AvgExp_selfsim_matrix
      cor_mat[is.na(cor_mat)]<-0
      col <- colorRampPalette(c("blue", "white", "red"))(200)
      p2 <- corrplot(cor_mat, method = "color", type = "upper", order = "hclust", tl.col = "black",col=col, tl.srt = 45)
      ann@results$self_similarity = list(similiarity_heatmap = p2)
      return(ann)
   }


}
