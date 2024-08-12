#' Create Self-Similarity Visualization
#'
#' This function creates a self-similarity visualization based on the provided
#' annotation object. It generates heatmaps for marker self-similarity or average
#' expression self-similarity, depending on the specified slot
#'
#' @param ann An object containing annotation data. It should have the required
#'             matrices (`Marker_selfsim_matrix` or `AvgExp_selfsim_matrix`) and
#'             result storage.
#' @param slot A character string indicating which slot to use. Should be either
#'             "Markers" or "AvgExp".
#' @param assay_db A character string specifying the assay database prefix. Default is "RNA".
#' @return The updated annotation object with self-similarity heatmaps stored in `ann@results$self_similarity`.
#' @importFrom ComplexHeatmap Heatmap rowAnnotation columnAnnotation
#' @importFrom stats cor
#' @importFrom grDevices colorRampPalette
#' @export
Create_SelfSimilarity_Viz <- function(ann, slot,assay_db="RNA") {
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
      ann@results$self_similarity$similiarity_heatmap_markers = p1
      return(ann)
   }else if(slot == "AvgExp"){
      cor_mat=ann@ann1$AvgExp_selfsim_matrix
      cor_mat[is.na(cor_mat)]<-0
      rownames(cor_mat) <- gsub(paste0("^",assay_db,"\\.g"), "", rownames(cor_mat))
      colnames(cor_mat) <- gsub(paste0("^",assay_db,"\\.g"), "", colnames(cor_mat))
      col <- colorRampPalette(c("white", "red"))(200)
      # p2 <- corrplot(cor_mat, method = "color", type = "upper", order = "hclust", tl.col = "black",col=col, tl.srt = 45)
      p2 <- Heatmap(cor_mat,
                    name = "Self Similarity Matrix",  # Optional: Set a name for the heatmap
                    col = col,  # Color palette
                    cluster_rows = TRUE,  # Hierarchical clustering for rows
                    cluster_columns = TRUE  # Hierarchical clustering for columns
                    #show_colorbar = TRUE
      )  # Display the colorbar

      ann@results$self_similarity$similiarity_heatmap_avgexp = p2
      return(ann)
   }


}
