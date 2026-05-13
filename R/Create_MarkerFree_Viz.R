#' Create Marker-free Visualization
#' 
#' Generates a heatmap representing correlation coefficients. 
#' Includes options for clustering and tiered significance markers (***, **, *).
#'
#' @param ann An object containing correlation data (`ann@results$marker_free$corr`) and p-values.
#' @param meta A data frame containing cell type metadata (optional).
#' @param cluster_rows Logical, whether to cluster query clusters (rows). Default is TRUE.
#' @param cluster_columns Logical, whether to cluster database cell types (columns). Default is TRUE.
#' @param add_significance Logical, whether to add tiered '*' to significant correlations. Default is FALSE.
#'
#' @return The `ann` object with the heatmap stored in `ann@results$marker_free$heatmap`.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation ht_opt
#' @importFrom grid gpar grid.text
#' @importFrom circlize colorRamp2
#' @importFrom dplyr %>% distinct
#' @importFrom stats na.omit
#'
#' @export

Create_MarkerFree_Viz <- function(ann,
                                  meta = NULL,
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  add_significance = FALSE){
   
   # ROBUST CLEANING: Standardizes names and removes noise words
   clean_string <- function(x) {
      x <- as.character(x)
      x <- gsub("[._/]", " ", x)
      x <- gsub("\\b(query|RNA|SCT|integrated|db)\\b", "", x, ignore.case = TRUE, perl = TRUE)
      x <- trimws(gsub("\\s+", " ", x))
      return(x)
   }
   
   # Set global heatmap options
   ComplexHeatmap::ht_opt(RESET = TRUE)
   ComplexHeatmap::ht_opt(legend_border = "black", heatmap_border = TRUE)
   
   # Set colors for heatmap
   col_fun = circlize::colorRamp2(c(-1, 0, 1), c("#05409e", "white", "#e64a02"))
   
   # Pull and clean the correlation matrix
   corr_df <- ann@results$marker_free$corr
   rownames(corr_df) <- clean_string(rownames(corr_df))
   colnames(corr_df) <- clean_string(colnames(corr_df))
   mat <- data.matrix(corr_df)
   
   # Define Tiered Significance Layer
   layer_logic <- NULL
   if (add_significance) {
      if (is.null(ann@results$marker_free$pvalue)) {
         stop("P-value matrix not found in ann@results$marker_free$pvalue. Please rerun CorrelateDS().")
      }
      
      # Prepare p-value matrix and ensure alignment with the correlation matrix
      p_df <- ann@results$marker_free$pvalue
      rownames(p_df) <- clean_string(rownames(p_df))
      colnames(p_df) <- clean_string(colnames(p_df))
      p_mat <- data.matrix(p_df[rownames(corr_df), colnames(corr_df)])
      
      # Vectorized layer function for tiered stars
      layer_logic <- function(j, i, x, y, w, h, fill) {
         p_values <- p_mat[cbind(i, j)]
         
         # Logic:
         # < 0.001: ***
         # < 0.01:  **
         # < 0.05:  *
         # >= 0.05: (nothing)
         stars <- cut(p_values, 
                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                      labels = c("﹡﹡﹡", "﹡﹡﹡", "﹡", ""),
                      right = FALSE) # intervals are [a, b)
         
         stars <- as.character(stars)
         stars[is.na(stars)] <- ""
         is_visible <- stars != ""
         
         if (any(is_visible)) {
            grid::grid.text(stars[is_visible], x[is_visible], y[is_visible], 
                            gp = grid::gpar(col = "black", fontsize = 12, fontface = "bold"))
         }
      }
   }
   
   if (!is.null(meta)) {
      meta$subclass_per <- clean_string(meta$subclass_per)
      annotation_data <- meta[meta$subclass_per %in% colnames(mat), ]
      annotation_data <- stats::na.omit(annotation_data[match(colnames(mat), annotation_data$subclass_per), ])
      
      neurotransmitter_unique <- annotation_data %>% dplyr::distinct(neurotransmitter, .keep_all = TRUE)
      neurotransmitter_colors <- setNames(neurotransmitter_unique$neurotransmitter_color, nm = neurotransmitter_unique$neurotransmitter)
      
      class_unique <- annotation_data %>% dplyr::distinct(class, .keep_all = TRUE)
      class_colors <- setNames(class_unique$class_color, class_unique$class)
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
         Neurotransmitter = annotation_data$neurotransmitter,
         Class = annotation_data$class,
         col = list(Class = class_colors, Neurotransmitter = neurotransmitter_colors)
      )
      
      ht1 <- ComplexHeatmap::Heatmap(mat,
                                     name = "R",
                                     col = col_fun,
                                     cluster_rows = cluster_rows,
                                     cluster_columns = cluster_columns,
                                     layer_fun = layer_logic,
                                     row_names_side = "left",
                                     row_names_gp = grid::gpar(fontsize = 11),
                                     column_names_side = "bottom",
                                     column_names_gp = grid::gpar(fontsize = 11),
                                     row_title = "query cluster",
                                     row_title_rot = 90,
                                     row_title_side = "left",
                                     row_title_gp = grid::gpar(fontface = "bold"),
                                     column_title_rot = 0,
                                     column_title_side = "bottom",
                                     column_title_gp = grid::gpar(fontface = "bold"),
                                     top_annotation = ha
      )
   } else {
      ht1 <- ComplexHeatmap::Heatmap(mat,
                                     name = "correlation",
                                     col = col_fun,
                                     cluster_rows = cluster_rows,
                                     cluster_columns = cluster_columns,
                                     layer_fun = layer_logic,
                                     row_names_side = "left",
                                     row_names_gp = grid::gpar(fontsize = 11),
                                     column_names_side = "bottom",
                                     column_names_gp = grid::gpar(fontsize = 11),
                                     row_title = "query cluster",
                                     row_title_rot = 90,
                                     row_title_side = "left",
                                     row_title_gp = grid::gpar(fontface = "bold"),
                                     column_title = "db",
                                     column_title_rot = 0,
                                     column_title_side = "bottom",
                                     column_title_gp = grid::gpar(fontface = "bold"))
   }
   
   ann@results$marker_free$heatmap = ht1
   return(ann)
}
