#' Creates a heatmap visualizing gene expression correlations
#'
#' Generates a heatmap representing correlation coefficients between
#' query and database expression profiles. Includes optional metadata annotation.
#'
#' @param ann An object containing correlation data (`ann@results$marker_free$corr`).
#' @param meta A data frame containing cell type metadata (optional).
#' @param ABC Character string indicating the brain region to filter by in the Allen Brain Atlas data (optional).
#'
#' @return The `ann` object with the heatmap stored in `ann@results$marker_free$heatmap`.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation ht_opt
#' @importFrom circlize colorRamp2
#' @importFrom dplyr %>% distinct na.omit
#'
#' @export

Create_MarkerFree_Viz <- function(ann,
                                  meta=NULL#,
   #                               ABC=NULL #deprecated option for specifying ABC region
                                  ){

   # set global heatmap options
   ht_opt(RESET = TRUE)
   ht_opt(legend_border = "black",
          heatmap_border = TRUE)


   # set colors for heatmap

   col_fun = colorRamp2(c(-1, 0, 1), c("#05409e", "white", "#e64a02"))

   ann3=ann@results$marker_free$corr

   # make data matrix
   rownames(ann3) <- gsub("^query\\.", "", rownames(ann3))
   colnames(ann3) <- gsub("^db\\.", "", colnames(ann3))

   mat <- data.matrix(ann3)
   if (!is.null(meta)) {
      # heatmap annotation for levels

      annotation_data <- meta[meta$subclass_per %in% colnames(mat), ]
      # if (!is.null(ABC)) { #deprecated option for specifying region of ABC metadata
      #    if (ABC=="ISOCTX") {#pre-prepared ABC data uses 10Xv3 data for ISOCTX
      #       annotation_data <- annotation_data[annotation_data$library_method %in% "10Xv3", ]
      #    }else{#pre-prepared ABC data uses 10Xv2 data for all other regions
      #       annotation_data <- annotation_data[annotation_data$library_method %in% "10Xv2", ]
      #    }
      # }

      annotation_data <- na.omit(annotation_data[match(colnames(mat), annotation_data$subclass_per), ])



      # set colors for annotation based on ABC-assigned colors

      neurotransmitter_unique <- annotation_data %>% distinct(neurotransmitter, .keep_all = TRUE)
      neurotransmitter_colors <- setNames(neurotransmitter_unique$neurotransmitter_color, nm = neurotransmitter_unique$neurotransmitter)

      class_unique <- annotation_data %>% distinct(class, .keep_all = TRUE)
      class_colors <- setNames(class_unique$class_color, class_unique$class)



      # create annotation

      ha <- HeatmapAnnotation(
         Neurotransmitter = annotation_data$neurotransmitter,
         Class = annotation_data$class,
         col = list(
            Class = class_colors,
            Neurotransmitter = neurotransmitter_colors
         )
      )


      ht1 <- Heatmap(mat,
                     name = "Pearson R",
                     col = col_fun,
                     cluster_rows = T,
                     cluster_columns = T,
                     row_names_side = "left",
                     row_names_gp = gpar(fontsize = 11),
                     column_names_side = "bottom",
                     column_names_gp = gpar(fontsize = 11),
                     row_title = "cluster",
                     row_title_rot = 90,
                     row_title_side = "left",
                     row_title_gp = gpar(fontface = "bold"),
                     #column_title = "ABC atlas",
                     column_title_rot = 0,
                     column_title_side = "bottom",
                     column_title_gp = gpar(fontface = "bold"),
                     top_annotation = ha
                     )
   }else{
      ht1 <- Heatmap(mat,
                     name = "Pearson correlation",
                     col = col_fun,
                     cluster_rows = F,
                     cluster_columns = T,
                     row_names_side = "left",
                     row_names_gp = gpar(fontsize = 11),
                     column_names_side = "bottom",
                     column_names_gp = gpar(fontsize = 11),
                     row_title = "query cluster",
                     row_title_rot = 90,
                     row_title_side = "left",
                     row_title_gp = gpar(fontface = "bold"),
                     column_title = "db",
                     column_title_rot = 0,
                     column_title_side = "bottom",
                     column_title_gp = gpar(fontface = "bold"))

   }



   ann@results$marker_free$heatmap=ht1
   return(ann)
}
