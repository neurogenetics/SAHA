#' Creates marker-based visualizations for cell type annotation
#'
#' Generates various visualizations based on marker-based cell type annotation results
#' stored in the `ann` object. Requires cell type metadata in the `meta`
#' data frame for adding cell type class information as facets.
#'
#' @param ann An SAHA analysis object containing annotation results (likely from `Annotate()` or similar).
#' @param meta A data frame containing cell type metadata.  Must include columns:
#'   `subclass_spa` (cell type names matching those in `ann`), `class` (broader cell type classification),
#'   `neurotransmitter` (even broader classification), `class_color` (hex code for class color),
#'   and `neurotransmitter_color` (hex code for neurotransmitter color). If NULL, faceting is disabled.
#'
#' @return The `ann` object, with the following plots added to the `results$marker_based` slot:
#'   \itemize{
#'     \item `markers_barplot`: A bar plot showing the total number of markers per cell type.
#'     \item `dotplot_all`: A dot plot showing all markers for each cell type and cluster.
#'     \item `dotplot_sig`: A dot plot showing only significant markers for each cell type and cluster.
#'     \item `dotplot_best`: A dot plot showing the most significant marker for each cell type and cluster.
#'   }
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% group_by, filter, slice_min, ungroup, mutate, select, as.data.frame
#' @importFrom stringr str_remove_all, str_replace, str_replace_all, if_else
#' @importFrom data.table full_join
#' @importFrom ggplot2 ggplot, geom_point, scale_size_continuous, labs, scale_color_manual, theme_bw, theme, element_text, facet_grid, geom_bar, geom_text, coord_cartesian
#' @importFrom ggplot2 strip_themed, elem_list_rect, guide_legend, guides # Added for facet colors
#' @importFrom ggh4x facet_grid2, strip_themed
#' @importFrom forcats fct_drop, fct_relevel
#'
#' @export


Create_MarkerBased_Viz <- function (ann, meta = NULL){
   # Check if meta is a data frame or NULL
   if (!is.null(meta)) {
      if (inherits(meta, "data.frame")) {
         # Meta is already a data frame, use it directly
         meta <- meta
      } else {
         # Meta is a path to a CSV file, read it into a data frame
         meta <- read.csv(meta)
      }
   } else {
      # Handle case where meta is NULL
      # warning("meta argument is NULL. Function might not work as expected.")
      # meta <- data.frame()  # Create an empty data frame
   }
   master_df = ann@ann2
   master_df$celltype <- gsub("-", " ", master_df$celltype)
   temp_met=meta[meta$subclass_spa %in% unique(master_df$celltype),]

   ###########ALL marker names need to have all spaces and NOT - !! - then this will work...

   master_df$class <- sapply(master_df$celltype, function(x) {
      match_idx <- match(x, temp_met$subclass_spa)
      if (length(match_idx) > 0) temp_met$class[match_idx] else NA
   })

   master_df$class_color <- sapply(master_df$celltype, function(x) {
      match_idx <- match(x, temp_met$subclass_spa)
      if (length(match_idx) > 0) temp_met$class_color[match_idx] else NA
   })

   master_df$nt <- sapply(master_df$celltype, function(x) {
      match_idx <- match(x, temp_met$subclass_spa)
      if (length(match_idx) > 0) temp_met$neurotransmitter[match_idx] else NA
   })

   master_df$nt_color <- sapply(master_df$celltype, function(x) {
      match_idx <- match(x, temp_met$subclass_spa)
      if (length(match_idx) > 0) temp_met$neurotransmitter_color[match_idx] else NA
   })


   # Create the ggplot plot with the ordered 'cluster' variable -- ALL

   plot_df = subset(master_df, cluster != "REF")
   plot_df=na.omit(plot_df)
   if("REF"%in%levels(plot_df$cluster)){
      plot_df <- plot_df %>%
         mutate(cluster = fct_drop(cluster, "REF"))
   }else{
      plot_df<-plot_df
   }
   if("0"%in%levels(plot_df$cluster)){
      plot_df$cluster=factor(plot_df$cluster,levels=(rev(c(0,1:(length(unique(plot_df$cluster))-1)))))
   }else{
      plot_df <- plot_df %>%
         mutate(cluster = fct_relevel(cluster, sort))
   }
   # Remove rows with NA values in the 'index' column
   plot_df <- plot_df %>%
      filter(!is.na(cluster))

   #Adding in ability to report out classes as colored facets
   unique_classes <- unique(plot_df$class) # Get unique class values
   unique_colors <- plot_df$class_color[match(unique_classes, plot_df$class)] # Match the colors to unique classes
   facet_colors <- setNames(unique_colors, unique_classes)
   unique_spaces <- sapply(seq_along(unique_classes), function(x) paste(rep(" ", x), collapse = ""))
   plot_df$class_space <- factor(plot_df$class, levels = unique_classes, labels = unique_spaces)

   ###########
   #PLOTS
   if (!is.null(meta)) {
      p1 <- ggplot(data = plot_df,
                   aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue, 10)))) +
         geom_point(aes(size = as.numeric(prop), color = sig, fill = class)) + # Use 'fill' for class in legend
         scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
         labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
         scale_color_manual(values = c("black", "red")) + # Keep original color scale for 'sig'
         scale_fill_manual(values = facet_colors) + # Map class to fill color for the legend
         theme_bw() +
         theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text.x = element_text(size = 8) # Keep space-based facet labels
         ) +
         labs(color = "Enrichment \n(cutoff p < 0.05)",
              size = "Proportion of Markers",
              alpha = "Enrichment \n(-log10(p))",
              fill = "Class") + # Add 'Class' to the fill legend
         facet_grid2(~class_space, scales = "free_x", space = "free",
                     strip = strip_themed(
                        background_x = elem_list_rect(fill = facet_colors) # Keep facet colors
                     )) +
         guides(fill = guide_legend(override.aes = list(color = facet_colors))) # Override fill legend colors

   }else{
      p1 <- ggplot(data = plot_df,
                   aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
         geom_point(aes(size = as.numeric(prop), color = sig)) +
         scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
         labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
         scale_color_manual(values = c("black", "red")) +
         theme_bw() +
         scale_y_discrete(limits = rev(levels(plot_df$cluster)))+
         theme(axis.text.x = element_text(angle = 90, hjust = 1),
               axis.text.y = element_text(size = 8),
               plot.title = element_text(size = 10),
               strip.text = element_text(size = 8))+
         labs(color = "Significant \n(p < 0.05)", size = "Proportion of Markers",alpha = "Enrichment \n(-log10(p))")
   }

   # Barplot showing total available markers and proportion covered
   p2 <- ggplot(data =plot_df, aes(x = celltype, y = total_marker)) +
      geom_bar(stat = "identity", fill = "Black") +
      labs(title = "Total Amount of # Markers in Database", y = "Set Size") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))


###############################################
# NEED TO ADD COLOR CODE TO P3 and P4 but turn around the if loops to match above....








   # Create the ggplot plot with the ordered 'cluster' variable -- ONLY SIGNFICANT without facet grid
   if (!is.null(meta) ) {
      p3 <- ggplot(data = plot_df,
                   aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue, 10)))) +
         geom_point(aes(size = as.numeric(prop), color = sig, fill = class)) + # Use 'fill' for class in legend
         scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
         labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
         scale_color_manual(values = c("white", "red")) + # Keep original color scale for 'sig'
         scale_fill_manual(values = facet_colors) + # Map class to fill color for the legend
         theme_bw() +
         theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text.x = element_text(size = 8) # Keep space-based facet labels
         ) +
         labs(color = "Enrichment \n(cutoff p < 0.05)",
              size = "Proportion of Markers",
              alpha = "Enrichment \n(-log10(p))",
              fill = "Class") + # Add 'Class' to the fill legend
         facet_grid2(~class_space, scales = "free_x", space = "free",
                     strip = strip_themed(
                        background_x = elem_list_rect(fill = facet_colors) # Keep facet colors
                     )) +
         guides(fill = guide_legend(override.aes = list(color = facet_colors))) # Override fill legend colors

   }else{
      p3 <- ggplot(data = plot_df,
                   aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
         geom_point(aes(size = as.numeric(prop), color = sig)) +
         scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
         labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
         scale_color_manual(values = c("white", "red")) +
         scale_y_discrete(limits = rev(levels(plot_df$cluster)))+
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1),
               axis.text.y = element_text(size = 8),
               plot.title = element_text(size = 10),
               strip.text = element_text(size = 8))+
         labs(color = "Significant \n(p < 0.05)", size = "Proportion of Markers",alpha = "Enrichment \n(-log10(p))")

   }



   # Create the ggplot plot with the ordered 'cluster' variable -- ONLY TOP SIGNFICANT without facet grid
   if (!is.null(meta) ) {
      p4 <- ggplot(data = plot_df,
                   aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue, 10)))) +
         geom_point(aes(size = as.numeric(prop), color = top_sig, fill = class)) + # Use 'fill' for class in legend
         scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
         labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
         scale_color_manual(values = c("white", "red")) + # Keep original color scale for 'sig'
         scale_fill_manual(values = facet_colors) + # Map class to fill color for the legend
         theme_bw() +
         theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text.x = element_text(size = 8) # Keep space-based facet labels
         ) +
         labs(color = "Enrichment \n(cutoff p < 0.05)",
              size = "Proportion of Markers",
              alpha = "Enrichment \n(-log10(p))",
              fill = "Class") + # Add 'Class' to the fill legend
         facet_grid2(~class_space, scales = "free_x", space = "free",
                     strip = strip_themed(
                        background_x = elem_list_rect(fill = facet_colors) # Keep facet colors
                     )) +
         guides(fill = guide_legend(override.aes = list(color = facet_colors))) # Override fill legend colors
   }else{
      p4 <- ggplot(data = plot_df,
                   aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
         geom_point(aes(size = as.numeric(prop), color = top_sig)) +
         scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
         labs(title = "The Most Significant Cell Type by Cluster", x = " ", y = "Cluster") +
         scale_color_manual(values = c("white","red")) +
         scale_y_discrete(limits = rev(levels(plot_df$cluster)))+
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1),
               axis.text.y = element_text(size = 8),
               plot.title = element_text(size = 10),
               strip.text = element_text(size = 8))+
         labs(color = "Significant \n(p < 0.05)", size = "Proportion of Markers",alpha = "Enrichment \n(-log10(p))")

   }



   ann@results$marker_based$markers_barplot=p2
   ann@results$marker_based$dotplot_all=p1
   ann@results$marker_based$dotplot_sig=p3
   ann@results$marker_based$dotplot_best=p4
   return(ann)
}
