#' Creates marker-based visualizations for cell type annotation
#'
#' Generates various visualizations based on marker-based cell type annotation results
#' stored in the `ann` object. Requires additional cell type metadata from the `meta`
#' data frame. Creates dot plots, bar plots, and facilitates faceting by cell type class.
#'
#' @param ann An SAHA analysis object containing annotation results.
#' @param meta A data frame or path to a CSV file containing cell type metadata.
#' @param facet Logical indicating whether to facet plots by cell type class.
#'
#' @return The `ann` object with additional elements in the `results$marker_based` slot:
#'   - `markers_barplot`: A bar plot showing the total number of markers per cell type.
#'   - `dotplot_all`: A dot plot showing all markers for each cell type and cluster.
#'   - `dotplot_sig`: A dot plot showing only significant markers for each cell type and cluster.
#'   - `dotplot_best`: A dot plot showing the most significant marker for each cell type and cluster.
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% group_by, filter, slice_min, ungroup, mutate, select, as.data.frame
#' @importFrom stringr str_remove_all, str_replace, str_replace_all, if_else
#' @importFrom data.table full_join
#' @importFrom ggplot2 ggplot, geom_point, scale_size_continuous, labs, scale_color_manual, theme_bw, theme, element_text, facet_grid, geom_bar, geom_text, coord_cartesian
#'
#' @export

Create_MarkerBased_Viz <- function (ann, meta = NULL, facet = F){
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
      warning("meta argument is NULL. Function might not work as expected.")
      meta <- data.frame()  # Create an empty data frame
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


   p1 <- ggplot(data = plot_df,
                aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
      geom_point(aes(size = as.numeric(prop), color = sig)) +
      scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
      scale_color_manual(values = c("black", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))+
      labs(color = "Enrichment \n(cutoff p < 0.05)", size = "Proportion of Markers",alpha = "Enrichment \n(-log10(p))")

   if (facet == TRUE & length(meta) > 0 ) {
      p1 <- p1+facet_grid(~class, scales="free_x", space = "free")
   }else{warning("You are attempting to facet the Marker-based plot, however there is a problem with the metadata file loaded. Please refer to the SAHA manual or open an issue on GitHub.")}

   # Barplot showing total available markers and proportion covered
   p2 <- ggplot(data =plot_df, aes(x = celltype, y = total_marker)) +
      geom_bar(stat = "identity", fill = "Black") +
      geom_text(aes(label = total_marker), vjust = -0.5, size = 3) +
      labs(title = "Total Amount of %s Markers in %s Database", y = "Proportion") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8)) +
      coord_cartesian(ylim = c(0, max(master_df$total_marker) * 1.05))



   # Create the ggplot plot with the ordered 'cluster' variable -- ONLY SIGNFICANT without facet grid
   p3 <- ggplot(data = plot_df,
                aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
      geom_point(aes(size = as.numeric(prop), color = sig)) +
      scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
      scale_color_manual(values = c("white", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))+
      labs(color = "Enrichment \n(cutoff p < 0.05)", size = "Proportion of Markers",alpha = "Enrichment \n(-log10(p))")

   if (facet == TRUE & length(meta) > 0 ) {
      p3 <- p3+facet_grid(~class, scales="free_x", space = "free")
   }



   # Create the ggplot plot with the ordered 'cluster' variable -- ONLY TOP SIGNFICANT without facet grid
   p4 <- ggplot(data = plot_df,
                aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
      geom_point(aes(size = as.numeric(prop), color = top_sig)) +
      scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(title = "The Most Significant Cell Type by Cluster", x = " ", y = "Cluster") +
      scale_color_manual(values = c("white","red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))+
      labs(color = "Enrichment \n(cutoff p < 0.05)", size = "Proportion of Markers",alpha = "Enrichment \n(-log10(p))")
   if (facet == TRUE & length(meta) > 0 ) {
      p4 <- p4+facet_grid(~class, scales="free_x", space = "free")
   }



   ann@results$marker_based$markers_barplot=p2
   ann@results$marker_based$dotplot_all=p1
   ann@results$marker_based$dotplot_sig=p3
   ann@results$marker_based$dotplot_best=p4
   return(ann)
}
