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

Create_MarkerBased_Viz <- function (ann, meta, facet){
   if (inherits(meta, "data.frame")) {
      # If query is already a data frame, use it directly
      meta <- meta
   } else {
      # If query is a path to a CSV file, read it into a data frame
      meta <- read.csv(meta)
   }
   master_df = ann@ann2
   #master_df = ann2 #######need to change to a list where master_df is added!!!!!!!!!!!!
   ##
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



   p1 <- ggplot(data = subset(master_df, cluster != "REF"),
                aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
      geom_point(aes(size = as.numeric(prop), color = sig)) +
      scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
      scale_color_manual(values = c("black", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))

   if (facet == TRUE) {
      p1 <- p1+facet_grid(~class, scales="free_x", space = "free")
   }

   # Barplot showing total available markers and proportion covered
   p2 <- ggplot(data = subset(master_df, cluster == "REF"), aes(x = celltype, y = total_marker)) +
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
   p3 <- ggplot(data = subset(master_df, cluster != "REF"),
                aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
      geom_point(aes(size = as.numeric(prop), color = sig)) +
      scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
      scale_color_manual(values = c("white", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))
   if (facet == TRUE) {
      p3 <- p3+facet_grid(~class, scales="free_x", space = "free")
   }



   # Create the ggplot plot with the ordered 'cluster' variable -- ONLY TOP SIGNFICANT without facet grid
   p4 <- ggplot(data = subset(master_df, cluster != "REF"),
                aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
      geom_point(aes(size = as.numeric(prop), color = top_sig)) +
      scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(title = "The Most Significant Cell Type by Cluster", x = " ", y = "Cluster") +
      scale_color_manual(values = c("white","red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 10),
            strip.text = element_text(size = 8))
   if (facet == TRUE) {
      p4 <- p4+facet_grid(~class, scales="free_x", space = "free")
   }



   ann@results$marker_based$markers_barplot=p2
   ann@results$marker_based$dotplot_all=p1
   ann@results$marker_based$dotplot_sig=p3
   ann@results$marker_based$dotplot_best=p4
   return(ann)
}
