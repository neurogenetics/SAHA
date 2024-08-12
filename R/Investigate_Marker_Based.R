#' Investigates marker overlap between query cluster and database cell types
#'
#' Generates visualizations or summaries to explore marker overlap between a specified
#' cluster in the query data and cell types within the database.
#'
#' @param ann An object containing `ann1$query` and `ann1$db` data frames.
#' @param query_cluster The cluster of interest in the query data.
#' @param db_cell_type A specific cell type within the database (optional).
#'   Set to "all" for a summary of cell types within the query cluster.
#' @param plot The type of visualization to generate (optional):
#'   - "venn": Creates a Venn diagram showing marker overlap.
#'   - "stacked": Creates a stacked bar chart of marker distribution.
#'
#' @return No return value (prints informative messages and plots).
#'
#' @export

Investigate_Marker_Based <- function(ann, query_cluster, db_cell_type = NULL, plot = NULL) {

   #1. Read in your queried data as a df or csv
   query_data <- ann@ann1$query
   marker_data <- ann@ann1$db

   #3. Formatting data
   # Calculate the number of all unique genes
   i = query_data[query_data$cluster == query_cluster, "gene"] # Genes for the specified cluster

   #4. Set up Venn
   if (plot == "venn") {

      # Filter marker_data for cluster_markers and db_cell_markers
      cluster_markers <- unique(query_data[query_data$gene %in% i, "gene"])
      db_cell_markers <- unique(marker_data[marker_data$cluster == db_cell_type, "SYMBOL"])

      # Create a list of sets
      sets <- list(`Query Cluster` = cluster_markers, `DB celltype` = db_cell_markers)

      # Generate the Euler diagram
      euler_venn <- euler(sets)

      # Open a new plotting window
      plot.new()

      # Adjust plot aesthetics
      print(plot(euler_venn,
                 labels = c(paste("Cluster", query_cluster), paste("DB:", db_cell_type)),
                 quantities = TRUE,
                 font.main = 14,
                 font.sub = 12,
                 shape = "circle",
                 fill = list(fill = c("pink", "gray"))))
   }

   #5. Set up Stacked
   if (plot == "stacked") {
      # Combine the two sets of data
      combined_data <- t(rbind(
         summary(factor(marker_data$cluster)),
         summary(factor(marker_data[marker_data$SYMBOL %in% i, "cluster"], levels = levels(factor(marker_data$cluster))))
      ))
      combined_data[,1]=combined_data[,1]-combined_data[,2]
      # Convert to data frame
      plot_data <- data.frame(
         celltype = factor(rownames(combined_data)),
         Database = as.numeric(combined_data[, 1]),
         Cluster = as.numeric(combined_data[, 2])
      )
      # Reshape the data into long format
      plot_data_long <- plot_data %>%
         tidyr::pivot_longer(cols = c(Database, Cluster), names_to = "variable", values_to = "value")
      # Plot stacked barplot
      print(ggplot(plot_data_long, aes(x = celltype, y = value, fill = fct_rev(variable))) +
               geom_bar(stat = "identity") +
               labs(title = "Stacked Barplot", y = "Count", x = "Cell Type", fill = "") +
               theme_bw() +
               theme(axis.text.x = element_text(angle = 90,hjust=1),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10),
                     strip.text = element_text(size = 8)) +
               scale_fill_manual(values = c("Cluster" = "red", "Database" = "black"), labels = c("Database" = "Database", "Cluster" = paste0("Cluster ", query_cluster))))
   }

   #6. Create if statements to help guide the user
   if (is.null(plot)){
      print("You can select between two different visualization methods, 'venn' and 'stacked'")
   }

   if (is.null(query_cluster)) {
      print("To use this function, you need to at least select one cluster for visualization.")
   }

   if (is.null(db_cell_type)) {
      print(paste("You can also select a cell type found within the cluster for further analysis. Try db_cell_type = '",marker_data$cluster[1],"'."))
   }

   if (!is.null(db_cell_type) && db_cell_type == "all") {
      cat("Summary of cell types found in cluster ", query_cluster, ".\n", sep = "")
      cat("Pick one of the cell types to continue. \n")
      print(summary(factor(marker_data[marker_data$official.gene.symbol %in% j, "cell.type"], levels(factor(marker_data$cell.type)))))
   }

   if (!is.null(db_cell_type) && db_cell_type != "all") {
      query_data_names_type <- marker_data[(marker_data$cluster == db_cell_type),"SYMBOL"][marker_data[(marker_data$cluster == db_cell_type),"SYMBOL"]%in% i]
      cat("Markers for", db_cell_type, "within Cluster", query_cluster, ":", paste(sort(query_data_names_type), collapse = ", "), "\n")
   }
}
