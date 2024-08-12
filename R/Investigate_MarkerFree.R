#' Compares gene expression between query cluster and database cell type
#'
#' Generates a scatter plot comparing gene expression levels between a specified query cluster
#' and a database cell type. Highlights top expressed genes and outliers.
#'
#' @param ann An object containing normalized gene expression data in `ann@results$marker_free$norm_merge`.
#' @param query_cluster The cluster of interest in the query data.
#' @param db_cell_type The cell type of interest in the database.
#'
#' @return A scatter plot visualizing the gene expression comparison.
#'
#' @export

Investigate_MarkerFree = function(ann,query_cluster, db_cell_type){
   # Extract gene expression data for the cluster of interest
   x_data <- ann@results$marker_free$norm_merge[,paste0("query.",query_cluster)]

   # Extract gene expression data for the best match cluster
   y_data <- ann@results$marker_free$norm_merge[,paste0("db.",db_cell_type)]


   # Combine the data into a data frame for plotting
   plot_data <- data.frame(
      Gene = rownames(ann@results$marker_free$norm_merge),
      Expression_in_Interest_Cluster = x_data,
      Expression_in_Best_Match_Cluster = y_data,
      Difference = abs(x_data - y_data)
   )

   top_10_genes <- plot_data[order(-plot_data$Expression_in_Best_Match_Cluster), ][1:10, ]

   top_5_outliers <- plot_data[order(-plot_data$Difference), ][1:5, ]

   # Create the scatter plot
   p1 = ggplot(plot_data, aes(x = Expression_in_Interest_Cluster, y = Expression_in_Best_Match_Cluster)) +
      geom_point() +
      geom_label_repel(data = top_10_genes, aes(label = Gene), max.overlaps = 20, fill = "lightblue") +
      geom_label_repel(data = top_5_outliers, aes(label = Gene), max.overlaps = 20, fill = "#f78a5c") +
      labs(
         title = paste("Gene Expression in", query_cluster, "vs", db_cell_type),
         x = paste("Expression in", query_cluster),
         y = paste("Expression in", db_cell_type)
      ) +
      theme_bw()
   print(p1)
}
