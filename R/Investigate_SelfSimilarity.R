#' Analyzes self-similarity between clusters in marker space
#'
#' Compares marker overlap between two specified clusters in the marker similarity matrix. Offers the option to return a data frame of shared features (optional).
#'
#' @param ann An object containing the `ann1$Marker_selfsim_matrix` data frame.
#' 
#' @param cluster1 The name of the first cluster to compare.
#' 
#' @param cluster2 The name of the second cluster to compare.
#' 
#' @param shared_df Logical indicating whether to return a data frame of shared features (TRUE). If FALSE (default), the function generates a Venn diagram.
#'
#' @return A Venn diagram or a data frame of shared features between cluste (depending on the `shared_df` argument).
#'
#' @export

Investigate_Self_Similarity <- function(ann, cluster1, cluster2,shared_df=NULL) {
   # Ensure data is numeric
   df <- as.matrix(ann@ann1$Marker_selfsim_matrix)
   df[is.na(df)] <- 0

   # Extract rows for the specified clusters
   cluster1_data <- df[paste(cluster1), ]
   cluster2_data <- df[paste(cluster2), ]

   # Find shared features (where both rows have 1)
   bound=cbind(cluster1_data,cluster2_data)
   bound=data.frame(bound)
   venn_input_1=rownames(bound[bound$cluster1_data==1,])
   venn_input_2=rownames(bound[bound$cluster2_data==1,])

   # Create a list of sets
   sets <- list(`Cluster_1` = venn_input_1, `Cluster_2` = venn_input_2)

   # Generate the Euler diagram
   euler_venn <- euler(sets)

   # Open a new plotting window
   plot.new()

   # Adjust plot aesthetics
   print(plot(euler_venn,
              labels = c(paste(cluster1), cluster2),
              quantities = TRUE,
              font.main = 14,
              font.sub = 12,
              shape = "circle",
              fill = list(fill = c("pink", "gray"))))
   if (!is.null(shared_df)) {
      bound_df=bound[rowSums(bound)>=1,]
      return(bound_df)
   }

}
