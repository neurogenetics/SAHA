#' Marker Richness Function
#'
#' This function calculates the richness of markers within clusters, either using all genes or a specified set of variable features.
#'
#' @param ann An object containing annotation data.
#' @param varfeat A vector of variable features (genes) to consider. If NULL, all genes are used.
#' @return A data frame with gene counts by cluster.
#' @importFrom dplyr group_by count mutate summarize
#' @export

Marker_Richness = function(ann, varfeat = NULL){
   if (is.null(varfeat)) {
      ann1=ann@ann1
      df = data.frame(ann1$query[,c("gene","cluster")])

      gene_count_by_cluster <- df %>%
         group_by(cluster) %>%  # Filter for rows where gene and cluster values match
         count()  # Count the number of rows after filtering

      print(gene_count_by_cluster)
   }else{
      ann1=ann@ann1
      df = data.frame(ann1$query[,c("gene","cluster")])

      # Count genes present in varfeat
      gene_counts <- df %>%
         mutate(in_varfeat = gene %in% varfeat) %>%  # Create a new column indicating presence
         group_by(cluster)%>%
         summarize(
            all_genes = n(),
            genes_in_varfeat=sum(in_varfeat)  # Count occurrences of each gene and presence flag
         )
      # Print the results
      print(gene_counts)
   }
}
