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
         group_by(cluster)   # Filter for rows where gene and cluster values match

      p1 = aggregate(gene~cluster, gene_count_by_cluster, length)
      p1 <- p1 %>%
         mutate(cluster = factor(cluster, levels = sort(levels(cluster)))) %>%
         arrange(cluster)

      print(p1)
      ggplot(p1,aes(x=cluster,y=gene))+
         geom_bar(stat="identity")+
         theme_bw() + theme(axis.text = element_text(hjust = 0.9),
    axis.text.x = element_text(angle = 90))
   }else{
      ann1=ann@ann1
      df = data.frame(ann1$query[,c("gene","cluster")])

      # Count genes present in varfeat
      gene_counts <- df %>%
         mutate(in_varfeat = gene %in% varfeat)

      p2 = aggregate(in_varfeat~cluster, gene_counts, sum)
      p2 <- p2 %>%
         mutate(cluster = factor(cluster, levels = sort(levels(cluster)))) %>%
         arrange(cluster)
      print(p2)
      # Print the results
      ggplot(p2,aes(x=cluster,y=in_varfeat))+
         geom_bar(stat="identity")+
         geom_hline(yintercept = length(varfeat),linetype = 3)+
         theme_bw() + theme(axis.text = element_text(hjust = 0.9),
                            axis.text.x = element_text(angle = 90))
   }
}
