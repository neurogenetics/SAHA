#' Analyzes marker-based enrichment for each query cluster
#'
#' This function analyzes marker enrichment for each cluster identified in the query data against the database clusters. It calculates the total number of markers per cluster, proportion of markers shared with each database cell type, and hypergeometric test p-values for enrichment significance.
#'
#' @param ann An object containing marker data in `ann1$query` and `ann1$db`.
#'
#' @return The input `ann` object with a new component `ann2` containing a data frame summarizing marker enrichment analysis results. The function also returns a named list containing the enrichment results (`master_df`), query marker data (`query_data`), and database marker data (`marker_data`).
#'
#' @importFrom dplyr %>% group_by, mutate, ungroup, data.frame
#' @importFrom stats dhyper, factor
#'
#' @export
Run_Marker_Based <- function(ann){
   #5. Making our Master Data Frame
   master_df <- data.frame(summary(factor(ann@ann1$db$cluster)))
   master_df$cluster="REF"
   colnames(master_df)[1]<-"total_marker"
   master_df$celltype <- rownames(master_df)

   #6. Loop through each cluster, find their markers, bind to masterdf
   for (i in c(1:length(unique(ann@ann1$query$cluster)))) {
      x=master_df[master_df$cluster=="REF",]
      x$cluster = unique(ann@ann1$query$cluster)[i]
      x$total_marker=SAHA_lookup_cluster(unique(ann@ann1$query$cluster)[i],ann)
      x$celltype = rownames(x)
      master_df=rbind(master_df,x)
   }

   #7. Find proportion
   master_df$prop="NA"
   for (j in c(1:length(unique(ann@ann1$query$cluster)))) {
      master_df[master_df$cluster==unique(ann@ann1$query$cluster)[j],"prop"]=master_df[master_df$cluster==unique(ann@ann1$query$cluster)[j],"total_marker"]/master_df[master_df$cluster=="REF","total_marker"]
   }

   #8. Find p value
   master_df$pvalue = 1

   j=1
   k=1
   suppressWarnings(
      for (j in c(1:length(unique(ann@ann1$query$cluster)))) {
         # Size of possible markers in a given cluster
         A = length(ann@ann1$query[ann@ann1$query$cluster == unique(ann@ann1$query$cluster)[j], "gene"])
         # Size of signature testing (number of genes in pangloa cell type)
         for (k in c(1:length(unique(master_df$celltype)))) {
            B = master_df[master_df$cluster=="REF",][k,"total_marker"]
            # Overlap
            t = master_df[master_df$cluster==unique(ann@ann1$query$cluster)[j],][k,"total_marker"]
            # Length of all possible cluster markers (not just the one testing)
            n = length(unique(ann@ann1$query$gene))
            master_df[master_df$cluster==unique(ann@ann1$query$cluster)[j],][k,"pvalue"]=sum(stats::dhyper(t:B, A, n - A, B))
         }
      })
   #8b. BH-correct p value
   master_df$padj = p.adjust(master_df$pvalue, method = p.adjust.methods, n = length(master_df$pvalue))

      #9. Conditional Facet
   #10. Default everything to F, not significant
   master_df$sig = "F"
   master_df[is.na(master_df$pvalue),"pvalue"]<-1
   master_df[master_df$padj <= 0.05,"sig"]="T" # Label only significant groups

   # Convert the 'cluster' variable to a factor with custom levels in ascending order
   master_df$cluster <- factor(master_df$cluster, levels = rev(unique(master_df$cluster)))

   # Create dataframe that only finds minimum pvalue for each cluster
   master_df <- master_df %>%
      group_by(cluster) %>%
      mutate(top_sig = ifelse(pvalue == min(pvalue) & sig == "T", "T", "F"),
             top_sig_prop = ifelse(top_sig == "F", 0, as.numeric(prop))) %>%
      ungroup()
   master_df <- data.frame(master_df)

   invisible(list(master_df = master_df, query_data = ann@ann1$query, marker_data = ann@ann1$db ))
   ann@ann2=master_df
   return(ann)
}
