#' Initializes marker sets for SAHA analysis
#'
#' This function filters marker genes from both the query and database datasets
#' in a SAHA object based on user-defined thresholds. It prepares the data for
#' downstream marker-based annotation by applying cutoffs for adjusted p-value,
#' fold change, sensitivity, and specificity.
#'
#' @param ann An object containing `query$Markers` and `db$Markers` data frames.
#' 
#' @param p_thresh Adjusted p-value threshold for marker selection.
#' 
#' @param FC_thresh Fold change threshold for marker selection.
#' 
#' @param sens_thresh Sensitivity threshold for marker selection.
#' 
#' @param spec_thresh Specificity threshold for marker selection.
#'
#' @return The input `ann` object with updated `ann1$query` and `ann1$db` components
#'   containing filtered marker sets.
#'
#' @importFrom dplyr subset
#'
#' @export
#'
Initialize_Markers <- function(ann, p_thresh = 0.05, FC_thresh = 1.5, sens_thresh = 0.25, spec_thresh = 0.75) {
   # Assuming your dataframe is named 'df' and has columns 'cluster' and 'log2FC'
   # Query dataset processing
   if (all(ann@query$Markers$avg_log2FC >= 0)) {
      ann_subset = subset(ann@query$Markers, p_val_adj < p_thresh & avg_log2FC > log(FC_thresh, base = 2) & pct.1 > sens_thresh & pct.2 < spec_thresh)
   } else {
      warning("Query dataset contains markers with negative or missing log2FC values. Proceeding with subset.\n")
      ann_subset = subset(ann@query$Markers, p_val_adj < p_thresh & avg_log2FC > log(FC_thresh, base = 2) & pct.1 > sens_thresh & pct.2 < spec_thresh)
   }

   # Summary of filtered query dataset
   query_summary <- aggregate(gene ~ cluster, ann_subset, length)
   cat(sprintf("Query dataset: %g clusters, median %g markers (range: %g to %g per cluster).\n",
               length(unique(ann_subset$cluster)), median(query_summary$gene), min(query_summary$gene), max(query_summary$gene)))

   # DB dataset processing
   if (is.null(ann@db$Markers$p_val_adj)) {
      db_subset <- ann@db$Markers
   } else {
      db_subset <- subset(ann@db$Markers, p_val_adj < p_thresh & avg_log2FC > log(FC_thresh, base = 2) & pct.1 > sens_thresh & pct.2 < spec_thresh)
   }

   # Summary of filtered DB dataset
   db_summary <- aggregate(gene ~ cluster, db_subset, length)
   cat(sprintf("DB dataset: %g clusters, median %g markers (range: %g to %g per cluster).\n",
               length(unique(db_subset$cluster)), median(db_summary$gene), min(db_summary$gene), max(db_summary$gene)))

   # Drop warning for clusters
   if (length(unique(ann@query$Markers$cluster)) != length(unique(ann_subset$cluster))) {
      warning("Some clusters were dropped from the query dataset based on specified thresholds.\n")
   }

   # Save filter parameters in ann
   ann@params$markers <- list(p_thresh = p_thresh, FC_thresh = FC_thresh, sens_thresh = sens_thresh, spec_thresh = spec_thresh)

   # Update ann with subsets
   ann@ann1$query <- ann_subset
   db_subset$cluster <- gsub("-", " ", db_subset$cluster)
   ann@ann1$db <- db_subset

   return(ann)
}
