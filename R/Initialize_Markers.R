#' Initializes marker sets for SAHA analysis
#'
#' Filters markers from query and database datasets based on specified thresholds.
#'
#' @param ann An object containing `query$Markers` and `db$Markers` data frames.
#' @param p_thresh Adjusted p-value threshold for marker selection.
#' @param FC_thresh Fold change threshold for marker selection.
#' @param sens_thresh Sensitivity threshold for marker selection.
#' @param spec_thresh Specificity threshold for marker selection.
#'
#' @return The input `ann` object with updated `ann1$query` and `ann1$db` components
#'   containing filtered marker sets.
#'
#' @importFrom dplyr subset
#'
#' @export

Initialize_Markers <- function(ann, p_thresh = 0.05, FC_thresh = 1.5, sens_thresh = 0.25, spec_thresh = 0.75){
   # Assuming your dataframe is named 'df' and has columns 'cluster' and 'log2FC'
   # QC reports
   if (all(ann@query$Markers$avg_log2FC >=0)) {
      ann_subset=subset(ann@query$Markers,(p_val_adj < p_thresh &
                                              avg_log2FC > log(FC_thresh,base = 2) &
                                              pct.1 > sens_thresh &
                                              pct.2 < spec_thresh))
   }else{

      cat("WARNING: Loaded query dataset contains markers with either negative or NA log2FC values. We gently discourage the use of negative or absent markers in cluster identification, however the pipeline will proceed regardless.\n")
      ann_subset=subset(ann@query$Markers,(p_val_adj < p_thresh &
                                              avg_log2FC > log(FC_thresh,base = 2) &
                                              pct.1 > sens_thresh &
                                              pct.2 < spec_thresh))
   }

   marker_summary <- aggregate(gene~cluster, ann_subset, length)

   cat(paste("Loaded query dataset contains", length(unique(ann_subset$cluster)),"unique clusters with a median of", summary(marker_summary$gene)[3], "markers per cluster and a range from",summary(marker_summary$gene)[1],"to",summary(marker_summary$gene)[6],"markers per cluster.\n"))

   db_subset = subset(ann@db$Markers,(p_val_adj < p_thresh &
                                         avg_log2FC > log(FC_thresh,base = 2) &
                                         pct.1 > sens_thresh &
                                         pct.2 < spec_thresh))


   db_marker_summary <- aggregate(gene~cluster, db_subset, length)

   cat(paste(
      "Loaded db dataset contains", length(unique(db_subset$cluster)),"unique clusters with a median of", summary(db_marker_summary$gene)[3], "markers per cluster and a range from",summary(db_marker_summary$gene)[1],"to",summary(db_marker_summary$gene)[6],"markers per cluster.\n"))
   cat("If these summaries match your expetaction for running the SAHA pipeline, no further action is needed. If you would like to re-define markers based on unique p_value or log2FC cutoffs, please re-run Initialize_SAHA() with arguments for thresholds. Otherwise, db markers can be specified using Tune_Markers(). \n")

   ##check how
   if (length(unique(ann@query$Markers$cluster))!=length(unique(ann_subset$cluster))) {
      cat("WARNING: by subsetting your original Markers dataframe with this function, some clusters have dropped out from analysis.")
   }
   #overwrite ann
   ann@ann1$query=ann_subset
   db_subset$cluster <- gsub("-", " ", db_subset$cluster)
   ann@ann1$db=db_subset
   return(ann)
}
