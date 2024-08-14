#' Auto-annotates query clusters with best matching database cell types
#'
#' This function offers three auto-annotation options based on the `data_type` argument:
#'   - "Markers": Uses the minimum p-value from the annotation data (`ann@ann2`)
#'     to identify the best matching cell type for each query cluster.
#'   - "AvgExp": Uses the maximum correlation coefficient from the average expression
#'     correlation matrix (`ann@results$marker_free$corr`) for auto-annotation.
#'   - "Both": Combines both marker-based and average expression approaches,
#'     providing marker-based and marker-free best matches and a consensus annotation.
#'
#' @param ann An object containing marker data (`ann@ann2`) and expression data (`ann@results$marker_free$corr`).
#' @param data_type The type of data to use for auto-annotation (optional):
#'   - "Markers" (default): Uses marker-based p-values.
#'   - "AvgExp": Uses average expression correlation.
#'   - "Both": Combines both marker and expression data.
#'
#' @return A data frame containing information about the auto-annotated cell types.
#'   The specific columns depend on the `data_type` argument.
#'
#' @importFrom dplyr %>% group_by, filter, slice_min, ungroup, mutate, select, as.data.frame
#' @importFrom stringr str_remove_all, str_replace, str_replace_all, if_else
#' @importFrom data.table full_join
#'
#' @export

AutoAnnotate = function(ann, data_type=NULL){
   if(is.null(data_type)){
      print("We gently discourage auto-annotation. If you wish to proceed please select which data_type to auto-annotate using: Markers, AvgExp, or Both")
   }else if(data_type=="Markers"){
      best_match.df <-  ann@ann2 %>%
         group_by(cluster) %>%
         filter(any(pvalue < 0.05)) %>%
         slice_min(pvalue) %>%
         ungroup() %>%
         mutate(celltype = ifelse(is.na(celltype), "UNKNOWN", celltype))%>%
         select(cluster,celltype, prop, pvalue)%>%
         as.data.frame()
      colnames(best_match.df)[2]="best_match"
      best_match.df$cluster=sort(as.numeric(unique(best_match.df$cluster)),decreasing = F)
      return(best_match.df)

   }else if(data_type=="AvgExp"){
      best_match.df <- data.frame(row.names = rownames(ann@results$marker_free$corr),
                                  best_match = apply(ann@results$marker_free$corr, 1, function(row) names(ann@results$marker_free$corr)[which.max(row)]),
                                  correlation = apply(ann@results$marker_free$corr, 1, max))
      colnames(best_match.df)[2]="best_match"
      best_match.df$cluster=gsub("^query\\.", "", rownames(best_match.df))
      best_match.df=best_match.df[,c(3,1,2)]
      best_match.df$cluster=sort(as.numeric(unique(best_match.df$cluster)),decreasing = F)
      return(best_match.df)
   }else if(data_type=="Both"){
      best_match.Markers <-  ann@ann2 %>%
         group_by(cluster) %>%
         filter(any(pvalue < 0.05)) %>%
         slice_min(pvalue) %>%
         ungroup() %>%
         mutate(celltype = ifelse(is.na(celltype), "UNKNOWN", celltype))%>%
         select(cluster,celltype, prop, pvalue)%>%
         as.data.frame()

      best_match.AvgExp <- data.frame(row.names = rownames(ann@results$marker_free$corr),
                                      best_match = apply(ann@results$marker_free$corr, 1, function(row) names(ann@results$marker_free$corr)[which.max(row)]),
                                      correlation = apply(ann@results$marker_free$corr, 1, max))
      best_match.AvgExp$cluster<- gsub(paste0("^query\\."), "", rownames(best_match.AvgExp))
      best_matches=invisible(full_join(best_match.AvgExp,y = best_match.Markers))
      best_matches <- best_matches %>%
         #select(best_match) %>%
         mutate(best_match_avg = str_replace(best_match, "^db\\.", "")) %>%
         mutate(best_match_avg = str_replace_all(best_match_avg, "\\.", " "))

      best_matches <- best_matches %>%
         mutate(consensus = if_else(str_remove_all(celltype, "\\s+") == str_remove_all(best_match_avg, "\\s+"), "MATCH", "DISAGREEMENT"))
      best_matches <- best_matches %>%
         mutate(final_output = if_else(consensus == "MATCH", celltype, "INCONCLUSIVE"))
      best_matches[is.na(best_matches)]<-"INCONCLUSIVE"
      best_matches <- best_matches %>%
         select(cluster, celltype,best_match_avg,consensus,final_output)
      colnames(best_matches)=c("cluster","marker_based","marker_free","consensus","best_match")
      best_matches$cluster=sort(as.numeric(unique(best_matches$cluster)),decreasing = F)
      return(best_matches)
   }
}
