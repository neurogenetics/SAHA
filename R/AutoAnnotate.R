#' Automatically annotate query clusters with best matching cell types
#'
#' This function assigns cell type annotations to clusters in a query dataset using one of three methods:
#' 
#' - `"Markers"`: Uses marker-based annotation from the `ann@ann2` data frame. It selects the cell type with the minimum p-value for each query cluster, filtering for p < 0.05. Clusters without significant matches are labeled as `"INCONCLUSIVE"`.
#'
#' - `"AvgExp"`: Uses a marker-free method by computing the highest Pearson correlation between average expression profiles (`ann@results$marker_free$corr`). This method assigns the reference cell type with the highest correlation for each cluster.
#'
#' - `"Both"`: Combines both methods above. It generates a consensus annotation by comparing the marker-based and correlation-based matches. If both methods agree (after whitespace removal), the consensus is labeled as `"MATCH"`; otherwise, it is labeled `"DISAGREEMENT"`. Final assignments are set to the matched cell type if consensus is reached, and `"INCONCLUSIVE"` otherwise.
#'
#' @param ann An annotation object containing both:
#'   - `ann@ann2`: A data frame with cluster-level marker enrichment results including columns `cluster`, `celltype`, `pvalue`, and `prop`.
#'   - `ann@results$marker_free$corr`: A matrix of average expression correlations between query and reference cell types.
#'
#' @param data_type Character string specifying the annotation mode. Options:
#'   - `"Markers"` (default): Perform annotation using marker enrichment data.
#'   - `"AvgExp"`: Use correlation of average expression.
#'   - `"Both"`: Use both methods and generate a consensus.
#'   If `NULL`, the function will print a warning and return nothing.
#'
#' @details Clusters labeled `"REF"` are excluded from marker-based analysis. If no significant p-values are found for a cluster,it is labeled `"INCONCLUSIVE"`. In consensus mode, whitespace is stripped for string comparison.
#'
#' @importFrom dplyr %>% group_by filter slice_min ungroup mutate select arrange as.data.frame
#' @importFrom stringr str_remove_all str_replace str_replace_all if_else
#' @importFrom data.table full_join
#'
#' @export

AutoAnnotate = function(ann, data_type=NULL){
  if(is.null(data_type)){
    print("We gently discourage auto-annotation. If you wish to proceed please select which data_type to auto-annotate using: Markers, AvgExp, or Both")
  }else if(data_type=="Markers"){
    best_match.df <-  ann@ann2 %>%
      group_by(cluster) %>%
      filter(any(pvalue < 0.05) & cluster != "REF") %>%
      slice_min(pvalue) %>%
      ungroup() %>%
      mutate(celltype = ifelse(is.na(celltype), "UNKNOWN", celltype))%>%
      select(cluster,celltype, prop, pvalue)%>%
      as.data.frame()
    unique_clust=unique(ann@ann2$cluster)
    filtered_clusters <- unique_clust[unique_clust != "REF"]
    if (!all(filtered_clusters%in%best_match.df$cluster)) {
      missing_clusters <- setdiff(filtered_clusters, best_match.df$cluster)
      missing_data <- data.frame(
        cluster = missing_clusters,
        celltype = "INCONCLUSIVE",
        prop = NA,
        pvalue = NA)
      best_match.df <- rbind(best_match.df, missing_data)
    }
    colnames(best_match.df)[2]="best_match"
    best_match.df <- best_match.df %>%
      arrange(cluster)
    return(best_match.df)
    
  }else if(data_type=="AvgExp"){
    best_match.df <- data.frame(row.names = rownames(ann@results$marker_free$corr),
                                best_match = apply(ann@results$marker_free$corr, 1, function(row) names(ann@results$marker_free$corr)[which.max(row)]),
                                correlation = apply(ann@results$marker_free$corr, 1, max))
    colnames(best_match.df)[2]="best_match"
    best_match.df$cluster=gsub("^query\\.", "", rownames(best_match.df))
    best_match.df=best_match.df[,c(3,1,2)]
    if (typeof(best_match.df$cluster)!="character") {
      best_match.df$cluster=sort(as.numeric(unique(best_match.df$cluster)),decreasing = F)
    }
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
    
    best_match.Markers$cluster <- as.character(best_match.Markers$cluster)
    best_match.AvgExp$cluster <- as.character(best_match.AvgExp$cluster)
    
    if (length(intersect(best_match.AvgExp$cluster, best_match.Markers$cluster)) == 0) {
      if (any(grepl("RNA\\.", best_match.AvgExp$cluster))) {
        best_match.AvgExp$cluster <- gsub("RNA\\.", "", best_match.AvgExp$cluster)
      }
      if (length(intersect(best_match.AvgExp$cluster, best_match.Markers$cluster)) == 0) {
        warning("No matching cluster names found in marker-based and marker-free analysis. Consider running separately or renaming clusters.") 
      }
      best_matches=suppressMessages(full_join(best_match.AvgExp,y = best_match.Markers))
      best_matches <- best_matches %>%
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
      return(best_matches)
    }else{
      best_matches=suppressMessages(full_join(best_match.AvgExp,y = best_match.Markers))
      best_matches <- best_matches %>%
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
      return(best_matches)
    }
  }
}
