#' Automatically annotate query clusters with best matching cell types 
#'
#' @param ann An annotation object.
#' @param data_type Character string specifying the mode ("Markers", "AvgExp", or "Both").
#' 
#' @importFrom dplyr %>% group_by filter slice_min ungroup mutate select arrange as.data.frame full_join coalesce
#' @importFrom stringr str_remove_all str_replace str_replace_all if_else 
#' @export
AutoAnnotate = function(ann, data_type=NULL){
   
   # ROBUST CLEANING: Removes noise words and punctuation regardless of position
   clean_string <- function(x) {
      x <- as.character(x)
      # 1. Replace all punctuation (dots, slashes, underscores) with spaces
      x <- gsub("[._/]", " ", x)
      # 2. Remove specific noise words as standalone words (using word boundaries \\b)
      # This catches "RNA", "query", "db", etc., even if nested
      x <- gsub("\\b(query|RNA|SCT|integrated|db)\\b", "", x, ignore.case = TRUE)
      # 3. Collapse multiple spaces and trim
      x <- trimws(gsub("\\s+", " ", x))
      return(x)
   }
   
   if(is.null(data_type)){
      message("Please select a data_type: Markers, AvgExp, or Both")
      return(NULL)
      
   } else if(data_type=="Markers"){
      best_match.df <-  ann@ann2 %>%
         group_by(cluster) %>%
         filter(any(pvalue < 0.05) & cluster != "REF") %>%
         slice_min(pvalue, n = 1, with_ties = FALSE) %>%
         ungroup() %>%
         mutate(celltype = ifelse(is.na(celltype), "UNKNOWN", celltype)) %>%
         select(cluster, celltype, prop, pvalue) %>%
         as.data.frame()
      
      # Standardize the cluster names for the final table
      best_match.df$cluster <- clean_string(best_match.df$cluster)
      colnames(best_match.df)[2] = "best_match"
      return(best_match.df %>% arrange(cluster))
      
   } else if(data_type=="AvgExp"){
      corr_mtx <- ann@results$marker_free$corr
      best_match.df <- data.frame(
         row.names = rownames(corr_mtx),
         best_match = apply(corr_mtx, 1, function(row) colnames(corr_mtx)[which.max(row)]),
         correlation = apply(corr_mtx, 1, max)
      )
      # Standardize the cluster names for the final table
      best_match.df$cluster <- clean_string(rownames(best_match.df))
      # Standardize the result column as well
      best_match.df$best_match <- clean_string(best_match.df$best_match)
      
      best_match.df <- best_match.df[,c("cluster", "best_match", "correlation")]
      return(best_match.df)
      
   } else if(data_type=="Both"){
      # 1. Prepare Marker-based results
      best_match.Markers <- ann@ann2 %>%
         group_by(cluster) %>%
         filter(pvalue < 0.05) %>%
         slice_min(pvalue, n = 1, with_ties = FALSE) %>%
         ungroup() %>%
         select(cluster, celltype) %>%
         as.data.frame()
      
      # 2. Prepare Expression-based results
      corr_mtx <- ann@results$marker_free$corr
      best_match.AvgExp <- data.frame(
         cluster = rownames(corr_mtx),
         best_match_raw = apply(corr_mtx, 1, function(row) colnames(corr_mtx)[which.max(row)]),
         stringsAsFactors = FALSE
      )
      
      # 3. Generate standardized Join IDs (This is what merges the 25 rows into 13)
      best_match.Markers$join_id <- clean_string(best_match.Markers$cluster)
      best_match.AvgExp$join_id <- clean_string(best_match.AvgExp$cluster)
      
      # 4. Perform the join
      best_matches <- suppressMessages(dplyr::full_join(best_match.AvgExp, best_match.Markers, by = "join_id"))
      
      # 5. Clean up the actual cell type assignments (removing "RNA", etc.)
      best_matches <- best_matches %>%
         mutate(
            marker_free_clean = clean_string(best_match_raw),
            marker_based_clean = clean_string(celltype),
            final_cluster = join_id
         )
      
      # 6. Consensus Logic
      best_matches <- best_matches %>%
         mutate(
            # Remove all spaces just for the final logical comparison
            norm_marker = str_remove_all(marker_based_clean, "\\s+"),
            norm_avg = str_remove_all(marker_free_clean, "\\s+"),
            # MATCH if they are the same after cleaning
            consensus = if_else(!is.na(norm_marker) & !is.na(norm_avg) & norm_marker == norm_avg, "MATCH", "DISAGREEMENT"),
            final_output = if_else(consensus == "MATCH", marker_based_clean, "INCONCLUSIVE")
         )
      
      # 7. Final Clean and Format
      best_matches[] <- lapply(best_matches, as.character)
      best_matches[is.na(best_matches)] <- "INCONCLUSIVE"
      
      res <- best_matches %>%
         select(final_cluster, marker_based_clean, marker_free_clean, consensus, final_output) %>%
         arrange(final_cluster)
      
      colnames(res) <- c("cluster", "marker_based", "marker_free", "consensus", "best_match")
      
      return(res)
   }
}
