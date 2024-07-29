# SAHA
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#################################################
#Rework so this is the quickstart option########
#################################################


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
      return(best_match.df)

   }else if(data_type=="AvgExp"){
      best_match.df <- data.frame(row.names = rownames(ann@results$marker_free$corr),
                                  best_match = apply(ann@results$marker_free$corr, 1, function(row) names(ann@results$marker_free$corr)[which.max(row)]),
                                  correlation = apply(ann@results$marker_free$corr, 1, max))
      colnames(best_match.df)[2]="best_match"
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
      best_matches=full_join(best_match.AvgExp,y = best_match.Markers)
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
      return(best_matches)
   }
}
