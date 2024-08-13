#' Semi-automated annotation of query clusters
#'
#' This function provides a semi-automated approach for annotating query clusters.
#' Users are prompted to manually assign names to clusters based on visualizations
#' of marker-based or marker-free data.
#'
#' @param ann An object containing SAHA analysis results.
#' @param data_type The type of data to use for annotation ("Markers", "AvgExp", or "Both").
#' @param refine A data frame containing previously refined annotations (optional).
#'
#' @return A data frame containing the original and new cluster names.
#'
#' @importFrom dplyr filter, arrange
#' @importFrom ggpubr ggarrange
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.grabExpr
#'
#' @export

SemiAutoAnnotate = function(ann,data_type=NULL,refine=NULL){
   #prompt through each one
   if(is.null(data_type)){
      print("Thank you for choosing the semi-automated approach. If you wish to proceed please select which data_type to auto-annotate using: Markers, AvgExp, or Both")
   }else if(data_type=="Markers"){

      temp=ann@results$marker_based$dotplot_all
      hand_names=data.frame(unique(temp$data$cluster))
      colnames(hand_names)[1]="old_names"
      hand_names$new_names=""
      if (!is.null(refine)) {
         todo <- hand_names %>%
            filter(old_names %in% unique(refine$cluster[refine$best_match=="INCONCLUSIVE"]))
         tokeep <- hand_names %>%
            filter(old_names %in% unique(refine$cluster[refine$best_match!="INCONCLUSIVE"]))
         tokeep$new_names=refine[refine$best_match!="INCONCLUSIVE","best_match"]
      }else{todo=hand_names}
      for (i in todo$old_names) {
         temp2=temp
         temp2$data=temp$data[temp$data$cluster==i,]
         print(temp2+theme(legend.position = "none"))
         x=readline(paste0("What would you like to name cluster ",i,": "))
         todo[todo$old_names==i,"new_names"]=x
      }
      if(!is.null(refine)){
         hand_names=rbind(todo,tokeep)
         hand_names %>%
            arrange(desc(old_names))
      }else{hand_names=todo}
      return(hand_names)


   }else if(data_type=="AvgExp"){
      ann3=ann@results$marker_free$corr
      rownames(ann3) <- gsub("^query\\.", "", rownames(ann3))
      colnames(ann3) <- gsub("^db\\.", "", colnames(ann3))
      hand_names=data.frame(rownames(ann3))
      colnames(hand_names)[1]="old_names"
      hand_names$new_names=""
      if (!is.null(refine)) {
         todo <- hand_names %>%
            filter(old_names %in% unique(refine$cluster[refine$best_match=="INCONCLUSIVE"]))
         tokeep <- hand_names %>%
            filter(old_names %in% unique(refine$cluster[refine$best_match!="INCONCLUSIVE"]))
         tokeep$new_names=refine[refine$best_match!="INCONCLUSIVE","best_match"]
      }else{todo=hand_names}
      for (i in todo$old_names) {
         temp = ann3[i,]
         mat <- data.matrix(temp)
         print(Heatmap(mat,
                       name = "Pearson correlation",
                       cluster_columns = F,
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 11),
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 11),
                       row_title = "query cluster",
                       row_title_rot = 90,
                       row_title_side = "left",
                       row_title_gp = gpar(fontface = "bold"),
                       column_title = "db",
                       column_title_rot = 0,
                       column_title_side = "bottom",
                       column_title_gp = gpar(fontface = "bold")))
         x=readline(paste0("What would you like to name cluster ",i,": "))
         todo[todo$old_names==i,"new_names"]=x
      }
      if(!is.null(refine)){
         hand_names=rbind(todo,tokeep)
         hand_names %>%
            arrange(desc(old_names))
      }else{hand_names=todo}
      return(hand_names)



   }else if(data_type=="Both"){
      ann2=ann@results$marker_based$dotplot_all
      hand_names1=data.frame(unique(ann2$data$cluster))
      colnames(hand_names1)[1]="old_names"
      ann3=ann@results$marker_free$corr
      rownames(ann3) <- gsub("^query\\.", "", rownames(ann3))
      colnames(ann3) <- gsub("^db\\.", "", colnames(ann3))
      hand_names2=data.frame(rownames(ann3))
      colnames(hand_names2)[1]="old_names"
      hand_names=merge(hand_names1,hand_names2)
      hand_names <- hand_names %>% arrange(desc(old_names))
      hand_names$new_names=""
      if (!is.null(refine)) {
         todo <- hand_names %>%
            filter(old_names %in% unique(refine$cluster[refine$best_match=="INCONCLUSIVE"]))
         tokeep <- hand_names %>%
            filter(old_names %in% unique(refine$cluster[refine$best_match!="INCONCLUSIVE"]))
         tokeep$new_names=refine[refine$best_match!="INCONCLUSIVE","best_match"]
      }else{todo=hand_names}
      for (i in todo$old_names) {
         temp2=ann2
         temp2$data=ann2$data[ann2$data$cluster==i,]
         p1=temp2+theme(legend.position = "none")
         temp = ann3[i,]
         mat <- data.matrix(temp)
         p2=Heatmap(mat,
                    name = "Pearson correlation",
                    cluster_columns = F,
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 11),
                    column_names_side = "bottom",
                    column_names_gp = gpar(fontsize = 11),
                    row_title = "query cluster",
                    row_title_rot = 90,
                    row_title_side = "left",
                    row_title_gp = gpar(fontface = "bold"),
                    column_title = "db",
                    column_title_rot = 0,
                    column_title_side = "bottom",
                    column_title_gp = gpar(fontface = "bold"))
         p2_grob <- grid.grabExpr(draw(p2))
         print(ggarrange(p1,p2_grob,nrow = 2))
         x=readline(paste0("What would you like to name cluster ",i,": "))
         todo[todo$old_names==i,"new_names"]=x
      }
      if(!is.null(refine)){
         hand_names=rbind(todo,tokeep)
         hand_names %>%
            arrange(desc(old_names))
      }else{hand_names=todo}
      return(hand_names)

   }

}
