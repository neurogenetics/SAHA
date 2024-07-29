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

SemiAutoAnnotate = function(ann,data_type=NULL){
   #prompt through each one
   if(is.null(data_type)){
      print("Thank you for choosing the semi-automated approach. If you wish to proceed please select which data_type to auto-annotate using: Markers, AvgExp, or Both")
   }else if(data_type=="Markers"){

      temp=ann@results$marker_based$all
      hand_names=data.frame(unique(temp$data$cluster))
      colnames(hand_names)[1]="old_names"
      hand_names$new_names=""
      for (i in hand_names$old_names) {
         temp2=temp
         temp2$data=temp$data[temp$data$cluster==i,]
         print(temp2+theme(legend.position = "none"))
         x=readline(paste0("What would you like to name cluster ",i,": "))
         hand_names[hand_names$old_names==i,"new_names"]=x
      }
      return(hand_names)


   }else if(data_type=="AvgExp"){
      ann3=ann@results$marker_free$corr
      rownames(ann3) <- gsub("^query\\.", "", rownames(ann3))
      colnames(ann3) <- gsub("^db\\.", "", colnames(ann3))
      hand_names=data.frame(rownames(ann3))
      colnames(hand_names)[1]="old_names"
      hand_names$new_names=""
      for (i in hand_names$old_names) {
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
         hand_names[hand_names$old_names==i,"new_names"]=x
      }
      return(hand_names)



   }else if(data_type=="Both"){
      ann2=ann@results$marker_based$all
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
      for (i in hand_names$old_names[1:2]) {
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
         hand_names[hand_names$old_names==i,"new_names"]=x
      }
      return(hand_names)

   }

}
