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

SAHA <- function(query,db,meta,data_type){
   ann=Create_SAHA_object(query = query,db = db,data_type = data_type)
   if (data_type=="Markers") {
      ann=Initialize_Markers(ann)
      ann=Tune_Markers(ann = ann,method = "absolute",method_value = 100,method_var = "avg_log2FC",set = "db")
      ann=Tune_Markers(ann = ann,method = "relative",method_value = 0.75,method_var = "avg_log2FC",set = "query")
      ann=Run_Marker_Based(ann)
      ann=Create_MarkerBased_Viz(ann,meta = meta,facet = TRUE)
      return(call_SAHA_plots(ann, plot_type = "Marker-based",data_type = "Markers"))
   }else if (data_type=="AvgExp") {
      ann=Initialize_MarkerFree(ann = ann)
      ann=Downsample(ann)
      ann=NormalizeDS(ann,assay_query = "RNA")
      ann=CorrelateDS(ann)
      ann=Create_MarkerFree_Viz(ann,facet = TRUE,meta = meta, ABC = TRUE, chemistry = "10Xv3")
      return(call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp"))
   }else{
      print("Something went wrong. Please read the documentation or consider running the full SAHA pipeline.")
   }
}
