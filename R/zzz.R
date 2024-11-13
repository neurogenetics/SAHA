# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
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

.onAttach <- function(...) {
   packageStartupMessage("Thank you for using SAHA. You have downloaded version prior to release (V1.0), as such it is in beta and may not work as intended. Check back at www.github.com/neurogenetics/SAHA/ periodically for updates. If you choose to use this in any preprints or publication, please cite using the zenodo DOI linked to this repo: 10.5281/zenodo.14040698.")
}

.onLoad <- function(...){
   setClass("SAHA",
            slots = c(
               data_type = "character",
               query = "list",
               db = "list",
               ann1 = "list",
               ann2 = "list",
               ann3 = "list",
               results = "list",
               params = "list"
            )
   )
   setMethod("show", "SAHA", function(object) {

      # Print the class type and data_type
      cat(is(object)[[1]], "\n",
          " Data Type:", object@data_type, "\n")

      # Avg exp data
      if (object@data_type == "AvgExp") {

         # conditionally store information

         #queried avg exp
         if (!is.null(object@query$AvgExp$X)) {
            q_ae <- length(object@query$AvgExp$X) # total gene entries
            q_clusters <- length(colnames(object@query$AvgExp))-1 # q clusters
         } else {
            q_ae <- length(rownames(object@query$AvgExp)) # total gene entries
            q_clusters <- length(colnames(object@query$AvgExp))-1 # q clusters
         }

         #db avg exp
         if (length(object@db$AvgExp) > 0) {
            db_ae <- length(rownames(object@db$AvgExp)) # db total gene entries
            db_anno <- length(colnames(object@db$AvgExp))-3 # db possible annotations
         }

         # Print results
         # reformat: 8533 query markers (31 clusters), 18971 database markers (19 possible annotations)
         cat(sprintf("  Description: %d avgexp query (%d clusters)\n               %d avgexp db (%d possible annotations)\n",
                     q_ae, q_clusters, db_ae, db_anno))

      }

      # Markers data
      if (object@data_type == "Markers") {
         #queried marker
         if (length(unique(object@query$Markers$cluster)) > 0) {
            q_clusters <- length(unique(object@query$Markers$cluster))  #queried clusters
            q_markers <- length(unique(object@query$Markers$gene)) #queried gene entries
         }

         #db marker
         if (length(unique(object@db$Markers$cluster)) > 0) {
            db_anno <- length(unique(object@db$Markers$cluster)) #database # of cell types
            db_markers <- length(unique(object@db$Markers$SYMBOL))#database gene entries
         }

         #print results
         #246 avgexp query (32 clusters), 32285 avgexp db (22 possible annotations)
         cat(sprintf("  Description: %d query markers (%d clusters)\n               %d database markers (%d possible annotations)\n",
                     q_markers, q_clusters, db_markers, db_anno))
      }

      # Both data
      if (object@data_type %in% c("Markers & AvgExp", "AvgExp & Markers")) {
         #queried avg exp
         if (!is.null(object@query$AvgExp$X)) {
            q_ae_gene <- length(object@query$AvgExp$X) #AvgExp gene entries
            q_ae_clusters <- length(colnames(object@query$AvgExp)) #AvgExp clusters
         } else {
            q_ae_gene <- length(rownames(object@query$AvgExp)) #AvgExp gene entries
            q_ae_clusters <- length(colnames(object@query$AvgExp)) #AvgExp gene clusters
         }

         #db avg exp
         if (length(object@query$AvgExp) > 0) {
            db_ae <- length(rownames(object@db$AvgExp)) # db total gene entries
            db_anno <- length(colnames(object@db$AvgExp)) # db possible annotations
         }

         #queried marker
         if (length(unique(object@query$Markers$cluster)) > 0) {
            q_db_clusters <- length(unique(object@query$Markers$cluster))  #queried clusters
            q_db_markers <- length(unique(object@query$Markers$gene)) #queried gene entries
         }

         #db marker
         if (length(unique(object@db$Markers$cluster)) > 0) {
            db_anno <- length(unique(object@db$Markers$cluster)) #database # of cell types
            db_markers <- length(unique(object@db$Markers$SYMBOL))#database gene entries
         }
         # put the results in two different lines
         cat(sprintf(
            "  Description: %d avgexp query (%d clusters)\n               %d avgexp db (%d possible annotations)\n               %d query markers (%d clusters)\n               %d database markers (%d possible annotations)\n",
            q_ae_gene, q_ae_clusters, db_ae, db_anno,
            q_db_markers, q_db_clusters, db_markers, db_anno
         ))
      }
      # Check which annotations are present and print them
      ann_list <- c()

      if (length(object@ann1)!= 0) {
         ann_list <- c(ann_list, "self-similarity")
      }
      if (length(object@ann2)!= 0) {
         ann_list <- c(ann_list, "marker-based")
      }
      if (length(object@ann3)!= 0) {
         ann_list <- c(ann_list, "marker-free")
      }
      if (length(ann_list) > 0) {
         cat("  Analysis:", paste(ann_list, collapse = ", "), "\n")
      }else {
         cat("  Analysis: none")
      }
   }
   )
}

## Make later to clean up anything if needed
#.onUnload
