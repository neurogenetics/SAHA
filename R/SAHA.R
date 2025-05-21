#' Conducts SAHA analysis on query and database data
#'
#' This function orchestrates the SAHA analysis pipeline based on the specified `data_type`. It includes steps for marker-based or marker-free analysis, depending on the input data.
#'
#' @param query A data frame or path to a CSV file containing query data (markers or average expression).
#' 
#' @param db A data frame or path to a CSV file containing database data (markers or average expression).
#' 
#' @param meta A data frame containing cell type metadata (required for marker-based analysis).
#' 
#' @param data_type The type of data being analyzed ("Markers" or "AvgExp").
#'
#' @return The function calls appropriate plotting functions based on the analysis type and returns the corresponding plots.
#'
#' @importFrom utils read.csv
#'
#' @export

SAHA <- function(query, db, meta, data_type) {
   capture.output({
      ann <- suppressMessages(Create_SAHA_object(query = query, db = db, data_type = data_type))

      if (data_type == "Markers") {
         ann <- suppressMessages(Initialize_Markers(ann))
         ann <- suppressMessages(Tune_Markers(ann = ann, method = "absolute", method_value = 100, method_var = "avg_log2FC", set = "db"))
         ann <- suppressMessages(Tune_Markers(ann = ann, method = "relative", method_value = 0.75, method_var = "avg_log2FC", set = "query"))
         ann <- suppressMessages(Run_Marker_Based(ann))
         ann <- suppressMessages(Create_MarkerBased_Viz(ann, meta = meta))
         return(call_SAHA_plots(ann, plot_type = "Marker-based", data_type = "Markers"))

      } else if (data_type == "AvgExp") {
         ann <- suppressMessages(Initialize_MarkerFree(ann = ann))
         ann <- suppressMessages(Downsample(ann))
         ann <- suppressMessages(NormalizeDS(ann, assay_query = "RNA"))
         ann <- suppressMessages(CorrelateDS(ann))
         ann <- suppressMessages(Create_MarkerFree_Viz(ann, meta = meta))
         return(call_SAHA_plots(ann, plot_type = "Marker-free", data_type = "AvgExp"))

      } else {
         cat("Something went wrong. Please read the documentation or consider running the full SAHA pipeline.")
      }
   })
}
