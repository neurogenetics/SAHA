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
  if (missing(query)) {
    stop("Argument 'query' is missing, with no default value")
  }
  if (missing(db)) {
    stop("Argument 'db' is missing, with no default value")
  }
  if (missing(meta)) {
    stop("Argument 'meta' is missing, with no default value")
  }
  if (missing(data_type)) {
    stop("Argument 'data_type' is missing, with no default value")
  }
  suppressMessages({
    suppressWarnings({
      capture.output({
        
        if (data_type == "Markers") {
          ann <- Create_SAHA_object(query = query, db = db, data_type = "Markers")
          ann <- Initialize_Markers(ann)
          ann <- Tune_Markers(ann = ann, method = "absolute", method_value = 100, method_var = "avg_log2FC", set = "db")
          ann <- Tune_Markers(ann = ann, method = "relative", method_value = 0.75, method_var = "avg_log2FC", set = "query")
          ann <- Run_Marker_Based(ann)
          ann <- Create_MarkerBased_Viz(ann, meta = meta)
          return(call_SAHA_plots(ann, plot_type = "Marker-based", data_type = "Markers"))
        } else if (data_type == "AvgExp") {
          ann <- Create_SAHA_object(query = query, db = db, data_type = "AvgExp")
          ann <- Initialize_MarkerFree(ann = ann)
          ann <- Downsample(ann)
          ann <- NormalizeDS(ann, assay_query = "RNA")
          ann <- CorrelateDS(ann)
          ann <- Create_MarkerFree_Viz(ann, meta = meta)
          return(call_SAHA_plots(ann, plot_type = "Marker-free", data_type = "AvgExp"))
        } else if (data_type == "Both") {
          
          markers_ann <- Create_SAHA_object(query = query_Markers, db = ISOCTX_Markers, data_type = "Markers")
          markers_ann <- Initialize_Markers(markers_ann)
          markers_ann <- Tune_Markers(markers_ann, method = "absolute", method_value = 100, method_var = "avg_log2FC", set = "db")
          markers_ann <- Tune_Markers(markers_ann, method = "relative", method_value = 0.75, method_var = "avg_log2FC", set = "query")
          markers_ann <- Run_Marker_Based(markers_ann)
          markers_ann <- Create_MarkerBased_Viz(markers_ann, meta = meta)
          markers_result <- call_SAHA_plots(markers_ann, plot_type = "Marker-based", data_type = "Markers")
    
          avgexp_ann <- Create_SAHA_object(query = query_AvgExp, db = ISOCTX_AvgExp, data_type = "AvgExp")
          avgexp_ann <- Initialize_MarkerFree(avgexp_ann)
          avgexp_ann <- Downsample(avgexp_ann)
          avgexp_ann <- NormalizeDS(avgexp_ann, assay_query = "RNA")
          avgexp_ann <- CorrelateDS(avgexp_ann)
          avgexp_ann <- Create_MarkerFree_Viz(avgexp_ann, meta = meta)
          avgexp_result <- call_SAHA_plots(avgexp_ann, plot_type = "Marker-free", data_type = "AvgExp")
          
          return(invisible(list(markers = markers_result, avgexp = avgexp_result)))
        } else {
          stop(paste("Invalid data_type:", data_type, "- must be 'Markers', 'AvgExp', or 'Both'"))
        }
      }, type = "output")
    })
  })
}