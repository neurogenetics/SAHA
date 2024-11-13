
#' Creates output dataframes from a Seurat object
#'
#' This function generates variable features, average expression, and/or marker genes
#' from a Seurat object, depending on the specified output option.
#'
#' @param obj A Seurat object.
#' @param output A character string specifying the desired output. Options are "Both" (default),
#'   "AvgExp", or "Markers".
#' @details
#'   - If output is "Both", the function returns variable features, average expression, and markers.
#'   - If output is "AvgExp", it returns variable features and average expression.
#'   - If output is "Markers", it returns variable features and markers.
#'
#' @export
#' @importFrom Seurat VariableFeatures AverageExpression FindAllMarkers
Seurat2SAHA <- function(obj, output="Both"){
   if (output=="Both") {
     varfeat = VariableFeatures(obj)
     avgexp =  suppressMessages(AverageExpression(neuron,assays = "RNA",slot="data"))
     markers = FindAllMarkers(obj, only.pos = TRUE)
     output = list(varfeat=varfeat, avgexp=avgexp, markers=markers)
   }else if (output == "AvgExp") {
      varfeat = VariableFeatures(obj)
      avgexp =  suppressMessages(AverageExpression(neuron,assays = "RNA",slot="data"))
      output = list(varfeat=varfeat, avgexp=avgexp)
   }else if (output == "Markers") {
      varfeat = VariableFeatures(obj)
      markers = FindAllMarkers(obj, only.pos = TRUE)
      output = list(varfeat=varfeat, markers=markers)
   }
}
