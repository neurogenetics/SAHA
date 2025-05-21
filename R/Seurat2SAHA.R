#' Creates output dataframes from a Seurat object
#'
#' This function generates variable features, average expression, and/or marker genes from a Seurat object, depending on the specified output option.
#'
#' @param obj A Seurat object.
#' 
#' @param output A character string specifying the desired output. Options are "Both" (default), "AvgExp", or "Markers".
#' @details
#'   - If output is "Both", the function returns variable features, average expression, and markers.
#'   - If output is "AvgExp", it returns variable features and average expression.
#'   - If output is "Markers", it returns variable features and markers.
#'
#' @export
#' @importFrom Seurat VariableFeatures AverageExpression FindAllMarkers
Seurat2SAHA <- function(obj, output="Both"){
   if (output=="Both") {
     cat("Calculating Variable Features: \n")
     varfeat = VariableFeatures(obj)
     cat(paste0("Storing ",length(varfeat)," variable features. \n"))
     cat("Calculating Average Expression: \n")
     avgexp =  suppressMessages(AverageExpression(obj,assays = "RNA",slot="data"))
     avgexp=data.frame(avgexp)
     cat(paste0("Storing expression of ",length(rownames(avgexp))," unique features across ",length(colnames(avgexp))," clusters. \n"))
     cat("Calculating Cluster Markers: \n")
     markers = FindAllMarkers(obj, only.pos = TRUE)
     cat(paste0("Storing "),length(unique(markers$gene)), " unique markers across ", length(unique(markers$cluster)), " clusters.")
     output = list(varfeat=varfeat, avgexp=avgexp, markers=markers)
   }else if (output == "AvgExp") {
      cat("Calculating Variable Features: \n")
      varfeat = VariableFeatures(obj)
      cat(paste0("Storing ",length(varfeat)," variable features. \n"))
      cat("Calculating Average Expression: \n")
      avgexp =  suppressMessages(AverageExpression(obj,assays = "RNA",slot="data"))
      avgexp=data.frame(avgexp)
      cat(paste0("Storing expression of ",length(rownames(avgexp))," unique features across ",length(colnames(avgexp))," clusters. \n"))
      output = list(varfeat=varfeat, avgexp=avgexp)
   }else if (output == "Markers") {
      cat("Calculating Variable Features: \n")
      varfeat = VariableFeatures(obj)
      cat(paste0("Storing ",length(varfeat)," variable features. \n"))
      cat("Calculating Cluster Markers: \n")
      markers = FindAllMarkers(obj, only.pos = TRUE)
      cat(paste0("Storing "),length(unique(markers$gene)), " unique markers across ", length(unique(markers$cluster)), " clusters.")
      output = list(varfeat=varfeat, markers=markers)
   }
   return(output)
}
