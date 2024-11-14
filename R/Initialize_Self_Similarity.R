#' Initialize Self-Similarity Matrices
#'
#' This function initializes self-similarity matrices based on the specified slot.
#'
#' @param ann An object containing annotation data.
#' @param slot A character string specifying the slot. It can be "AvgExp" or "Markers".
#' @return An updated object with the self-similarity matrix stored in the appropriate slot.
#' @importFrom dplyr %>%
#' @importFrom tidyr spread
#' @importFrom tibble column_to_rownames
#' @export

Initialize_Self_Similiarity <- function(ann,slot,custom_ds=NULL){
   if (slot=="AvgExp") {
      temp=ann@query$AvgExp
      if ("X"%in%colnames(temp)) {
         temp <- temp %>%
            column_to_rownames("X")
      }
      if (!is.null(custom_ds)) {
         temp <- temp[rownames(temp)%in%custom_ds,]
      }
      cor_mat <- cor(temp)
      # Print the resulting dataframe
      ann@ann1$AvgExp_selfsim_matrix=cor_mat
   }else if(slot=="Markers"){
      temp=ann@query$Markers
      df = data.frame(temp[,c("gene","cluster")])
      df$incidence = 1

      # Reshape using spread (assuming values for each id and time combination)
      wide_data <- df %>%
         spread(gene, incidence)  # Spreads "value" based on "time" for each "id"

      # Print the resulting wide dataframe
      wide_data <- wide_data %>%
         column_to_rownames("cluster")

      # Print the resulting dataframe
      ann@ann1$Marker_selfsim_matrix=wide_data
   }

   # Save the method parameters as a tuning entry in the history
   selfsim_entry <- list(
      slot = slot
   )


   # Append to tuning history, creating it if it doesn't exist
   if (is.null(ann@params$markers$selfsim_analysis)) {
      ann@params$markers$selfsim_analysis <- list(selfsim_entry)
   } else {
      ann@params$markers$selfsim_analysis <- append(ann@params$markers$selfsim_analysis, list(selfsim_entry))
   }
   return(ann)  # Return the updated object at the end
}
