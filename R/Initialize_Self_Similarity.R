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

Initialize_Self_Similiarity <- function(ann,slot){
   if (slot=="AvgExp") {
      temp=ann@query$AvgExp
      temp <- temp %>%
         column_to_rownames("X")
      cor_mat <- cor(temp)
      # Print the resulting dataframe
      ann@ann1$AvgExp_selfsim_matrix=cor_mat

      return(ann)
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
      return(ann)
   }

}
