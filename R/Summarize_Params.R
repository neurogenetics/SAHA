#' Collects parameter settings for SAHA analysis
#'
#' This function gathers parameter settings for various analyses (e.g., marker-based, self-similarity)
#' and stores them in a structured data frame within the `ann` object.
#'
#' @param ann An object containing analysis parameters within `ann@params`.
#'
#' @return The input `ann` object, updated with a data frame of parameters in `ann@results$params`.
#'
#' @importFrom dplyr rbind
#'
#' @export
Summarize_Params <- function(ann) {
   # Initialize an empty data frame for params
   params_data <- data.frame(Analysis = character(), Input = character(), Function = character(), Parameter = character(), Value = character())

   # Initialize_Markers parameters with Analysis = "Marker-Based"
   init_params <- data.frame(
      Analysis = "Marker-Based",
      Input = NA,  # Not specific to Query/DB
      Function = "Initialize_Markers",
      Parameter = c("p_thresh", "FC_thresh", "sens_thresh", "spec_thresh"),
      Value = c(ann@params$markers$p_thresh, ann@params$markers$FC_thresh,
                ann@params$markers$sens_thresh, ann@params$markers$spec_thresh)
   )

   # Self_Similarity parameters with Analysis = "Self_Sim"
   self_similarity_params <- data.frame(
      Analysis = "Self_Sim",
      Input = NA,
      Function = "Create_SelfSimilarity_Viz",
      Parameter = "assay_markers",
      Value = ann@params$markers$assay_markers
   )

   # Tune_Markers parameters with conditional naming for DB and Analysis = "Marker-Based"
   if (!is.null(ann@params$markers$tuning_history)) {
      for (i in seq_along(ann@params$markers$tuning_history)) {
         tune_params <- ann@params$markers$tuning_history[[i]]
         # Determine the prefix based on Input type
         prefix <- ifelse(tune_params$set == "db", "db_", "query_")

         method_params <- data.frame(
            Analysis = "Marker-Based",
            Input = ifelse(tune_params$set == "query", "Query", "DB"),
            Function = "Tune_Markers",
            Parameter = paste0(prefix, c("set", "method", "method_var", "method_value")),
            Value = c(tune_params$set, tune_params$method, tune_params$method_var, tune_params$method_value)
         )
         # Append each set of parameters for Tune_Markers
         params_data <- rbind(params_data, method_params)
      }
   }

   # Combine all parameters
   params_data <- rbind(params_data, init_params, self_similarity_params)

   # Clear and add the params to ann@results
   ann@params$summary <- NULL
   ann@params$summary <- params_data

   return(ann)
}
