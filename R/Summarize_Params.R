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

   init_params <- NULL
   init_selfsim_data <- NULL
   self_similarity_params <- NULL
   ds_params <- NULL


   if (!is.null(ann@params$markers$p_thresh)) {
      # Initialize_Markers parameters with Analysis = "Marker-Based"
      init_params <- data.frame(
         Analysis = "Marker-Based",
         Input = "Query",
         Function = "Initialize_Markers",
         Parameter = c("p_thresh", "FC_thresh", "sens_thresh", "spec_thresh"),
         Value = c(ann@params$markers$p_thresh, ann@params$markers$FC_thresh,
                   ann@params$markers$sens_thresh, ann@params$markers$spec_thresh)
      )
   }

   # Initialize_SelfSimilarity parameters
   if (!is.null(ann@params$markers$selfsim_analysis)) {
      for (i in seq_along(ann@params$markers$selfsim_analysis)) {
         sim_params <- ann@params$markers$selfsim_analysis[[i]]
         # Determine the prefix based on Input type
         prefix <- "Query_" # for future use: #ifelse(sim_params$set == "db", "db_", "query_")

         method_params <- data.frame(
            Analysis = "Self_Sim",
            Input = "Query", # for future use: #ifelse(sim_params$set == "db", "db_", "query_")
            Function = "Initialize_SelfSimilarity",
            Parameter = paste0(prefix, c("slot")),
            Value = c(sim_params$slot)
         )
         # Append each set of parameters for Tune_Markers
         init_selfsim_data <- rbind(init_selfsim_data, method_params)
      }
   }

   if (!is.null(ann@params$markers$assay_markers)) {
      # Self_Similarity parameters with Analysis = "Self_Sim"
      self_similarity_params <- data.frame(
         Analysis = "Self_Sim",
         Input = "Query",
         Function = "Create_SelfSimilarity_Viz",
         Parameter = "assay_markers",
         Value = ann@params$markers$assay_markers
      )
   }

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


   if (!is.null(ann@params$marker_free$downsample) &&
       !is.null(ann@params$marker_free$length_ds) &&
       !is.null(ann@params$marker_free$norm_method) &&
       !is.null(ann@params$marker_free$corr_method)
						) {
      # Downsample parameters
      ds_params <- data.frame(
         Analysis = "Marker-Free",
         Input = "Query&DB",
         Function = "Downsample",
         Parameter = c("downsample", "length_ds","norm_method","corr_method"),
         Value = c(ann@params$marker_free$downsample, ann@params$marker_free$length_ds,ann@params$marker_free$norm_method,ann@params$marker_free$corr_method)
      )
   }

   # Combine all parameters

   # Combine only the non-null data frames
   params_list <- list(params_data, init_params, init_selfsim_data, self_similarity_params, ds_params)
   params_data <- do.call(rbind, Filter(Negate(is.null), params_list))

   # Clear and add the params to ann@results
   ann@params$summary <- NULL
   ann@params$summary <- params_data

   return(ann)
}
