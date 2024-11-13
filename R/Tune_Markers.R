#' Tunes marker set based on specified criteria
#'
#' This function filters markers within a specified set (query or database) based on
#' absolute or relative thresholding of a given variable.
#'
#' @param ann An object containing marker data in `ann1$query` or `ann1$db`.
#' @param method The tuning method ("absolute" or "relative").
#' @param method_value The threshold value for the tuning method.
#' @param method_var The variable to apply the threshold to.
#' @param set The set to tune ("query" or "db").
#'
#' @return The input `ann` object with the tuned marker set in `ann1[[set]]`.
#'
#' @importFrom dplyr %>% group_by, top_n, filter, ungroup, sym
#'
#' @export
Tune_Markers = function(ann, method, method_value, method_var, set) {
   if (is.null(ann@ann1[[set]])) { #and statement here, if method %in% colnames(ann@ann1[[set]])
      warning("Either Initialize_Markers() wasn't run or the specificed method does not exist for the set. Please refer to documentation.")
   } else {
      if (method == "absolute") {
         tune_subset <- ann@ann1[[set]] %>%
            group_by(cluster) %>%
            top_n(method_value, wt = !!sym(method_var)) %>%
            ungroup()

      } else if (method == "relative") {
         tune_subset <- ann@ann1[[set]] %>%
            group_by(cluster) %>%
            filter(!!sym(method_var) >= quantile(!!sym(method_var), probs = method_value)) %>%
            ungroup()

      } else {
         print("Method not supported. Please check documentation ?Tune_Markers().")
         return(ann)
      }

      # Save the method parameters as a tuning entry in the history
      tuning_entry <- list(
         method = method,
         method_value = method_value,
         method_var = method_var,
         set = set
      )

      # Append to tuning history, creating it if it doesn't exist
      if (is.null(ann@params$markers$tuning_history)) {
         ann@params$markers$tuning_history <- list(tuning_entry)
      } else {
         ann@params$markers$tuning_history <- append(ann@params$markers$tuning_history, list(tuning_entry))
      }

      # Assign tuned markers to the object
      ann@ann1[[set]] <- data.frame(tune_subset)
   }

   # Print summary information
   input_summary <- aggregate(gene ~ cluster, slot(ann,set)$Markers, length)
   output_summary <- aggregate(gene ~ cluster, ann@ann1[[set]], length)
   cat(paste0("By running Tune_Markers() you have lowered the median markers per cluster from ",
              summary(input_summary$gene)[3], " to ", summary(output_summary$gene)[3], "."))

   return(ann)
}
