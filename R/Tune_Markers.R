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

Tune_Markers = function(ann, method, method_value, method_var,set){
   if (is.null(ann@ann1[[set]])) {
      print("WARNING: Markers for the given set have not been initialized. Please run Initialize_Markers() prior to tuning with this funciton.")
   }else{
      if (method == "absolute") {
         tune_subset <- ann@ann1[[set]] %>%
            group_by(cluster)%>%
            top_n(method_value, wt = !!sym(method_var))%>%
            ungroup()

      }else if (method == "relative") {
         tune_subset <- ann@ann1[[set]] %>%
            group_by(cluster) %>%
            filter(!!sym(method_var) >= quantile(!!sym(method_var), probs = method_value)) %>%
            ungroup()
      }else{
         print("Method not supported. Please check documentation ?Tune_Markers().")
      }
   }
 #  input_summary <- ann@ann1[[set]]%>%
 #     group_by(cluster) %>%
 #     summarise(n_markers =n())

#   output_summary <- data.frame(tune_subset)%>%
#      group_by(cluster) %>%
#      summarise(n_markers =n())
#   print(paste0("By running Tune_Markers() you have lowered the median markers per cluster from ", summary(input_summary$n_markers)[3], " to ",summary(output_summary$n_markers)[3],"."))
   ann@ann1[[set]]=data.frame(tune_subset)
   return(ann)


}
