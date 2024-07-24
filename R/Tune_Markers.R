# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

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
   input_summary <- ann@ann1[[set]]%>%
      group_by(cluster) %>%
      summarise(n_markers =n())

   output_summary <- data.frame(tune_subset)%>%
      group_by(cluster) %>%
      summarise(n_markers =n())
   print(paste0("By running Tune_Markers() you have lowered the median markers per cluster from ", summary(input_summary$n_markers)[3], " to ",summary(output_summary$n_markers)[3],"."))
   ann@ann1[[set]]=data.frame(tune_subset)
   return(ann)


}
