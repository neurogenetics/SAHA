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
