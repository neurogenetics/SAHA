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

Investigate_Self_Similarity <- function(df, cluster1, cluster2,shared_df=NULL) {
   # Ensure data is numeric
   df <- as.matrix(df)
   df[is.na(df)] <- 0

   # Extract rows for the specified clusters
   cluster1_data <- df[paste(cluster1), ]
   cluster2_data <- df[paste(cluster2), ]

   # Find shared features (where both rows have 1)
   bound=cbind(cluster1_data,cluster2_data)
   bound=data.frame(bound)
   venn_input_1=rownames(bound[bound$cluster1_data==1,])
   venn_input_2=rownames(bound[bound$cluster2_data==1,])

   # Create a list of sets
   sets <- list(`Cluster_1` = venn_input_1, `Cluster_2` = venn_input_2)

   # Generate the Euler diagram
   euler_venn <- euler(sets)

   # Open a new plotting window
   plot.new()

   # Adjust plot aesthetics
   print(plot(euler_venn,
              labels = c(paste(cluster1), cluster2),
              quantities = TRUE,
              font.main = 14,
              font.sub = 12,
              shape = "circle",
              fill = list(fill = c("pink", "gray"))))
   if (!is.null(shared_df)) {
      return(bound)
   }

}
