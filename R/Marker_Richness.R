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

MarkerRichness = function(ann, varfeat = NULL){
   if (is.null(varfeat)) {
      ann1=ann@ann1
      df = data.frame(ann1$query[,c("gene","cluster")])

      gene_count_by_cluster <- df %>%
         group_by(cluster) %>%  # Filter for rows where gene and cluster values match
         count()  # Count the number of rows after filtering

      print(gene_count_by_cluster)
   }else{
      ann1=ann@ann1
      df = data.frame(ann1$query[,c("gene","cluster")])

      # Count genes present in varfeat
      gene_counts <- df %>%
         mutate(in_varfeat = gene %in% varfeat) %>%  # Create a new column indicating presence
         group_by(cluster)%>%
         summarize(
            all_genes = n(),
            genes_in_varfeat=sum(in_varfeat)  # Count occurrences of each gene and presence flag
         )
      # Print the results
      print(gene_counts)
   }
}
