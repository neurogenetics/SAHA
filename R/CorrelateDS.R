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
CorrelateDS <- function(ann){
   # compute Pearson coefficients b/w query and ABC celltype expression profiles
   db_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("db", colnames(ann@results$marker_free$norm_merge))]
   query_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("query", colnames(ann@results$marker_free$norm_merge))]

   correlation.df <- matrix(NA, nrow = length(query_columns), ncol = length(db_columns))
   rownames(correlation.df) <- query_columns
   colnames(correlation.df) <- db_columns

   for (i in seq_along(query_columns)) {
      for (j in seq_along(db_columns)) {
         cor_test <- cor.test(ann@results$marker_free$norm_merge[[query_columns[i]]], ann@results$marker_free$norm_merge[[db_columns[j]]], method = "pearson")
         correlation.df[i, j] <- cor_test$estimate
      }
   }

   correlation.df <- as.data.frame(correlation.df)
   ann@results$marker_free$corr=correlation.df
   return(ann)
}
