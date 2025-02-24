#' Computes correlations between query and database expression profiles
#'
#' This function calculates correlation coefficients between all query and
#' database cell type expression profiles stored in the `ann@results$marker_free$norm_merge` data
#' frame within the `ann` object. The results are stored in a data frame named
#' `corr` within `ann@results$marker_free`.
#'
#' @param ann An object containing normalized gene expression data (`ann@results$marker_free$norm_merge`).
#' @param corr_method Character string indicating which correlation coefficient is to be computed.
#'        Options are "pearson" (default), "kendall", or "spearman".
#'
#' @return The modified `ann` object with the correlation data frame stored within `ann@results$marker_free$corr`.
#'
#' @importFrom stats cor.test
#'
#' @export
CorrelateDS <- function(ann, corr_method = "pearson"){
   # compute Pearson coefficients b/w query and ABC celltype expression profiles
   db_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("db", colnames(ann@results$marker_free$norm_merge))]
   query_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("query", colnames(ann@results$marker_free$norm_merge))]

   correlation.df <- matrix(NA, nrow = length(query_columns), ncol = length(db_columns))
   rownames(correlation.df) <- query_columns
   colnames(correlation.df) <- db_columns

   for (i in seq_along(query_columns)) {
      for (j in seq_along(db_columns)) {
         suppressWarnings(cor_test <- cor.test(ann@results$marker_free$norm_merge[[query_columns[i]]], ann@results$marker_free$norm_merge[[db_columns[j]]], method = corr_method))
         correlation.df[i, j] <- cor_test$estimate
      }
   }

   correlation.df <- as.data.frame(correlation.df)
   ann@results$marker_free$corr=correlation.df
   ann@params$marker_free$corr_method <- corr_method
   return(ann)
}
