#' Computes correlations between query and database expression profiles
#'
#' This function calculates correlation coefficients and p-values between all query and
#' database cell type expression profiles. Results are stored in `corr` and `pvalue` slots.
#'
#' @param ann An object containing normalized gene expression data (`ann@results$marker_free$norm_merge`).
#' 
#' @param corr_method Character string indicating which correlation coefficient is to be computed.
#'        Options are "pearson" (default), "kendall", or "spearman".
#'
#' @return The modified `ann` object with correlation and p-value data frames.
#'
#' @importFrom stats cor.test
#'
#' @export
CorrelateDS <- function(ann, corr_method = "pearson"){
   # Identify columns for database and query
   db_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("db", colnames(ann@results$marker_free$norm_merge))]
   query_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("query", colnames(ann@results$marker_free$norm_merge))]
   
   # Initialize matrices for estimates and p-values
   correlation.mat <- matrix(NA, nrow = length(query_columns), ncol = length(db_columns))
   pvalue.mat <- matrix(NA, nrow = length(query_columns), ncol = length(db_columns))
   
   rownames(correlation.mat) <- rownames(pvalue.mat) <- query_columns
   colnames(correlation.mat) <- colnames(pvalue.mat) <- db_columns
   
   # Loop through each pair to run cor.test
   for (i in seq_along(query_columns)) {
      for (j in seq_along(db_columns)) {
         # Run test and suppress warnings (e.g., for constant vectors)
         suppressWarnings(cor_test <- stats::cor.test(
            ann@results$marker_free$norm_merge[[query_columns[i]]], 
            ann@results$marker_free$norm_merge[[db_columns[j]]], 
            method = corr_method
         ))
         
         # Save estimate and p-value
         correlation.mat[i, j] <- cor_test$estimate
         pvalue.mat[i, j] <- cor_test$p.value
      }
   }
   
   # Store results in the ann object
   ann@results$marker_free$corr <- as.data.frame(correlation.mat)
   ann@results$marker_free$pvalue <- as.data.frame(pvalue.mat) # New Slot!
   
   ann@params$marker_free$corr_method <- corr_method
   
   return(ann)
}
