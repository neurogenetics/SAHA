#' Looks up marker overlap between query cluster and database cell types
#'
#' This function calculates the number of markers from a specified query cluster that overlap with each cell type in the database.
#'
#' @param cluster The name of the query cluster.
#' 
#' @param ann An object containing marker data in `ann1$query` and `ann1$db`.
#'
#' @return A numeric vector representing the number of markers from the query cluster that overlap with each database cell type.
#'
#' @importFrom dplyr filter
#'
#' @export
SAHA_lookup_cluster <- function(cluster,ann) {
   k=ann@ann1$query[ann@ann1$query$cluster==cluster,"gene"]
   return(summary(factor(ann@ann1$db[ann@ann1$db$gene %in% k, "cluster"], levels(factor(ann@ann1$db$cluster)))))
}
