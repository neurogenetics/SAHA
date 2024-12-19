#' Initializes data for marker-free analysis
#'
#' This function performs pre-processing steps on the query and database components
#' within an `ann` object to prepare them for marker-free analysis in the SAHA workflow.
#' These steps include:
#'
#' * Removing potential column names that might interfere with analysis (`"X"`, `"X.1"`)
#' * Setting gene symbols as row names for the database data frame (if `SYMBOL` column exists)
#' * Removing duplicate genes and genes with missing symbols in the database data frame
#'
#' @param ann An object containing `query$AvgExp` and `db$AvgExp` data frames.
#'
#' @return The input `ann` object with pre-processed data frames for marker-free analysis.
#'   The function modifies the `query$AvgExp` and `db$AvgExp` slots within `ann`.
#'
#' @importFrom dplyr %>% select
#'
#' @export
Initialize_MarkerFree <- function(ann){
   #####CLEANING NEEDS TO BE PERFORMED PRIOR TO SAHA
   if ("X"%in%colnames(ann@query$AvgExp)) {
      rownames(ann@query$AvgExp)=ann@query$AvgExp$X
      ann@query$AvgExp = ann@query$AvgExp %>%
         select(-X)
   }
   if ("X.1"%in%colnames(ann@query$AvgExp)) {
      ann@query$AvgExp = ann@query$AvgExp %>%
         select(-X.1)
   }
   if ("X"%in%colnames(ann@db$AvgExp)) {
      ann@db$AvgExp = ann@db$AvgExp %>%
         select(-X)
   }
   if ("X.1"%in%colnames(ann@db$AvgExp)) {
      ann@db$AvgExp = ann@db$AvgExp %>%
         select(-X.1)
   }
   # if (all(!duplicated(ann@db$AvgExp$SYMBOL))) {
   #    rownames(ann@db$AvgExp)=ann@db$AvgExp$SYMBOL
   #    ann@db$AvgExp = ann@db$AvgExp %>%
   #       select(-SYMBOL, -gene)
   # }else{
   #    ann@db$AvgExp=ann@db$AvgExp[!duplicated(ann@db$AvgExp$SYMBOL),]
   #    ann@db$AvgExp=ann@db$AvgExp[!is.na(ann@db$AvgExp$SYMBOL),]
   #    rownames(ann@db$AvgExp)=ann@db$AvgExp$SYMBOL
   #    ann@db$AvgExp = ann@db$AvgExp %>%
   #       select(-SYMBOL, -gene)
   # }
   return(ann)
}

