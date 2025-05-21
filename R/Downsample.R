#' Downsamples query and database to shared genes
#'
#' This function reduces the size of the query and database components within an `ann` object by retaining only genes expressed in both datasets.
#'
#' @param ann An object containing `query$AvgExp` and `db$AvgExp` data frames.
#'
#' @return The input `ann` object with an added `ann3` component containing downsampled data frames for query and database.
#'
#' @export
Downsample <- function(ann, custom_ds=NULL){
   #3. Downsample database and query to mutually expressed genes
   if (is.null(custom_ds)) {
     query_genes=rownames(ann@query$AvgExp)
   }else{
     intended_query_genes=custom_ds
     overlap <- intersect(rownames(ann@db$AvgExp), rownames(ann@query$AvgExp))
     query_genes <- intersect(overlap, intended_query_genes)
   }
   db_ds=ann@db$AvgExp[rownames(ann@db$AvgExp)%in%query_genes,]
   query_ds=ann@query$AvgExp[rownames(ann@query$AvgExp)%in%rownames(db_ds),]

   #4. Create list of data.frames
   ann3=list("query"=query_ds, "db"=db_ds)
   #5. Check to see that all rownames of query are in rownames of db
   if (all(rownames(ann3$query)%in%rownames(ann3$db))) {
      cat(paste("Downsampled query and database contain",length(rownames(ann3$query)),"genes."))
   }else{
      cat("Something went wrong! It appears there are either no shared genes between query and db. Please check downsampling manually to ensure that symbols (gene names) are in the same format and that query and db share common genes for SAHA flavor 3.")
   }

   ann@params$marker_free <- list(downsample=ifelse(is.null(custom_ds), "query&db",deparse(substitute(custom_ds))),length_ds= ifelse(is.null(custom_ds), length(rownames(ann3$query)), paste0(length(rownames(ann3$query)),"/",length(custom_ds))))

   #6. Return the dataframe
   ann@ann3=ann3
   return(ann)
}
