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
Downsample <- function(ann){
   #3. Downsample database and query to mutually expressed genes
   query_genes=rownames(ann@query$AvgExp)
   db_ds=ann@db$AvgExp[rownames(ann@db$AvgExp)%in%query_genes,]
   query_ds=ann@query$AvgExp[rownames(ann@query$AvgExp)%in%rownames(db_ds),]
   #4. Create list of data.frames
   ann3=list("query"=query_ds, "db"=db_ds)
   #5. Check to see that all rownames of query are in rownames of db
   if (all(rownames(ann3$query)%in%rownames(ann3$db))) {
      print(paste("Downsampled query and database contain",length(rownames(ann3$query)),"genes."))
   }else{
      print("Something went wrong! It appears there are either no shared genes between query and db. Please check downsampling manually to ensure that symbols (gene names) are in the same format and that query and db share common genes for SAHA flavor 3.")
   }
   #6. Return the dataframe
   ann@ann3=ann3
   return(ann)
}
