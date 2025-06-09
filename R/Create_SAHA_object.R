#' Creation of SAHA object for analysis
#'
#' The `SAHA` class is used to store and manage data throughout the SAHA analysis workflow. It holds information about query and database data, matched markers for different analyses, and final analysis results.
#'
#' @slot query A list containing data frames for both marker (Markers) and average expression (AvgExp) data of the query.
#' 
#' @slot db A list containing data frames for both marker (Markers) and average expression (AvgExp) data of the database.
#' 
#' @slot ann1 A list to store data frames containing matched query and database markers for self-similarity analyses (populated later).
#' 
#' @slot ann2 A list to store data frames containing matched query and database markers for marker-based analyses (populated later).
#' 
#' @slot ann3 A list to store data frames containing matched query and database markers for marker-free analyses (populated later).
#' 
#' @slot results A list to store objects containing final analysis results required for visualizations (populated later).
#' 
#' @slot data_type A character string indicating the currently loaded data type ("Markers", "AvgExp", or "Markers & AvgExp").
#'
#' @export

Create_SAHA_object <- function(query, db,data_type,existing=NULL){
   #### IF statement that allows users to either load Marker or AvgExp data
   if (data_type == "Markers") {
      #loading query data
      cat("You have selected to load Query and DB marker data.\n")
      #1. Read in your queried data as a df or csv
      if (inherits(query, "data.frame")) {
         # If query is already a data frame, use it directly
         query_marker_data <- query
      } else {
         # If query is a path to a CSV file, read it into a data frame
         query_marker_data <- read.csv(query)
      }


      #2. Read in your db data as a df or csv
      if (inherits(db, "data.frame")) {
         # If query is already a data frame, use it directly
         db_marker_data <- db
      } else {
         # If query is a path to a CSV file, read it into a data frame
         db_marker_data <- read.csv(db)
      }
   }else{
      cat("You have selected to load Query and DB AvgExp data.\n")
      #1. Read in your queried data as a df or csv
      if (inherits(query, "data.frame")) {
         # If query is already a data frame, use it directly
         query_AvgExp_data <- query
      } else {
         # If query is a path to a CSV file, read it into a data frame
         query_AvgExp_data <- read.csv(query)
      }

      #2. Read in your db data as a df or csv
      if (inherits(db, "data.frame")) {
         # If query is already a data frame, use it directly
         db_AvgExp_data <- db
      } else {
         # If query is a path to a CSV file, read it into a data frame
         db_AvgExp_data <- read.csv(db)
      }
   }
   if (data_type == "Markers") {
      temp_ann <- new("SAHA", data_type=data_type, query = list(Markers=data.frame(query_marker_data), AvgExp=data.frame(NULL)),db = list(Markers=data.frame(db_marker_data), AvgExp=data.frame(NULL)),params=list(NULL))
   }else if (data_type =="AvgExp"){
      temp_ann <- new("SAHA", data_type=data_type, query = list(Markers=data.frame(NULL), AvgExp=data.frame(query_AvgExp_data)),db = list(Markers=data.frame(NULL), AvgExp=data.frame(db_AvgExp_data)),params=list(NULL))
   }

   if (data_type == "Markers") {
      #3A. Check to see that all rownames of query are in rownames of db
      cat(paste("Loaded database contains", length(unique(temp_ann@db$Markers$cluster)),"unique cell types.\n"))
      #print the name of the loaded dataset
      cat(paste("Loaded query dataset contains", length(unique(temp_ann@query$Markers$cluster)),"unique clusters.\n"))
      matched_genes <- intersect(temp_ann@query$Markers$gene, temp_ann@db$Markers$gene)
      if (length(matched_genes) == 0) {
         warning("\nSAHA did not detect any matching genes between your query and your database of interest.\nYour data has still been processed.\nPlease refer to the troubleshooting vignette.\n")
      }
   }else if (data_type =="AvgExp"){
      #3B. Check to see that all rownames of query are in rownames of db
      cat(paste("Loaded database contains", length(colnames(temp_ann@db$AvgExp)),"unique cell types.\n"))
      #print the name of the loaded dataset
      cat(paste("Loaded query dataset contains", length(colnames(temp_ann@query$AvgExp)),"unique clusters.\n"))
      }
}

   if (is.null(existing)) {
      ann=temp_ann
      return(ann)
   }else{
      # Check if 'existing' is a valid SAHA object
      if (!inherits(existing, "SAHA")) {
         warning("\n'existing' object is not a valid SAHA object. Proceed with a new SAHA object.\n")
         ann = temp_ann
         return(ann)
      }

      if (data_type=="Markers") {
         existing@query$Markers = data.frame(temp_ann@query$Markers)
         existing@db$Markers = data.frame(temp_ann@db$Markers)
         existing@data_type = "AvgExp & Markers"
      }else if(data_type=="AvgExp"){
         existing@query$AvgExp = data.frame(temp_ann@query$AvgExp)
         existing@db$AvgExp = data.frame(temp_ann@db$AvgExp)
         existing@data_type = "AvgExp & Markers"
      }
      return(existing)
   }

}
