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

Create_SAHA_object <- function(query, db,data_type,existing=NULL){
   #SAHA Object will be a class made up of 6 components
   # 1. query = a list of df's containing either Markers or AvgExp
   # 2. db = a list of df's containing either Markers or AvgExp
   # 3. ann1 = a list of df's containing matched query and db markers for self-similarity analyses
   # 4. ann2 = a list of df's containing matched query and db markers for marker-based analyses
   # 5. ann3 = a list of df's containing matched query and db markers for marker-free analyses
   # 6. results = a list of objects containing final analyses  required to make plots
   if (is.null(existing)){
      setClass("SAHA",
               slots = c(
                  data_type = "character",
                  query = "list",
                  db = "list",
                  ann1 = "list",
                  ann2 = "list",
                  ann3 = "list",
                  results = "list"
               )
      )
      setMethod("show", "SAHA", function(object) {
         cat(is(object)[[1]], "\n",
             "  data_type: ", object@data_type, "\n",
             "  Your query contains:  ", length(unique(object@query$Markers$cluster))," clusters.", "\n",
             "  Your database contains:  ", length(unique(object@db$Markers$cluster))," possible annotations.", "\n",
             sep = ""
         )
      })
   }else{print("Merging new markers into existing SAHA object.")}


   #################
   # Code written in Gemini for printing different slot types w/ lapply. I feel okay about it so work with it to adapt to printign diferent slot types.....
   #################
   #setMethod("set", "yourClassName", function(object, ..., verbose = TRUE) {
   # ... standard set method logic ...

   #  if (verbose) {
   #    slots <- slotNames(object)
   #    lapply(slots, function(slot) {
   #      value <- slot(object, slot)
   #      if (!is.null(value) && length(value) > 0) {
   #        cat(paste0("Slot ", slot, ":\n"))
   #
   #        # Determine slot type and apply appropriate summary function
   #        if (is.data.frame(value)) {
   #          print(head(value))  # Example: print first few rows
   #       } else if (is.list(value) && all(sapply(value, is.data.frame))) {
   # Handle list of data frames
   #          lapply(value, function(df) print(head(df)))
   #        } else if (is.character(value)) {
   #          cat(paste0("Value: ", value))
   #        } else {
   #          print(summary(value))  # Default summary for other types
   #        }
   #      }
   #    })
   #  }

   #  object
   #})

   #### IF statement that allows users to either load Marker or AvgExp data
   if (data_type == "Markers") {
      #loading query data
      print("You have selected to load Query and DB marker data.")
      #1. Read in your queried data as a df or csv
      if (inherits(query, "data.frame")) {
         # If query is already a data frame, use it directly
         query_marker_data <- query
      } else {
         # If query is a path to a CSV file, read it into a data frame
         query_marker_data <- read.csv(query)
      }

      #loading db data ########### THIS NEEDS TO BE A CHOICE TO ALLOW USERS TO LOAD PRE-LOADED PACKAGE DATA!!!!
      #2. Read in your db data as a df or csv
      if (inherits(db, "data.frame")) {
         # If query is already a data frame, use it directly
         db_marker_data <- db
      } else {
         # If query is a path to a CSV file, read it into a data frame
         db_marker_data <- read.csv(db)
      }
   }else{
      print("You have selected to load Query and DB AvgExp data.")
      #1. Read in your queried data as a df or csv
      if (inherits(query, "data.frame")) {
         # If query is already a data frame, use it directly
         query_AvgExp_data <- query
      } else {
         # If query is a path to a CSV file, read it into a data frame
         query_AvgExp_data <- read.csv(query)
      }

      #loading db data ########### THIS NEEDS TO BE A CHOICE TO ALLOW USERS TO LOAD PRE-LOADED PACKAGE DATA!!!!
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
      temp_ann <- new("SAHA", data_type=data_type, query = list(Markers=data.frame(query_marker_data), AvgExp=data.frame(NULL)),db = list(Markers=data.frame(db_marker_data), AvgExp=data.frame(NULL)))
   }else if (data_type =="AvgExp"){
      temp_ann <- new("SAHA", data_type=data_type, query = list(Markers=data.frame(NULL), AvgExp=data.frame(query_AvgExp_data)),db = list(Markers=data.frame(NULL), AvgExp=data.frame(db_AvgExp_data)))
   }

   if (data_type == "Markers") {
      #3A. Check to see that all rownames of query are in rownames of db
      print(paste("Loaded database contains", length(unique(temp_ann@db$Markers$cluster)),"unique cell types."))
      #print the name of the loaded dataset
      print(paste("Loaded query dataset contains", length(unique(temp_ann@query$Markers$cluster)),"unique clusters."))
   }else if (data_type =="AvgExp"){
      #3B. Check to see that all rownames of query are in rownames of db
      print(paste("Loaded database contains", length(colnames(temp_ann@db$AvgExp)),"unique cell types."))
      #print the name of the loaded dataset
      print(paste("Loaded query dataset contains", length(colnames(temp_ann@query$AvgExp)),"unique clusters."))
   }

   #the else could use an inherit that confirms that existing is a SAHA object....
   if (is.null(existing)) {
      ann=temp_ann
      return(ann)
   }else{
      if (data_type=="Markers") {
         existing@query$Markers = data.frame(temp_ann@query$Markers)
         existing@db$Markers = data.frame(temp_ann@db$Markers)
         existing@data_type = "AvgExp & Markers"
      }else if(data_type=="AvgExp"){
         existing@query$AvgExp = data.frame(temp_ann@query$AvgExp)
         existing@db$AvgExp = data.frame(temp_ann@db$AvgExp)
         existing@data_type = "Markers & AvgExp"
      }
      return(existing)
   }


}
