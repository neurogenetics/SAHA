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
Initialize_MarkerFree <- function(ann){
   #####CLEANING NEEDS TO BE PERFORMED PRIOR TO SAHA
   if ("X"%in%colnames(ann@query$AvgExp)) {
      rownames(ann@query$AvgExp)=ann@query$AvgExp$X
      ann@query$AvgExp = ann@query$AvgExp %>%
         select(-X)
   }
   if (all(!duplicated(ann@db$AvgExp$SYMBOL))) {
      rownames(ann@db$AvgExp)=ann@db$AvgExp$SYMBOL
      ann@db$AvgExp = ann@db$AvgExp %>%
         select(-SYMBOL, -gene, -X)
   }else{
      ann@db$AvgExp=ann@db$AvgExp[!duplicated(ann@db$AvgExp$SYMBOL),]
      ann@db$AvgExp=ann@db$AvgExp[!is.na(ann@db$AvgExp$SYMBOL),]
      rownames(ann@db$AvgExp)=ann@db$AvgExp$SYMBOL
      ann@db$AvgExp = ann@db$AvgExp %>%
         select(-SYMBOL, -gene, -X)
   }
   return(ann)
}
