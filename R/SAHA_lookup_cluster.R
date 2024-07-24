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
SAHA_lookup_cluster <- function(cluster,ann) {
   k=ann@ann1$query[ann@ann1$query$cluster==cluster,"gene"]
   return(summary(factor(ann@ann1$db[ann@ann1$db$SYMBOL %in% k, "cluster"], levels(factor(ann@ann1$db$cluster)))))
}
