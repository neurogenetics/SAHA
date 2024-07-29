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
NormalizeDS <- function(ann,assay_db = "RNA",assay_query){
   #1. Pivot to a new table - query
   piv_query <- ann@ann3$query %>%
      tibble::rownames_to_column(var = "Row") %>%
      pivot_longer(cols = -Row, names_to = "cluster", values_to = "Expression")

   # average and make percentiles

   grouped_data_query <- split(piv_query, piv_query$Row)

   zscore_query <- lapply(grouped_data_query, function(gene_data) {
      gene_data$z_score <- scale(gene_data$Expression)
      return(gene_data)
   })

   avgexp_query.df <- do.call(rbind, zscore_query)
   avgexp_query.df$z_score <- as.numeric(avgexp_query.df$z_score)
   avgexp_query.df$percentile <- pnorm(avgexp_query.df$z_score)

   #2. Pivot to a new table - db
   piv_db <- ann@ann3$db %>%
      tibble::rownames_to_column(var = "Row") %>%
      pivot_longer(cols = -Row, names_to = "celltype", values_to = "Expression")

   # average and make percentiles

   grouped_data_db <- split(piv_db, piv_db$Row)

   zscore_db <- lapply(grouped_data_db, function(gene_data) {
      gene_data$z_score <- scale(gene_data$Expression)
      return(gene_data)
   })

   avgexp_db.df <- do.call(rbind, zscore_db)
   avgexp_db.df$z_score <- as.numeric(avgexp_db.df$z_score)
   avgexp_db.df$percentile <- pnorm(avgexp_db.df$z_score)

   #3. making comparison b/w query and ABC expression profiles for major cell types

   avgexp_db.df$celltype <- gsub(paste0("^",assay_db,"\\.g"), "", avgexp_db.df$celltype) ###########!!!!!!!!!!!need to print "Assay" into here!!
   avgexp_query.df$cluster <- gsub(paste0("^",assay_query,"\\.g"), "", avgexp_query.df$cluster) ###########!!!!!!!!!!!need to print "Assay" into here!!


   avgexp_db.df <- avgexp_db.df[c("Row", "celltype", "percentile")]
   avgexp_query.df <- avgexp_query.df[c("Row", "cluster", "percentile")]



   # rearrange dataframes

   avgexp_db.df <- avgexp_db.df %>%
      pivot_longer(
         cols = starts_with("percentile"),
         names_to = "platform",
         values_to = "expression"
      ) %>%
      mutate(platform = "db") %>%
      unite("platform_celltype", platform, celltype, sep = ".") %>%
      pivot_wider(
         names_from = platform_celltype,
         values_from = expression
      )


   avgexp_query.df <- avgexp_query.df %>%
      pivot_longer(
         cols = starts_with("percentile"),
         names_to = "platform",
         values_to = "expression"
      ) %>%
      mutate(platform = "query") %>%
      unite("platform_cluster", platform, cluster, sep = ".") %>%
      pivot_wider(
         names_from = platform_cluster,
         values_from = expression
      )


   merged.df <- merge(avgexp_db.df, avgexp_query.df, by = "Row", all.x = T, all.y = T)

   merged.df <- column_to_rownames(merged.df, var = "Row")
   ann@results$marker_free$norm_merge=merged.df
   #this need to be made into a list
   return(ann)

}
