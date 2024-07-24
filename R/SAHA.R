# SAHA
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
#################################################
#Rework so this is the quickstart option########
#################################################


SAHA <- function(query, database = "Panglao", nameoforgan = NULL, mode = "db", speciesquery = NULL, speciesdatabase = NULL) {

   #1. Read in your queried data as a df or csv
   if (inherits(query, "data.frame")) {
      # If query is already a data frame, use it directly
      query_data <- query
   } else {
      # If query is a path to a CSV file, read it into a data frame
      query_data <- read.csv(query)
   }

   #2. Gene lookup
   hmf <- read.csv(paste0(path,"20230911_HumanMouseFly_gene_nomenclature.csv"))
   na.omit(hmf[hmf$human_symbol.x == "MAPT", ])
   na.omit(hmf[hmf$mouse_symbol == "Mapt", ])
   na.omit(hmf[hmf$fruit.fly_symbol == "tau", ])

   #3. Decide which mode to proceed in.
   if (mode == "1v1") {
   } else {
      #a)Read in database of choice
      #i) PanglaoDB
      if (database == "Panglao") {
         data <- read.delim(paste0(path,"PanglaoDB_markers_27_Mar_2020.tsv"), sep = "\t", header = TRUE)
      }
      #ii) Custom
      if (database != "Panglao" & database != "cellmarker2.0")  {
         if (inherits(database, "data.frame")) {
            # If query is already a data frame, use it directly
            data <- database
         } else {
            # If query is a path to a CSV file, read it into a data frame
            data=read.csv(database)
         }
      }
      #iii) Subset Custom or Panglao by a given organ
      if (!is.null(nameoforgan)) {
         marker_data <- subset(data, organ == nameoforgan)
      } else {
         marker_data <- data
      }

      #iv). Cellmarker2.0
      if (database == "cellmarker2.0") {
         # Cell Marker2.0 all cell type download 2023101
         # http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html
         data=data.frame(read_xlsx(paste0(path,"Cell_marker_All.xlsx"), col_names = TRUE))
         if (!is.null(nameoforgan)) {
            marker_data <- subset(data, organ == nameoforgan)
         } else {
            marker_data <- data
         }
         colnames(marker_data)[colnames(marker_data) == "tissue_type"] <- "species"
         # Change column name for downstream workflow
         colnames(marker_data)[colnames(marker_data) == "cell_name"] <- "cell.type"
         colnames(marker_data)[colnames(marker_data) == "Symbol"] <- "official.gene.symbol"
      }


      #4. Function to look up any cluster
      SAHA_lookup_cluster <- function(cluster) {
         # Need to convert i to human genes for pangloa
         i=query_data[query_data$cluster==cluster,"gene"]
         # Convert genes to human orthologs
         j=unique(na.omit(hmf[hmf$mouse_symbol%in%i,"human_symbol.x"]))
         return(summary(factor(marker_data[marker_data$official.gene.symbol %in% j, "cell.type"], levels(factor(marker_data$cell.type)))))
      }

      #5. Making our Master Data Frame
      master_df <- data.frame(summary(factor(marker_data$cell.type)))
      master_df$cluster="REF"
      colnames(master_df)[1]<-"total_marker"
      master_df$celltype <- rownames(master_df)

      #6. Loop through each cluster, find their markers, bind to masterdf
      for (i in c(1:length(unique(query_data$cluster)))) {
         x=master_df[master_df$cluster=="REF",]
         x$cluster = unique(query_data$cluster)[i]
         x$total_marker=SAHA_lookup_cluster(unique(query_data$cluster)[i])
         x$celltype = rownames(x)
         master_df=rbind(master_df,x)
      }

      #7. Find proportion
      master_df$prop="NA"
      for (j in c(1:length(unique(query_data$cluster)))) {
         master_df[master_df$cluster==unique(query_data$cluster)[j],"prop"]=master_df[master_df$cluster==unique(query_data$cluster)[j],"total_marker"]/master_df[master_df$cluster=="REF","total_marker"]
      }

      #8. Find p value
      master_df$pvalue = 1
      for (j in c(1:length(unique(query_data$cluster)))) {
         # Size of possible markers in a given cluster
         A = length(query_data[query_data$cluster == unique(query_data$cluster)[j], "gene"])
         # Size of signature testing (number of genes in pangloa cell type)
         for (k in c(1:length(unique(master_df$celltype)))) {
            B = master_df[master_df$cluster=="REF",][k,"total_marker"]
            # Overlap
            t = master_df[master_df$cluster==unique(query_data$cluster)[j],][k,"total_marker"]
            # Length of all possible cluster markers (not just the one testing)
            n = length(unique(rownames(query_data)))
            master_df[master_df$cluster==unique(query_data$cluster)[j],][k,"pvalue"]=sum(stats::dhyper(t:B, A, n - A, B))
         }
      }

      #9. Conditional Facet
      if (database == "Panglao") {
         # Label the cell types
         master_df <- master_df %>%
            mutate(cell_group = case_when(
               celltype %in% c("Adrenergic neurons", "Cajal-Retzius cells", "Cholinergic neurons",
                               "Dopaminergic neurons", "GABAergic neurons", "Glutaminergic neurons",
                               "Glycinergic neurons", "Immature neurons", "Interneurons",
                               "Motor neurons","Neurons", "Noradrenergic neurons","Purkinje neurons",
                               "Pyramidal cells", "Retinal ganglion cells", "Serotonergic neurons",
                               "Trigeminal neurons") ~ "Neurons",
               celltype %in% c("Astrocytes", "Bergmann glia", "Glial cell", "Microglia", "Oligodendrocytes",
                               "Radial glia cells", "Satellite glial cells", "Schwann cells", "Tanycytes") ~ "Glia",
               celltype %in% c("Choroid plexus cells", "Ependymal cells", "Meningeal cells") ~ "EC",
               celltype %in% c("Anterior pituitary gland cells", "Neuroendocrine cells", "Pinealocytes") ~ "EN",
               celltype %in% c("Neural stem/precursor cells", "Neuroblasts", "Oligodendrocyte progenitor cells",
                               "Retinal ganglion cells") ~ "Precursor",
               TRUE ~ "Other"  # Default category for other cell types not specified
            )) %>%
            ungroup()

         # Reorder levels of celltype
         master_df$celltype <- fct_relevel(master_df$celltype, "Neurons", "Neuroendocrine cells")
      }

      #10. Default everything to F, not significant
      master_df$sig = "F"
      master_df[master_df$pvalue <= 0.05,"sig"]="T" # Label only significant groups

      # Convert the 'cluster' variable to a factor with custom levels in ascending order
      master_df$cluster <- factor(master_df$cluster, levels = rev(unique(master_df$cluster)))

      # Create dataframe that only finds minimum pvalue for each cluster
      master_df <- master_df %>%
         group_by(cluster) %>%
         mutate(top_sig = ifelse(pvalue == min(pvalue) & sig == "T", "T", "F"),
                top_sig_prop = ifelse(top_sig == "F", 0, as.numeric(prop))) %>%
         ungroup()
      master_df <- data.frame(master_df)

      invisible(list(master_df = master_df, query_data = query_data, marker_data = marker_data ))
   }
}

