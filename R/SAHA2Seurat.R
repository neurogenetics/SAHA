#' SAHA2Seurat
#'
#' Adds annotations from a specified annotation database to a Seurat object. The function renames identities in the Seurat object based on cluster annotations and can add the annotations to the object's metadata if specified.
#'
#' @param obj A Seurat object in which identities will be renamed according to the annotation database.
#' 
#' @param anno_db A data frame containing annotation information with at least two columns:
#'        'cluster' for cluster IDs 
#'        'best_match' for annotation labels.
#'        
#' @param anno A string specifying the column in `anno_db` to use for renaming identities; defaults to "best_match".
#' 
#' @param anno2meta A string specifying the name of the metadata field in which to store the annotations; defaults to "SAHA_anno". If NULL, annotations are only set as Idents(obj).
#'
#' @return A Seurat object with updated identities and optional metadata field with annotations.
#' @importFrom Seurat RenameIdents Idents
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object and `annotations_db` is a data frame
#' # containing cluster annotations:
#' seurat_obj <- SAHA2Seurat(seurat_obj, annotations_db)
SAHA2Seurat <- function(obj, anno_db, anno = "best_match", anno2meta = "SAHA_anno") {
  
  # Figure out which column names are present
  if (all(c("cluster", "best_match") %in% colnames(anno_db))) {
    cluster_col <- "cluster"
    label_col   <- "best_match"
  } else if (all(c("old_names", "new_names") %in% colnames(anno_db))) {
    cluster_col <- "old_names"
    label_col   <- "new_names"
  } else {
    stop("anno_db must contain either columns 'cluster' and 'best_match' or 'old_names' and 'new_names'.")
  }
  
  # Build named vector for RenameIdents
  annotations <- anno_db[, c(cluster_col, label_col)]
  new.names   <- annotations[[label_col]]
  names(new.names) <- annotations[[cluster_col]]
  
  obj <- RenameIdents(obj, new.names)
  
  if (is.null(anno2meta)) {
    cat("SAHA annotations written to Idents(obj).")
  } else {
    obj[[anno2meta]] <- Idents(obj)
    cat(paste0("SAHA annotations written to Idents(obj) and to obj@meta.data$", anno2meta, "."))
  }
  
  return(obj)
}
