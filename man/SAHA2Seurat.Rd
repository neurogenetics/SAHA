\name{SAHA2Seurat}
\alias{SAHA2Seurat}
\title{Add SAHA Annotations from Database to Seurat Object}
\description{
Adds annotations from a specified annotation database to a Seurat object by renaming identities based on cluster annotations. It can also add these annotations to the object's metadata if specified.
}
\usage{
SAHA2Seurat(obj, anno_db, anno = "best_match", anno2meta = "SAHA_anno")
}
\arguments{
  \item{obj}{A \code{Seurat} object in which identities will be renamed according to the annotation database.}

  \item{anno_db}{A \code{data.frame} containing annotation information with at least two columns:
    \code{cluster} for cluster IDs and \code{best_match} for annotation labels.}

  \item{anno}{A \code{character} string specifying the column in \code{anno_db} to use for renaming identities; defaults to \code{"best_match"}.}

  \item{anno2meta}{A \code{character} string specifying the name of the metadata field in which to store the annotations;
    defaults to \code{"SAHA_anno"}. If set to \code{NULL}, annotations are only set as \code{Idents(obj)}.}
}
\value{
Returns the modified \code{Seurat} object with updated identities and, if specified, a metadata field containing annotations.
}
\details{
This function allows for the addition of custom annotations to a \code{Seurat} object by renaming its identities based on a user-provided annotation database.
If \code{anno2meta} is specified, the annotations are also added to the metadata of the Seurat object, making it accessible in both \code{Idents(obj)} and \code{obj@meta.data}.
}
\examples{
# Example usage of the SAHA2Seurat function
# Assuming \code{seurat_obj} is a Seurat object and \code{annotations_db} is a data frame
# containing cluster annotations with columns \code{cluster} and \code{best_match}:
seurat_obj <- SAHA2Seurat(seurat_obj, annotations_db)
}
\seealso{
\code{\link[Seurat]{Idents}}, \code{\link[Seurat]{RenameIdents}}
}

