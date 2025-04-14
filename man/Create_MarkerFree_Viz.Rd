\name{Create_MarkerFree_Viz}
\alias{Create_MarkerFree_Viz}
\title{Marker-free annotation: Creating marker-free annotation heatmap}
\description{Generates a heatmap of correlation of gene expression profiles between query dataset and database.}
\usage{
Create_MarkerFree_Viz(ann, meta = NULL, ABC = NULL)
}
\arguments{
  \item{ann}{A SAHA object which has been correlated using CorrelateDS() for marker-free annotation.}
  \item{meta}{Optional dataframe containing metadata about database cell types. If provided, it enables additional annotations on the heatmap. The dataframe requires several specific columns:
    \itemize{
      \item "subclass_per": names of the cell types from the database (words must be separated by periods).
      \item "neurotransmitter": highest-level (broadest) classification for each cell type.
      \item "class": more specific level of classification (between subclass and neurotransmitter).
      \item "neurotransmitter_color": hex code for the color of the "neurotransmitter" classification.
      \item "class_color": hex code for the color of the "class" classification.
    }
  }
}
\value{
Returns the modified \code{ann} object with the heatmap stored in \code{ann@results$marker_free$heatmap}.
}
\details{
This function generates a heatmap representing correlation coefficients between query and database expression profiles. If metadata is provided, the heatmap includes annotations for neurotransmitter and class. If ABC is specified, it filters the data based on the specified brain region in the Allen Brain Atlas.

The heatmap is created using the ComplexHeatmap package. It clusters both rows and columns, and uses a color scale from blue (negative correlation) through white (no correlation) to orange (positive correlation).
}
\examples{
# Basic usage without metadata
Create_MarkerFree_Viz(ann)

# With metadata
Create_MarkerFree_Viz(ann, meta = cell_type_metadata)

}
