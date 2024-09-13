#' Dummy query marker data
#'
#' Data contains
#' @format A dataframe with ...:
#'  \describe{
#'  \item{X}{Running index of gene*cluster ID}
#'  \item{p_val}{Uncorrected p value for gene in cluster-of-interest vs. ALL}
#'  \item{avg_log2FC}{Average log2FC of expression for gene in cluster-of-interest vs. ALL}
#'  \item{pct.1}{Percent of cells in cluster-of-interest to express gene}
#'  \item{pct.2}{Percent of cells not in cluster-of-interest to express gene}
#'  \item{p_val_adj}{Adjusted p value for gene in cluster-of-interest vs. ALL}
#'  \item{cluster}{Cluster-of-interest}
#'  \item{gene}{ENSMBL ID for Gene}
#'  \item{SYMBOL}{Gene Symbol}}
#' @source {https://portal.brain-map.org/atlases-and-data/rnaseq}
#'
#' @examples
#' data(x) #will be lazy loading
"x"

#' Dummy db marker data
#'
#' Data contains
#' @format A dataframe with ...:
#'  \describe{
#'  \item{X}{Running index of gene*cluster ID}
#'  \item{p_val}{Uncorrected p value for gene in cluster-of-interest vs. ALL}
#'  \item{avg_log2FC}{Average log2FC of expression for gene in cluster-of-interest vs. ALL}
#'  \item{pct.1}{Percent of cells in cluster-of-interest to express gene}
#'  \item{pct.2}{Percent of cells not in cluster-of-interest to express gene}
#'  \item{p_val_adj}{Adjusted p value for gene in cluster-of-interest vs. ALL}
#'  \item{cluster}{Cluster-of-interest}
#'  \item{gene}{ENSMBL ID for Gene}
#'  \item{SYMBOL}{Gene Symbol}}
#' @source {https://portal.brain-map.org/atlases-and-data/rnaseq}
#'
#' @examples
#' data(y) #will be lazy loading
"y"
