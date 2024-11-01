
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SAHA

<!-- badges: start -->

<!-- badges: end -->

<u>S</u>emi-<u>a</u>utomated <u>H</u>and <u>A</u>nnotation for Single
Cell and Spatial Datasets

<img src="files/SAHA_logo_wbg.png" width="188"/>

<b>Beta Prior to Release:</b> If you chose to use, please check back regularly for updates and improvements. <b>Version 1.0</b> will indicate stable release and contain links to publication or preprint.

## Useful links
[Installation](https://github.com/neurogenetics/SAHA/edit/main/README.md#installation)

[Getting started with SAHA](https://github.com/neurogenetics/SAHA/wiki/SAHA-Main-Workflow)

## Why Use SAHA to annotate your single cell or spatial data?
SAHA is a user-friendly package with simple meta-data level input resulting in easy-to-understand exploration of cell annotation in single cell and spatial RNA sequencing datasets. Using either marker gene or average experession dataframes, the user has the option to compare the similarity of their <i>unannotated</i> clusters to <b>any</b> annotations in the literature.

<b>IMPORTANTLY:</b> This approach differs from label-transfer and integration-based approaches in that it does not require a large Seurat, SCE, Anndata, or any common data type in which the computational burden of the pipeline scales with the size of an object. Most (if not all) data fed into SAHA will be able to run on a laptop without the need for high-RAM devices or supercomputers.

Whether you can run every single cell package ever released or just finished your first single cell vignette, SAHA offers the options to <b>Create</b>, <b>Investigate</b>, or <b>Refine</b> cell annotations based on:
* Self-similarity: likeness of clusters within your query data
* Marker-based: presence of markers in a known database
* Marker-free: expression of <b>non-marker genes</b> in your data

## Installation

### Dependencies
* [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
* [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
* dplyr
* methods
* [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)
* [eulerr](https://cran.r-project.org/web/packages/eulerr/index.html)
* reshape2
* corrplot
* ggpubr
* ape
* vegan
* dendextend
   
### Install development version from github using devtools

You can install the development version of SAHA like so:

``` r
library(devtools)
devtools::install_github("neurogenetics/SAHA")#main package
devtools::install_github("neurogenetics/SAHAdata")#data package containing Allen Brain Cell Atlas db inputs
```
Upon stable release, instructions for CRAN installation will be posted here.

