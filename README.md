<!-- README.md is generated from README.Rmd. Please edit that file -->

# SAHA

<!-- badges: start -->
<!-- badges: end -->

<u>S</u>emi-<u>a</u>utomated <u>H</u>and <u>A</u>nnotation for Single
Cell and Spatial Datasets

<img src="files/SAHA_logo_wbg.png" width="188" />

## Installation

You can install the development version of SAHA like so:

``` r
library(devtools)
devtools::install_github("neurogenetics/SAHA")
```

## Quickstart

This is a basic example which shows you how to solve a common problem:

``` r
library(SAHA)
#> Thank you for using SAHA. If you enjoy this package, please consider citing Acri et al., (XXXX) bioRxiv.
data("ABC_meta")
data("CB_Markers")#Allen Brain Atlas Markers from Cerebellum
data("CB_AvgExp")#Allen Brain Atlas AvgExp from Cerebellum
#My_Markers #your marker data
#My_AvgExp #your avgexp data
## Quickstart to Marker-based dotplot
#SAHA(query = My_markers, db = CB_Markers,meta = meta, data_type = "Markers")
#Quickstart to Marker-free heatmap
#SAHA(query = My_AvgExp, db = CB_AvgExp,meta = meta, data_type = "AvgExp")
```
