<!-- badges: start -->
[![Travis build status](https://travis-ci.org/whtns/seuratTools.svg?branch=master)](https://travis-ci.org/whtns/seuratTools)
<!-- badges: end -->

# seuratTools

Howdy y'all! this is a set of convenience functions for interacting with [Seurat](https://github.com/satijalab/seurat) objects. There are functions for:
1. Creating Seurat objects from Stringtie files
2. [Integration and Labe lTransfer](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html)
3. Louvain Clustering at a Range of Resolutions 
4. Cell cycle state regression and labeling 
5. RNA velocity calculation with [Velocyto.R](https://velocyto.org/)


## Installation

You can install the released version of seuratTools from [our github](https://github.com/whtns/seuratTools) with:

``` r
devtools::install_package("whtns/seuratTools")
```

## Site

You can view documentation on the [seuratTools website](https://whtns.github.io/seuratTools)

## How To 

### subset by csv

![subset by csv](README_docs/subset_by_csv.gif)

### add custom metadata

![add custom metadata](README_docs/add_arbitrary_metadata.gif)

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```


