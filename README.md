
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Chevreul

This package includes a set of Shiny apps for exploring single cell RNA
datasets processed with
<a href="https://github.com/satijalab/object" target="_blank" rel="noopener noreferrer">Seurat</a>

A demo using a human gene transcript dataset from Shayler et al. (link)
is available
<a href="http://cobrinik-1.saban-chla.usc.edu:8080/app/objectApp" target="_blank" rel="noopener noreferrer">here</a>

There are also convenient functions for:

- Clustering and Dimensional Reduction of Raw Sequencing Data.
- <a href="https://satijalab.org/seurat/archive/v3.0/integration.html" target="_blank" rel="noopener noreferrer">Integration
  and Label Transfer</a>
- Louvain Clustering at a Range of Resolutions
- Cell cycle state regression and labeling
- RNA velocity calculation with
  <a href="https://velocyto.org/" target="_blank" rel="noopener noreferrer">Velocyto.R</a>
  and
  <a href="https://scvelo.readthedocs.io/" target="_blank" rel="noopener noreferrer">scvelo</a>

> \[!WARNING\] Chevreul was designed for full-length smart-seq based
> single cell data. Default settings may not be appropriate for droplet
> (10x) data, though most can be adjusted. Keep in mind [best
> practices](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
> regarding normalization, dimensional reduction, etc. when using.

## Installation

You can install the released version of chevreul from
<a href="https://github.com/whtns/chevreul" target="_blank" rel="noopener noreferrer">github</a>
with:

### Install locally and run in three steps:

You can install Chevreul locally using the following steps:

``` r
install.packages("devtools")
devtools::install_github("whtns/chevreul")
chevreul::create_project_db()
```

You can also customize the location of the app using these steps:

``` r
devtools::install_github("whtns/chevreul")
chevreul::create_project_db(destdir = "/your/path/to/app")
```

## Getting Started

First, load Chevreul and all other packages required

``` r
library(chevreul)
library(Seurat)
library(tidyverse)
library(ggraph)
```

## TLDR

Chevreul provides a single command to:

- construct a Seurat object

- filter genes by minimum expression and ubiquity

- normalize and scale expression by any of several methods packaged in
  Seurat

## Run clustering on a single object

By default clustering will be run at ten different resolutions between
0.2 and 2.0. Any resolution can be specified by providing the resolution
argument as a numeric vector.

``` r
clustered_object <- clustering_workflow(human_gene_transcript_object,
    experiment_name = "object_hu_trans",
    organism = "human"
)
```

## Get a first look at a processed dataset using an interactive shiny app

``` r
minimalSeuratApp(human_gene_transcript_object)
```

## Set up a object

We start with a gene by cell matrix of count/UMI values and a table of
cell metadata

``` r
human_count[1:5, 1:5]

head(human_meta)
```

We can then create a object in the usual manner using
`CreatSeuratObject` function

``` r
myobject <- CreateSeuratObject(human_count, assay = "gene", meta.data = human_meta)
```

## Preprocess the object

Chevreul includes a handy function to preprocess the data that handles
normalization and scaling required for downstream analysis. If needed,
parameters can be specified by the user.

``` r
myobject <- object_preprocess(myobject)
```

This single function includes object sub-functions that normalizes,
identifies highly variable features and scales the data

## Perform dimension reduction

Chevreul also implements a standardized dimension reduction step to
select variable features at a user-specified threshold and perform PCA,
tSNE, and UMAP. The default assay the dimension reduction is being run
on is “gene”.

``` r
myobject <- object_reduce_dimensions(myobject, assay = "RNA")
```

## Community detection by clustering

Clustering analysis is performed via Louvain(default) or alternative
algorithms available in Seurat. Clustering is performed at a range of
resolutions with default value ranging from 0.2 to 2 and pca reduction

``` r
object <- object_cluster(object = Dim_Red_object, resolution = seq(0.2, 2, by = 0.2))
```

This function produces clustering analysis via two steps performed using
two different sub-functions

## Split included dataset based on collection technology

Chevreul includes a function, `SplitObject`, which is capable of
splitting the dataset into subsets based on a single attribute indicated
by the split.by argument

``` r
split_human <- SplitObject(human_gene_transcript_object, split.by = "dataset")
```

In this example the `split_human` object consists of a list of subsetted
objects which are split based on batch

## Run object batch integration on ‘child’ projects

When joint analysis of 2 or more datasets is to be performed
`integration_workflow` function can be used, which takes in a list of
objects as input and returns an integrated object

``` r
integrated_object <- integration_workflow(split_human)
```

## View analysis details

``` r
Misc(integrated_object, "experiment") %>%
    tibble::enframe() %>%
    knitr::kable()
```
