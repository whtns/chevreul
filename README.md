
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Chevreul

This package includes a set of Shiny apps for exploring single cell RNA
datasets processed with
<a href="https://github.com/satijalab/seurat" target="_blank" rel="noopener noreferrer">Seurat</a>

A demo using a human gene transcript dataset from the Seurat team is
available
<a href="http://cobrinik-1.saban-chla.usc.edu:8080/app/0seuratApp" target="_blank" rel="noopener noreferrer">here</a>

There are also convenient functions for:

-   Clustering and Dimensional Reduction of Raw Sequencing Data.
-   <a href="https://satijalab.org/seurat/archive/v3.0/integration.html" target="_blank" rel="noopener noreferrer">Integration
    and Label Transfer</a>
-   Louvain Clustering at a Range of Resolutions
-   Cell cycle state regression and labeling
-   RNA velocity calculation with
    <a href="https://velocyto.org/" target="_blank" rel="noopener noreferrer">Velocyto.R</a>
    and
    <a href="https://scvelo.readthedocs.io/" target="_blank" rel="noopener noreferrer">scvelo</a>

## Installation

You can install the released version of chevreul from
<a href="https://github.com/whtns/chevreul" target="_blank" rel="noopener noreferrer">github</a>
with:

### Install locally and run in three steps:

``` r
install.packages("devtools")
devtools::install_github("whtns/chevreul")
chevreul::create_project_db()
```

### Install locally (custom location!) and run in three steps:

``` r
devtools::install_github("whtns/chevreul")
chevreul::create_project_db(destdir='/your/path/to/app')
```

## Getting Started

``` r
library(chevreul)
library(Seurat)
library(tidyverse)
library(ggraph)
```

## TLDR

chevreul provides a single command to:

-   construct a Seurat object

-   filter genes by minimum expression and ubiquity

-   normalize and scale expression by any of several methods packaged in
    Seurat

## Run clustering on a single seurat object

By default clustering will be run at ten different resolutions between
0.2 and 2.0. Any resolution can be specified by providing the resolution
argument as a numeric vector.

``` r
clustered_seu <- clustering_workflow(human_gene_transcript_seu, experiment_name = "seurat_pancreas", organism = "human")
```

## Get a first look at a processed dataset using an interactive shiny app

``` r
minimalSeuratApp(clustered_seu)
```

## Set up a seurat object

We start with a gene by cell matrix of count/UMI values and a table of
cell metadata

``` r
human_count[1:5, 1:5]

head(human_meta)
```

we can create a seurat object in the usual manner

``` r
myseu <- CreateSeuratObject(human_count, assay = "gene", meta.data = human_meta)
```

## Preprocess the seurat object

Chevreul includes a handy function to preprocess the data that handles
normalization and scaling required for downstream analysis.
Preprocessing is performed using existing Seurat functions. If needed,
parameters can be specified by the user.

``` r
myseu <- seurat_preprocess(myseu)
```

This single function includes sub-functions that normalizes, identifies
highly variable features and scales the data:

-   The normalizing sub-function performs log normalization using a
    default scale factor of 10,000.

``` r
preprocess_seu <- NormalizeData(myseu, verbose = FALSE)
```

-   After normalization, subset of features that exhibit high
    cell-to-cell variation in the dataset are identified. By default,
    2,000 features per dataset are returned by this function.

``` r
preprocess_seu <- FindVariableFeatures(preprocess_seu, selection.method = "vst", 
    verbose = FALSE)
```

-   Finally, the data is scaled by applying linear transformation. This
    step shifts the gene expression, so that the mean expression across
    cells is 0 and scales the gene expression, so that the variance
    across cells is 1.

``` r
pre_process_seu <- ScaleData(preprocess_seu)
```

## Perform dimension reduction

Chevreul also implements a standardized dimension reduction step to
select variable features at a user-specified threshold and perform PCA,
tSNE, and UMAP. The default assay the dimension reduction is being run
on is “gene”.

``` r
myseu <- seurat_reduce_dimensions(myseu, assay = "RNA") 
```

This function includes existing seurat functions which performs
dimension reduction techniques.

-   Perform PCA: Runs a PCA dimensionality reduction.

``` r
Dim_Red_seu <- RunPCA(myseu,features = VariableFeatures(myseu), 
      do.print = FALSE)
```

-   Perform tSNE: Runs t-SNE dimensionality reduction on selected
    features

``` r
Dim_Red_seu <- RunTSNE(Dim_Red_seu, dims = 1:30)
```

-   Perform UMAP: Runs the Uniform Manifold Approximation and Projection
    (UMAP) dimensional reduction technique.

``` r
Dim_Red_seu <- RunUMAP(Dim_Red_seu, dims = 1:30)
DimPlot(Dim_Red_seu, reduction = "umap")
```

## Community detection by clustering

Clustering analysis is performed via Louvain(default) or alternative
algorithms available in Seurat. Clustering is performed at a range of
resolutions with default value ranging from 0.2 to 2 and pca reduction

``` r
 seu <- seurat_cluster(seu = Dim_Red_seu, resolution = seq(0.2, 2, by = 0.2) )
```

This function produces clustering analysis via two steps performed using
two different sub-functions

-   `FindNeighbours`: This function computes the nearest neighbors for a
    given dataset using k-nearest neighbor algorithm.

-   `FindClusters`: The output from FindNeighbours is then used to
    compute

## Split included dataset based on collection technology

Chevreul includes a function, `SplitObject`, which is capable of
splitting the dataset into subsets based on a single attribute indicated
by the split.by argument

``` r
split_panc8 <- SplitObject(panc8, split.by = "dataset")
```

Here, the split_panc8 object consists of a list of subsetted objects
which are split based on batch

## Run seurat batch integration on ‘child’ projects

Multiple seurat objects can be integrated using a Chevreul function
which takes in a list of seurat objects for all batches as an argument
and returns an integrated seurat object containg merged batches

``` r
integrated_seu <- integration_workflow(split_panc8)
```

## launch app to inspect

``` r
minimalSeuratApp(integrated_seu)
```

## view analysis details

``` r
Misc(integrated_seu, "experiment") %>% 
  tibble::enframe() %>% 
  knitr::kable()
```
