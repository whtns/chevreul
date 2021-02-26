
<!-- README.md is generated from README.Rmd. Please edit that file -->

# seuratTools

<!-- badges: start -->
<!-- badges: end -->
<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/whtns/seuratTools.svg?branch=master)](https://travis-ci.org/whtns/seuratTools)
<!-- badges: end -->

This project is not associated with the Seurat team; It will likely
[change
names](https://www.njtierney.com/post/2017/10/27/change-pkg-name/) in
the future

# Seurat Tools

This package includes a set of Shiny apps for exploring single cell RNA
datasets processed with [Seurat](https://github.com/satijalab/seurat)
objects.

A demo using a pancreas dataset from the Seurat team is availble
[here](http://cobriniklab.saban-chla.usc.edu:3838/seuratTools_demo/)

There are also convenient functions for: 1. Clustering and Dimensional
Reduction of Raw Sequencing Data 2. [Integration and Label
Transfer](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html)
3. Louvain Clustering at a Range of Resolutions 4. Cell cycle state
regression and labeling 5. RNA velocity calculation with
[Velocyto.R](https://velocyto.org/) and
[scvelo](https://scvelo.readthedocs.io/)

## Installation

You can install the released version of seuratTools from
[github](https://github.com/whtns/seuratTools) with:

### Install locally and run in three steps:

    install.packages("devtools")
    devtools::install_github("whtns/seuratTools")
    seuratTools::create_project_db()

### Install locally (custom location!) and run in three steps:

    devtools::install_github("whtns/seuratTools")
    seuratTools::create_project_db(destdir='/your/path/to/app')

## Site

You can view documentation on the [seuratTools
website](https://whtns.github.io/seuratTools)

## How To

### subset by csv

![subset by csv](README_docs/subset_by_csv.gif)

### add custom metadata

![add custom metadata](README_docs/add_arbitrary_metadata.gif)

## Getting Started

    library(seuratTools)
    library(Seurat)
    library(tidyverse)
    library(ggraph)

### view included dataset

    panc8

### run clustering on a single seurat object

By default clustering will be run at ten different resolutions between
0.2 and 2.0. Any resolution can be specified by providing the resolution
argument as a numeric vector.

    clustered_seu <- clustering_workflow(panc8, experiment_name = "seurat_pancreas", organism = "human")

    minimalSeuratApp(clustered_seu)

### split included dataset based on collection technology

    split_panc8 <- SplitObject(panc8, split.by = "dataset")

### run seurat batch integration on ‘child’ projects

    integrated_seu <- integration_workflow(split_panc8)

### launch app to inspect


    minimalSeuratApp(integrated_seu)

### view analysis details

    Misc(integrated_seu, "experiment") %>% 
      tibble::enframe() %>% 
      knitr::kable()
