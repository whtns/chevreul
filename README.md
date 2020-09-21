<!-- badges: start -->
[![Travis build status](https://travis-ci.org/whtns/seuratTools.svg?branch=master)](https://travis-ci.org/whtns/seuratTools)
<!-- badges: end -->

# Seurat Tools

This is a set of convenience functions for interacting with [Seurat](https://github.com/satijalab/seurat) objects. 
There are functions for:
1. Clustering and Dimensional Reduction of Raw Sequencing Data
2. [Integration and Label Transfer](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html)
3. Louvain Clustering at a Range of Resolutions 
4. Cell cycle state regression and labeling 
5. RNA velocity calculation with [Velocyto.R](https://velocyto.org/)


## Installation

You can install the released version of seuratTools from [our github](https://github.com/whtns/seuratTools) with:

Install locally and run in three steps:

```
devtools::install_github("whtns/seuratTools")
seuratTools::create_project_db()
```

# Install locally (custom location!) and run in three steps:
```
devtools::install_github("whtns/seuratTools")
seuratTools::create_project_db(destdir='/your/path/to/app')
seuratTools::seuratApp(proj_dir, feature_types = c("gene", "transcript"))
```


## Site

You can view documentation on the [seuratTools website](https://whtns.github.io/seuratTools)

## How To 

### subset by csv

![subset by csv](README_docs/subset_by_csv.gif)

### add custom metadata

![add custom metadata](README_docs/add_arbitrary_metadata.gif)


```
library(seuratTools, lib.loc = "~/rpkgs/devel_install/")
library(tidyverse)
library(Seurat)
library(ggraph)
library(iheatmapr)
library(formattable)
```

# view included dataset 

```
seurat_pancreas_reduced

glimpse(seurat_pancreas_reduced)
```

# run clustering on a single seurat object

By default clustering will be run at ten different resolutions between 0.2 and 2.0. Any resolution can be specified by providing the resolution argument as a numeric vector.

```

clustered_seu <- clustering_workflow(seurat_pancreas_reduced, experiment_name = "seurat_pancreas", organism = "human")
```

```
minimalSeuratApp(clustered_seu)
```

## split included dataset based on collection technology 

```
batches <- seurat_pancreas_reduced %>%
 purrr::map(Seurat::SplitObject, split.by = "dataset") %>%
 purrr::transpose()

names(batches)

glimpse(batches)

```

# run seurat batch integration on 'child' projects

```
integrated_seu <- integration_workflow(batches)
```

# launch app to inspect

```

minimalSeuratApp(integrated_seu)

```

# view analysis details

```
integrated_seu$gene@misc$experiment %>% 
  tibble::enframe() %>% 
  knitr::kable()
```
