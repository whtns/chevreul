# Overview

SeuratTools is a shiny app for exploratory data analysis of single cell sequencing data after processing via [Seurat](https://satijalab.org/seurat/). It is intended to allow scientists without extensive coding experience to interact with their sequencing data. 

seuratTools combines several methods of single cell analysis without the need for extensive coding, including:

1. The ability to integrate datasets as implemented in [seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) within a gui environment
2. The ability to perform RNA velocity analysis 
3. Pseudotime analysis via monocle3

There are currently several shiny apps for single cell EDA via iSEE, scclustviz, and Cerebro. SeuratTools is oriented toward analysis of full-length sequencing data and is designed for labs working with sequencing data composed of multiple batches. 
