# Overview

SeuratTools is a shiny app for exploratory data analysis of single cell sequencing data after processing via [Seurat](https://satijalab.org/seurat/). It is designed to allow scientists without extensive coding experience to interact with their single cell sequencing data. 

seuratTools brings together several tools for single cell analysis, including:

1. Batch integration as implemented in [seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)

There are currently several shiny apps for single cell exploratory data analysis including [iSEE](https://bioconductor.org/packages/release/bioc/html/iSEE.html), [scclustviz](https://baderlab.github.io/scClustViz/), and [Cerebro](https://github.com/romanhaa/Cerebro). SeuratTools is distinct from these applications in that it is oriented toward analysis of full-length sequencing data via Smart-seq or similar technologies. In addition, seuratTools is designed for labs working with single cell data composed of multiple batches for whome integration is a primary concern.

seuratTools is built as an R package with functions for processing a single cell dataset from a summarized count/uMI matrix. The core of the tool are two shiny apps designed for 1) rapid visualization of a single dataset or 2) on-demand comparison of multiple datasets. The latter application provides an expanded set of functions for batch integration and project management. 

Documentation of all features is described below and in handy help icons throughout the app itself. 
